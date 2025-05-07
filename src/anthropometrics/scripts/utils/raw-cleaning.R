#!/usr/bin/env Rscript

##
#
# This script takes phenos from the raw database, runs initial cleaning, id mapping, and puts everything in a new database. 
#
##

print(paste0(Sys.time(), " - Load and clean raw QC values"))


## Command line input

# Command line arguments
args <- commandArgs(TRUE)

preg_id_linkage_raw_table_path <- args[1]
child_id_linkage_raw_table_path <- args[2]
mother_id_linkage_raw_table_path <- args[3]
father_id_linkage_raw_table_path <- args[4]
genomics_fam_file_path <- args[5]
unrelated_children_id_path <- args[6]
variables_mapping_table <- args[7]
ids_mapping_table <- args[8]
raw_phenotypes_tables_folder <- args[9]
tablesFolder <- args[10]
qcFolder <- args[11]
project_number <- args[12]


##
#
# Debug Marc - do not uncomment
# args to run standalone
# 
# preg_id_linkage_raw_table_path <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/raw/linkage/PDB315_SV_INFO_V12_20250131.gz"
# child_id_linkage_raw_table_path <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/raw/linkage/PDB315_MoBaGeneticsTot_Child_20221228.gz"
# mother_id_linkage_raw_table_path <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/raw/linkage/PDB315_MoBaGeneticsTot_Mother_20221228.gz"
# father_id_linkage_raw_table_path <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/raw/linkage/PDB315_MoBaGeneticsTot_Father_20221228.gz"
# genomics_fam_file_path <- "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam"
# unrelated_children_id_path <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/id/children_id_unrelated"
# variables_mapping_table <- "src/anthropometrics/scripts/resources/variable_mapping"
# ids_mapping_table <- "src/anthropometrics/scripts/resources/identifiers"
# raw_phenotypes_tables_folder <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/raw/phenotypes"
# kostUngdom_raw_table_path <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/raw/phenotypes/PDB315_Kosthold_ungdom_v12.gz"
# ungdomsskjema_barn_raw_table_path <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/raw/phenotypes/PDB315_Ungdomsskjema_Barn_v12_standard.gz"
# tablesFolder <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11"
# qcFolder <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_25-02-11/qc"
# project_number <- 315
# #
##

# Libraries

library(stringr)
library(janitor)
library(dplyr)
library(purrr)
library(glue)


# Housekeeping

if (!dir.exists(tablesFolder)) {
  
  dir.create(
    path = tablesFolder,
    showWarnings = T,
    recursive = T
  )
  
}

if (!dir.exists(qcFolder)) {
  
  dir.create(
    path = qcFolder,
    showWarnings = T,
    recursive = T
  )
  
}


## Parameters

# The variable mapping
source("src/anthropometrics/scripts/utils/variables_mapping.R")

identifiers_mapping <- read.table(
  file = ids_mapping_table,
  header = T
) %>% 
  mutate(
    moba_identifier = str_replace_all(moba_identifier, "\\{project_number\\}", as.character(project_number))
  )
id_moba_to_project <- identifiers_mapping$project_identifier
names(id_moba_to_project) <- identifiers_mapping$moba_identifier


variable_mapping <- read.table(
  file = variables_mapping_table,
  header = T
)
variable_moba_to_project <- variable_mapping$project_variable
names(variable_moba_to_project) <- variable_mapping$moba_variable

new_variables <- list()


# Cache for the standard values of the normal distribution
anchorUp <- pnorm(1)
anchorCenter <- pnorm(0)
anchorDown <- pnorm(-1)
margin99 <- qnorm(0.99)

# Pregnancy duration cut-offs (inclusive)
pregnancy_duration_min <- 154
pregnancy_duration_max <- 308
pregnancy_term_min <- 259
pregnancy_term_max <- 294

#Batch order to prioritize genotypes when samples have been measured in duplicates across batches
batchOrder <- c("NORMENT_JAN21", "NORMENT_MAI21", "NORMENT_SEP20_R996R1029", "NORMENT-JAN20", "NORMENT-FEB18", "NORMENT-JUN15", "NORMENT-MAY16", "NORMENT-JAN15", "ROTTERDAM1", "ROTTERDAM2", "HARVEST", "TED", "PDB1382_R875_R876")


## Functions


#' Counts the number of non-missing values in each column of the given data frame.
#' 
#' @param df the data frame to inspect
#' 
#' @return a data frame with the number of missing values for each column
getNValues <- function(df) {
  
  colNames <- names(df)
  nValues <- numeric(length(colNames))
  
  for (i in 1:length(colNames)) {
    
    colName <- colNames[i]
    nValue <- sum(!is.na(df[[colName]]))
    nValues[i] <- nValue
    
  }
  
  return(
    data.frame(
      label = colNames,
      nValues,
      stringsAsFactors = F
    )
  )
}

# Project specific variables

preg_id_column <- paste0("preg_id_", project_number)
mother_id_column <- paste0("m_id_", project_number)
father_id_column <- paste0("f_id_", project_number)


# Load identifiers

print(paste0(Sys.time(), "    Loading identifiers"))

trioIdDF <- read.table(
  file = preg_id_linkage_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  select(
    preg_id = !!sym(preg_id_column),
    mother_id = !!sym(mother_id_column),
    father_id = !!sym(father_id_column)
  )

childIdDF <- read.table(
  file = child_id_linkage_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  select(
    rank_siblings = barn_nr,
    preg_id = !!sym(preg_id_column),
    sentrix_id,
    batch,
    sampletype
  ) %>% 
  mutate(
    child_id = paste0(preg_id, "_", rank_siblings)
  ) %>% 
  filter(
    preg_id %in% trioIdDF$preg_id
  )

if (sum(is.na(childIdDF$sentrix_id)) > 0) {
  
  stop("Missing child sentrix id")
  
}

motherIdDF <- read.table(
  file = mother_id_linkage_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  select(
    mother_id = !!sym(mother_id_column),
    sentrix_id,
    batch,
    sampletype
  ) %>% 
  filter(
    mother_id %in% trioIdDF$mother_id
  )

if (sum(is.na(motherIdDF$sentrix_id)) > 0) {
  
  stop("Missing mother sentrix id")
  
}


fatherIdDF <- read.table(
  file = father_id_linkage_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  select(
    father_id = !!sym(father_id_column),
    sentrix_id,
    batch,
    sampletype
  ) %>% 
  filter(
    father_id %in% trioIdDF$father_id
  )

if (sum(is.na(fatherIdDF$sentrix_id)) > 0) {
  
  stop("Missing father sentrix id")
  
}

idDFs <- list(childIdDF, motherIdDF, fatherIdDF)
idDFLabels <- c("children", "mothers", "fathers")
idColumns <- c("child_id", "mother_id", "father_id")


# Check that all batches are supported

for (i in 1:length(idDFs)) {
  
  idDF <- idDFs[[i]]
  label <- idDFLabels[i]
  
  if (sum(!idDF$batch %in% batchOrder) > 0) {
    
    missingBatches <- unique(idDF$batch[!idDF$batch %in% batchOrder])
    
    stop(paste0("Batch not supported for ", label, ": ", paste(missingBatches, collapse = ", ")))
    
  }
}


# Remove duplicates

for (i in 1:length(idDFs)) {
  
  idDF <- idDFs[[i]]
  label <- idDFLabels[i]
  idColumn <- idColumns[i]
  
  print(glue("{Sys.time()} - Removing duplicates in {label}"))
  
  nBefore <- nrow(idDF)
  
  sampleOccurrenceDF <- as.data.frame(
    table(idDF[[idColumn]]),
    stringsAsFactors = F
  ) %>% 
    clean_names() %>% 
    filter(
      freq > 1
    )
  
  for (duplicate in sampleOccurrenceDF$var1) {
    
    duplicateIds <- idDF %>% 
      filter(
        !!sym(idColumn) == duplicate
      ) %>% 
      mutate(
        batchFactor = factor(batch, levels = batchOrder)
      ) %>% 
      arrange(
        batchFactor
      )
    
    idDF <- idDF %>% 
      filter(
        !!sym(idColumn) != duplicate | sentrix_id == duplicateIds$sentrix_id[1]
      )
    
  }
  
  if (length(unique(idDF[[idColumn]])) != length(unique(idDF$sentrix_id))) {
    
    stop("Duplicate sentrix_ids found")
    
  }
  
  if (nrow(idDF) != length(unique(idDF[[idColumn]]))) {
    
    stop("Duplicate ids found")
    
  }
  
  nAfter <- nrow(idDF)
  
  print(paste0("Duplicates removed for ", label, " : ", nBefore - nAfter, " (", round(100 * (nBefore - nAfter)/nBefore, 1), "%)."))
  
  idDFs[[i]] <- idDF
  
}

childIdDF <- idDFs[[1]]
motherIdDF <- idDFs[[2]]
fatherIdDF <- idDFs[[3]]


# Load the fam file

print(paste0(Sys.time(), "    Loading familial relationship from fam file"))

famDF <- read.table(
  file = genomics_fam_file_path,
  sep = " ",
  header = F,
  stringsAsFactors = F
)

names(famDF) <- c("family_id", "sentrix_id", "father_sentrix_id", "mother_sentrix_id", "genetic_sex", "pheno")

famDF <- famDF %>% 
  select(
    sentrix_id, father_sentrix_id, mother_sentrix_id, genetic_sex
  ) %>% 
  filter(
    !is.na(sentrix_id)
  )

if (nrow(famDF) != length(unique(famDF$sentrix_id))) {
  
  stop("Duplicate Sentrix Id in the fam file.")
  
}

famDF$mother_sentrix_id[! famDF$mother_sentrix_id %in% motherIdDF$sentrix_id] <- NA
famDF$father_sentrix_id[! famDF$father_sentrix_id %in% fatherIdDF$sentrix_id] <- NA

famDF <- famDF %>% 
  left_join(
    motherIdDF %>% 
      select(
        mother_sentrix_id = sentrix_id,
        mother_id
      ),
    by = "mother_sentrix_id"
  ) %>% 
  left_join(
    fatherIdDF %>% 
      select(
        father_sentrix_id = sentrix_id,
        father_id
      ),
    by = "father_sentrix_id"
  )


# Load the values from the raw phenotype tables, keep only the variables mapped and identifiers

print(paste0(Sys.time(), " - Loading raw phenotypes"))

raw_tables <- list()

for (table_name in unique(variable_mapping$moba_table)) {
  
  print(paste0(Sys.time(), "     ", table_name))
  
  file_name <- paste0("PDB", project_number, "_", table_name, ".gz")
  table_file <- file.path(raw_phenotypes_tables_folder, file_name)
  
  if (!file.exists(table_file)) {
    
    stop(paste0("File `", table_file, "` corresponding to table `", table_name, "` not found."))
    
  }
  
  raw_table <- read.table(
    file = table_file,
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )
  
  variables_to_extract <- variable_mapping$moba_variable[variable_mapping$moba_table == table_name]
  
  missing_variables <- variables_to_extract[!variables_to_extract %in% names(raw_table)]
  
  if (length(missing_variables) > 0) {
    
    print(glue("The following variables were not found in `{table_name}`:"))
    
    for (missing_variable in missing_variables) {
      
      print(missing_variable)
      
    }
    
    # stop("Error: Variables mismatch")
    
    variables_to_extract <- variables_to_extract[variables_to_extract %in% names(raw_table)]
    
  }
  
  ids_to_extract <- identifiers_mapping$moba_identifier[identifiers_mapping$moba_table == table_name]
  
  missing_ids <- ids_to_extract[!ids_to_extract %in% names(raw_table)]
  
  if (length(missing_ids) > 0) {
    
    print(glue("The following ids were not found in `{table_name}`:"))
    
    for (missing_id in missing_ids) {
      
      print(missing_id)
      
    }
    
    # stop("Error: Missing id")
    
    ids_to_extract <- ids_to_extract[ids_to_extract %in% names(raw_table)]
    
  }
  
  names(variables_to_extract) <- variable_moba_to_project[variables_to_extract]
  names(ids_to_extract) <- id_moba_to_project[ids_to_extract]
  columns_to_extract <- c(variables_to_extract, ids_to_extract)
  
  raw_table <- raw_table %>% 
    select(
      all_of(
        c(ids_to_extract, variables_to_extract)
      )
    )
  
  raw_tables[[table_name]] <- raw_table
  
}


# Make a table of identifiers

print(paste0(Sys.time(), " - Making a table of identifiers"))

rawPheno <- childIdDF %>% 
  filter(
    sentrix_id %in% famDF$sentrix_id
  ) %>% 
  select(
    preg_id, rank_siblings, child_id, child_sentrix_id = sentrix_id, child_batch = batch
  ) %>% 
  left_join(
    famDF %>% 
      rename(child_sentrix_id = sentrix_id),
    by = "child_sentrix_id"
  ) %>% 
  left_join(
    motherIdDF %>% 
      select(
        mother_sentrix_id = sentrix_id,
        mother_batch = batch
      ),
    by = "mother_sentrix_id"
  ) %>% 
  left_join(
    fatherIdDF %>% 
      select(
        father_sentrix_id = sentrix_id,
        father_batch = batch
      ),
    by = "father_sentrix_id"
  )


# Merge with the pheno tables

for (table_name in names(raw_tables)) {
  
  raw_table <- raw_tables[[table_name]]
  ids_to_extract <- identifiers_mapping$project_identifier[identifiers_mapping$moba_table == table_name]
  
  if (length(ids_to_extract) == 0) {
    
    stop(glue("No identifiers found for table `{table_name}`."))
    
  }
  
  columns_to_select <- c(ids_to_extract, names(raw_table)[!names(raw_table) %in% names(rawPheno)])
  
  missing_ids <- ids_to_extract[!ids_to_extract %in% names(raw_table)]
  
  if (length(missing_ids) > 0) {
    
    stop(glue("Missing identifiers `{missing_ids}` in raw_table `{table_name}`."))
    
  }
  
  raw_table <- raw_table %>% 
    select(
      all_of(
        columns_to_select
      )
    )
  
  missing_ids <- ids_to_extract[!ids_to_extract %in% names(rawPheno)]
  
  if (length(missing_ids) > 0) {
    
    stop(glue("Missing identifiers `{missing_ids}` in merged table rawPheno."))
    
  }
  
  rawPheno <- rawPheno %>% 
    left_join(
      raw_table,
      by = ids_to_extract
    )
  
}


# Sanity check

if (nrow(rawPheno) != length(unique(rawPheno$child_id))) {
  
  stop("Duplictes were introduced when merging data frames.")
  
}

print(glue("Phenotypes loaded:"))
print(glue("- Children with phenotypes: {nrow(rawPheno)}"))
print(glue("- Children genotyped: {sum(!is.na(rawPheno$child_sentrix_id))}"))
print(glue("- Mothers genotyped linked to a child: {length(unique(rawPheno$mother_sentrix_id))}"))
print(glue("- Fathers genotyped linked to a child: {length(unique(rawPheno$father_sentrix_id))}"))


# Combination of variables

rawPheno$diabetes_3y <- NA

if ("diabetes_no_3y" %in% names(rawPheno) && "diabetes_yes_3y" %in% names(rawPheno) && "diabetes_previous_3y" %in% names(rawPheno) ) {
  
  rawPheno$diabetes_3y[rawPheno$diabetes_no_3y == 1] <- 0
  rawPheno$diabetes_3y[rawPheno$diabetes_yes_3y == 1 | rawPheno$diabetes_yes_3y == 1] <- 1
  
}

new_variables[["child_health"]] <- c(new_variables[["child_health"]], "diabetes_3y")


rawPheno$underweight_3y <- NA

if ("gained_too_little_weight_no_3y" %in% names(rawPheno) && "gained_too_little_weight_yes_3y" %in% names(rawPheno) && "gained_too_little_weight_previous_3y" %in% names(rawPheno)) {
  
  rawPheno$underweight_3y[rawPheno$gained_too_little_weight_no_3y == 1] <- 0
  rawPheno$underweight_3y[rawPheno$gained_too_little_weight_yes_3y == 1 | rawPheno$gained_too_little_weight_previous_3y == 1] <- 1
  
}

new_variables[["child_health"]] <- c(new_variables[["child_health"]], "underweight_3y")


rawPheno$overweight_3y <- NA

if ("gained_too_much_weight_no_3y" %in% names(rawPheno) && "gained_too_much_weight_yes_3y" %in% names(rawPheno) && "gained_too_much_weight_previous_3y" %in% names(rawPheno)) {
  
  rawPheno$overweight_3y[rawPheno$gained_too_much_weight_no_3y == 1] <- 0
  rawPheno$overweight_3y[rawPheno$gained_too_much_weight_yes_3y == 1 | rawPheno$gained_too_much_weight_previous_3y == 1] <- 1
  
}

new_variables[["child_health"]] <- c(new_variables[["child_health"]], "overweight_3y")


rawPheno$length_7y <- ifelse(!is.na(rawPheno$height_7y_cm), rawPheno$height_7y_cm, rawPheno$height_7y_m * 100)

new_variables[["child_anthropometrics_raw"]] <- c(new_variables[["child_anthropometrics_raw"]], "length_7y")


hospitalization_columns <- c("hospitalized_prolonged_nausea_vomiting", "hospitalized_bleeding", "hospitalized_amniotic_fluid_leakage", "hospitalized_threatening_preterm_labour", "hospitalized_high_blood_pressure", "hospitalized_pre_eclampsia", "hospitalized_other")
suffixes <- c("_0_4w", "_5_8w", "_9_12w", "_13_16w", "_17_20w", "_21_24w", "_25_28w", "_after_29w")

for (hospitalization_column in hospitalization_columns) {
  
  for (suffix in suffixes) {
    
    week_column <- paste0(hospitalization_column, suffix)
    
    rawPheno[[hospitalization_column]] <- ifelse(!is.na(rawPheno[[week_column]]) & rawPheno[[week_column]] == 1, 1, rawPheno[[hospitalization_column]])
    
  }
  
  rawPheno$hospitalized_30w <- ifelse(!is.na(rawPheno[[hospitalization_column]]) & rawPheno[[hospitalization_column]] == 1, "Yes", rawPheno$hospitalized_30w)
  
}

new_variables[["pregnancy"]] <- c(new_variables[["pregnancy"]], hospitalization_columns)


# Correct units for columns to be merged

print(paste(Sys.time(), " Correcting units"))

lineIndexes <- !is.na(rawPheno$weight_15_18m_1) & rawPheno$weight_15_18m_1 < 50 # Wrong unit
rawPheno$weight_15_18m_1[lineIndexes] <- rawPheno$weight_15_18m_1[lineIndexes] * 1000
lineIndexes <- !is.na(rawPheno$weight_15_18m_1) & rawPheno$weight_15_18m_1 < 500 # Wrong decimal separator
rawPheno$weight_15_18m_1[lineIndexes] <- rawPheno$weight_15_18m_1[lineIndexes] * 100
lineIndexes <- !is.na(rawPheno$weight_15_18m_1) & rawPheno$weight_15_18m_1 < 5000 # Wrong decimal separator
rawPheno$weight_15_18m_1[lineIndexes] <- rawPheno$weight_15_18m_1[lineIndexes] * 10

lineIndexes <- !is.na(rawPheno$weight_15_18m_2) & rawPheno$weight_15_18m_2 > 40 # missing decimal separator
rawPheno$weight_15_18m_2[lineIndexes] <- rawPheno$weight_15_18m_2[lineIndexes]/10
lineIndexes <- !is.na(rawPheno$weight_15_18m_2) & rawPheno$weight_15_18m_2 < 2 # Wrong decimal separator
rawPheno$weight_15_18m_2[lineIndexes] <- rawPheno$weight_15_18m_2[lineIndexes] * 10

rawPheno$weight_15_18m_2 <- rawPheno$weight_15_18m_2 * 1000 # Convert to grams

rawPheno <- rawPheno %>% 
  mutate(
    umbilical_cord_length = ifelse(umbilical_cord_length >= 250, umbilical_cord_length / 10, umbilical_cord_length), # Wrong decimal separator
    umbilical_cord_length = ifelse(umbilical_cord_length <= 10, umbilical_cord_length * 10, umbilical_cord_length) # Wrong decimal separator
  )


# Convert to number

print(paste(Sys.time(), " Converting string input to number"))

rawPheno$age_birth <- 0

rawPheno$weight_14c <- str_remove_all(rawPheno$weight_14c, " KG")
rawPheno$weight_14c <- str_remove_all(rawPheno$weight_14c, "MINDRE ENN ")
rawPheno$weight_14c <- str_remove_all(rawPheno$weight_14c, "MER ENN ")
rawPheno$weight_14c <- as.numeric(rawPheno$weight_14c)

rawPheno$height_14c <- str_remove_all(rawPheno$height_14c, " CM")
rawPheno$height_14c <- str_remove_all(rawPheno$height_14c, "LAVERE ENN ")
rawPheno$height_14c <- str_remove_all(rawPheno$height_14c, "HØYERE ENN ")
rawPheno$height_14c <- as.numeric(rawPheno$height_14c)

rawPheno$weight_mother_14m <- str_remove_all(rawPheno$weight_mother_14m, " KG")
rawPheno$weight_mother_14m <- str_remove_all(rawPheno$weight_mother_14m, "MINDRE ENN ")
rawPheno$weight_mother_14m <- str_remove_all(rawPheno$weight_mother_14m, "MER ENN ")
rawPheno$weight_mother_14m <- as.numeric(rawPheno$weight_mother_14m)

rawPheno$height_mother_14m <- str_remove_all(rawPheno$height_mother_14m, " CM")
rawPheno$height_mother_14m <- str_remove_all(rawPheno$height_mother_14m, "LAVERE ENN ")
rawPheno$height_mother_14m <- str_remove_all(rawPheno$height_mother_14m, "HØYERE ENN ")
rawPheno$height_mother_14m <- as.numeric(rawPheno$height_mother_14m)

new_variables[["child_anthropometrics_raw"]] <- c(new_variables[["child_anthropometrics_raw"]], c("age_birth"))


# Merge columns

print(paste(Sys.time(), " Combining phenotypes"))

if ("mother_height_self" %in% names(rawPheno) & "mother_height" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      mother_height = ifelse(is.na(mother_height), mother_height_self, mother_height)
    )
  
}

if ("mother_weight_end" %in% names(rawPheno) & "mother_weight_end_self" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      mother_weight_end = ifelse(is.na(mother_weight_end), mother_weight_end_self, mother_weight_end)
    )
  
}

if ("father_height" %in% names(rawPheno) & "father_height_self" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      father_height = ifelse(is.na(father_height), father_height_self, father_height)
    )
  
}

if ("father_weight" %in% names(rawPheno) & "father_weight_self" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      father_weight = ifelse(is.na(father_weight), father_weight_self, father_weight)
    )
  
}

if ("formula_colett_0m" %in% names(rawPheno) & "formula_colett_omega3_0m" %in% names(rawPheno) & "formula_nan_0m" %in% names(rawPheno) & "formula_nan_ha1_0m" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      formula_0m = ifelse(!is.na(formula_colett_0m) | !is.na(formula_colett_omega3_0m) | !is.na(formula_nan_0m) | !is.na(formula_nan_ha1_0m), 1, NA)
    )
  
}

if ("formula_colett_1m" %in% names(rawPheno) & "formula_colett_omega3_1m" %in% names(rawPheno) & "formula_nan_1m" %in% names(rawPheno) & "formula_nan_ha1_1m" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      formula_1m = ifelse(!is.na(formula_colett_1m) | !is.na(formula_colett_omega3_1m) | !is.na(formula_nan_1m) | !is.na(formula_nan_ha1_1m), 1, NA)
    )
  
}

if ("formula_colett_2m" %in% names(rawPheno) & "formula_colett_omega3_2m" %in% names(rawPheno) & "formula_nan_2m" %in% names(rawPheno) & "formula_nan_ha1_2m" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      formula_2m = ifelse(!is.na(formula_colett_2m) | !is.na(formula_colett_omega3_2m) | !is.na(formula_nan_2m) | !is.na(formula_nan_ha1_2m), 1, NA)
    )
  
}

if ("formula_colett_3m" %in% names(rawPheno) & "formula_colett_omega3_3m" %in% names(rawPheno) & "formula_nan_3m" %in% names(rawPheno) & "formula_nan_ha1_3m" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      formula_3m = ifelse(!is.na(formula_colett_3m) | !is.na(formula_colett_omega3_3m) | !is.na(formula_nan_3m) | !is.na(formula_nan_ha1_3m), 1, NA)
    )
  
}

if ("formula_colett_4m" %in% names(rawPheno) & "formula_colett_omega3_4m" %in% names(rawPheno) & "formula_nan_4m" %in% names(rawPheno) & "formula_nan_ha1_4m" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      formula_4m = ifelse(!is.na(formula_colett_4m) | !is.na(formula_colett_omega3_4m) | !is.na(formula_nan_4m) | !is.na(formula_nan_ha1_4m), 1, NA)
    )
  
}

if ("formula_colett_5m" %in% names(rawPheno) & "formula_colett_omega3_5m" %in% names(rawPheno) & "formula_nan_5m" %in% names(rawPheno) & "formula_nan_ha1_5m" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      formula_5m = ifelse(!is.na(formula_colett_5m) | !is.na(formula_colett_omega3_5m) | !is.na(formula_nan_5m) | !is.na(formula_nan_ha1_5m), 1, NA)
    )
  
}

if ("formula_colett_6m" %in% names(rawPheno) & "formula_colett_omega3_6m" %in% names(rawPheno) & "formula_nan_6m" %in% names(rawPheno) & "formula_nan_ha1_6m" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      formula_6m = ifelse(!is.na(formula_colett_6m) | !is.na(formula_colett_omega3_6m) | !is.na(formula_nan_6m) | !is.na(formula_nan_ha1_6m), 1, NA)
    )
  
}

rawPheno <- rawPheno %>%
  mutate(
    age_16m = ifelse(is.na(age_15_18m_1), age_15_18m_1, age_15_18m_2),
    length_16m = ifelse(is.na(length_15_18m_1), length_15_18m_1, length_15_18m_2),
    weight_16m = ifelse(is.na(weight_15_18m_1), weight_15_18m_1, weight_15_18m_2)
  )

new_variables[["child_anthropometrics_raw"]] <- c(new_variables[["child_anthropometrics_raw"]], c("age_16m", "length_16m", "weight_16m"))


# Get a list of unrelated kids

print(paste(Sys.time(), " Loading identifiers for unrelated children"))

unrelatedDF <- read.table(
  file = unrelated_children_id_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

rawPheno <- rawPheno %>%
  mutate(
    unrelated_children = ifelse(child_sentrix_id %in% unrelatedDF$sentrix_id, 1, 0)
  )

print(glue("Unrelated children:"))
print(glue("- Children genotyped: {sum(!is.na(rawPheno$child_sentrix_id) & rawPheno$unrelated_children == 1)}"))
print(glue("- Mothers genotyped linked to a child: {length(unique(rawPheno$mother_sentrix_id[rawPheno$unrelated_children == 1]))}"))
print(glue("- Fathers genotyped linked to a child: {length(unique(rawPheno$father_sentrix_id[rawPheno$unrelated_children == 1]))}"))


# Set sex columns according to genotypes and registry, and make a consensus

print(paste(Sys.time(), " Setting consensus sex between genotypes and registry"))

rawPheno <- rawPheno %>%
  rename(
    registry_sex = sex
  ) %>% 
  mutate(
    sex = 0,
    sex = ifelse(registry_sex == "Male", 1, sex),
    sex = ifelse(registry_sex == "Female", 2, sex),
    sex = ifelse(!is.na(genetic_sex) & genetic_sex != 0, genetic_sex, sex)
  )

new_variables[["ids"]] <- c(new_variables[["ids"]], c("registry_sex", "genetic_sex"))


# Check that all longitudinal columns are available

missing_column <- age_columns[!age_columns %in% names(rawPheno)]

if (length(missing_column) > 0) {
  
  stop(glue("Missing age column: {missing_column}"))
  
}

missing_column <- weight_columns[!weight_columns %in% names(rawPheno)]

if (length(missing_column) > 0) {
  
  stop(glue("Missing weight column: {missing_column}"))
  
}

missing_column <- length_columns[!length_columns %in% names(rawPheno)]

if (length(missing_column) > 0) {
  
  stop(glue("Missing length/height column: {missing_column}"))
  
}

missing_column <- head_circumference_columns[!head_circumference_columns %in% names(rawPheno)]

if (length(missing_column) > 0) {
  
  stop(glue("Missing height circumference column: {missing_column}"))
  
}


# Pregnancy duration

print(paste(Sys.time(), " Pregnancy duration"))

rawPheno <- rawPheno %>% 
  mutate(
    pregnancy_duration_ultrasound = ifelse(!is.na(pregnancy_duration_ultrasound) & pregnancy_duration_ultrasound <= pregnancy_duration_max & pregnancy_duration_ultrasound >= pregnancy_duration_min, pregnancy_duration_ultrasound, NA),
    pregnancy_duration_mens = ifelse(!is.na(pregnancy_duration_mens) & pregnancy_duration_mens <= pregnancy_duration_max & pregnancy_duration_mens >= pregnancy_duration_min, pregnancy_duration_mens, NA),
    pregnancy_duration = ifelse(!is.na(pregnancy_duration_ultrasound), pregnancy_duration_ultrasound, pregnancy_duration_mens),
    pregnancy_duration_term = ifelse(!is.na(pregnancy_duration) & pregnancy_duration >= pregnancy_term_min & pregnancy_duration <= pregnancy_term_max, 1, 0),
    pregnancy_duration_preterm = ifelse(!is.na(pregnancy_duration) & pregnancy_duration < pregnancy_term_min, 1, 0)
  )

term_pregnancies <- sum(!is.na(rawPheno$child_sentrix_id) & rawPheno$pregnancy_duration_term == 1)
n_genotyped <- sum(!is.na(rawPheno$child_sentrix_id))

new_variables[["pregnancy"]] <- c(new_variables[["pregnancy"]], c("pregnancy_duration_preterm", "pregnancy_duration_term"))

print(glue("{term_pregnancies} genotyped children with delivery at term ({round(term_pregnancies / n_genotyped * 100)} %)"))


# Set zeros to NA for longitudinal values

print(paste(Sys.time(), " Formatting"))

for (column in c(weight_columns, length_columns, head_circumference_columns)) {
  
  rawPheno[!is.na(rawPheno[[column]]) & rawPheno[[column]] == 0, column] <- NA
  
}


# Convert months to days, g to kg

rawPheno <- rawPheno %>% 
  mutate(
    
    age_5y = age_5y / 12 * 365.25,
    age_7y = age_7y / 12 * 365.25,
    age_8y = age_8y / 12 * 365.25,
    age_answering_q_14c = age_answering_q_14c / 12 * 365.25,
    
    weight_birth = weight_birth / 1000,
    weight_6w = weight_6w / 1000,
    weight_3m = weight_3m / 1000,
    weight_6m = weight_6m / 1000,
    weight_8m = weight_8m / 1000,
    weight_1y = weight_1y / 1000,
    weight_16m = weight_16m / 1000
    
  )


# Keep track of the number of values for the different phenotypes

number_values_tables <- list()

print(paste(Sys.time(), " Storing the number of values"))

n_list <- list()

for (column in c(weight_columns, length_columns, head_circumference_columns)) {
  
  n_df <- data.frame(
    phenotype = column,
    n_values = sum(!is.na(rawPheno[[column]])),
    n_genotyped = sum(!is.na(rawPheno[[column]]) & !is.na(rawPheno$child_sentrix_id)),
    n_genotyped_unrelated = sum(!is.na(rawPheno[[column]]) & !is.na(rawPheno$child_sentrix_id) & rawPheno$unrelated == 1),
    stringsAsFactors = F
  )
  
  n_list[[length(n_list) + 1]] <- n_df
  
}

number_values_tables$number_values <- do.call("rbind", n_list)


# Exclusion of outliers using five SD

print(paste(Sys.time(), " Exclusion of extreme outliers"))

values <- rawPheno

n_list <- list()

for (column in c(weight_columns, length_columns, head_circumference_columns)) {
  
  mean_value <- mean(values[[column]][!is.na(values[[column]]) & !is.na(values$child_sentrix_id)])
  sd_value <- sd(values[[column]][!is.na(values[[column]]) & !is.na(values$child_sentrix_id)])
  
  toExclude <- !is.na(values[[column]]) & (values[[column]] < mean_value - 5 * sd_value | values[[column]] > mean_value + 5 * sd_value)
  values[[column]][toExclude] <- NA
  
  n_df <- data.frame(
    phenotype = column,
    n_outliers = sum(toExclude),
    n_outliers_genotyped = sum(toExclude & !is.na(values$child_sentrix_id)),
    n_outliers_genotyped_unrelated = sum(toExclude & !is.na(values$child_sentrix_id) & values$unrelated == 1),
    stringsAsFactors = F
  )
  
  n_list[[length(n_list) + 1]] <- n_df
  
}

number_values_tables$number_outliers <- do.call("rbind", n_list)


# Impute late breast feeding from earlier values

breastFeedingColumns <- c("breastmilk_first_week", "breastmilk_0m", "breastmilk_1m", "breastmilk_2m", "breastmilk_3m", "breastmilk_4m", "breastmilk_5m", "breastmilk_6m",
                          "breastmilk_6_8m", "breastmilk_9_11m", "breastmilk_12_14m", "breastmilk_15_18m")

n_list <- list()

for (column in breastFeedingColumns) {
  
  if (column %in% names(values)) {
    
    n_df <- data.frame(
      phenotype = column,
      n_breast_feeding = sum(!is.na(values[[column]])),
      n_breast_feeding_genotyped = sum(!is.na(values[[column]]) & !is.na(values$child_sentrix_id)),
      n_breast_feeding_genotyped_unrelated = sum(!is.na(values[[column]]) & !is.na(values$child_sentrix_id) & values$unrelated == 1),
      stringsAsFactors = F
    )
  }
}

number_values_tables$breast_feeding_raw <- do.call("rbind", n_list)

for (i in 2:length(breastFeedingColumns)) {
  
  columnI <- breastFeedingColumns[i]
  
  if (columnI %in% names(values)) {
    
    feeding <- !is.na(values[[columnI]]) & values[[columnI]] == 1
    
    for (j in 1:(i - 1)) {
      
      columnJ <- breastFeedingColumns[j]
      
      if (columnJ %in% names(values)) {
        
        values[feeding, columnJ] <- 1
        
      }
    }
  }
}

n_list <- list()

for (column in breastFeedingColumns) {
  
  if (column %in% names(values)) {
    
    n_df <- data.frame(
      phenotype = column,
      n_breast_feeding = sum(!is.na(values[[column]])),
      n_breast_feeding_genotyped = sum(!is.na(values[[column]]) & !is.na(values$child_sentrix_id)),
      n_breast_feeding_genotyped_unrelated = sum(!is.na(values[[column]]) & !is.na(values$child_sentrix_id) & values$unrelated == 1),
      stringsAsFactors = F
    )
  }
}

number_values_tables$breast_feeding_imputed <- do.call("rbind", n_list)


# Impute late diabetes from earlier values

diabetesColumns <- c("diabetes_3y", "diabetes_7y", "diabetes_8y")

n_list <- list()

for (column in diabetesColumns) {
  
  if (column %in% names(values)) {
    
    n_df <- data.frame(
      phenotype = column,
      n_diabetes = sum(!is.na(values[[column]]) & values[[column]] == 1),
      n_diabetes_genotyped = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$child_sentrix_id)),
      n_diabetes_unrelated = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$child_sentrix_id) & values$unrelated == 1),
      stringsAsFactors = F
    )
  }
}

number_values_tables$diabetes_raw <- do.call("rbind", n_list)

for (i in 2:length(diabetesColumns)) {
  
  columnI <- diabetesColumns[i]
  
  if (columnI %in% names(values)) {
    
    diabetes0 <- !is.na(values[[columnI]]) & values[[columnI]] == 0
    
    for (j in 1:(i - 1)) {
      
      columnJ <- diabetesColumns[j]
      
      if (columnJ %in% names(values)) {
        
        values[diabetes0, columnJ] <- 0
        
      }
    }
  }
}

for (i in 1:(length(diabetesColumns) - 1)) {
  
  columnI <- diabetesColumns[i]
  
  if (columnI %in% names(values)) {
    
    diabetes1 <- !is.na(values[[columnI]]) & values[[columnI]] == 1
    
    for (j in (i + 1):length(diabetesColumns)) {
      
      columnJ <- diabetesColumns[j]
      
      if (columnJ %in% names(values)) {
        
        values[diabetes1, columnJ] <- 1
        
      }
    }
  }
}

n_list <- list()

for (column in diabetesColumns) {
  
  if (column %in% names(values)) {
    
    n_df <- data.frame(
      phenotype = column,
      n_diabetes = sum(!is.na(values[[column]]) & values[[column]] == 1),
      n_diabetes_genotyped = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$child_sentrix_id)),
      n_diabetes_unrelated = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$child_sentrix_id) & values$unrelated == 1),
      stringsAsFactors = F
    )
  }
}

number_values_tables$diabetes <- do.call("rbind", n_list)


# Variables combination

print(paste(Sys.time(), " Creating variable combinations"))

values <- values %>%
  mutate(
    mother_median_height = as.numeric(
      pmap(
        list(mother_height, mother_height_3y, mother_height_5y, mother_height_8y),
        ~median(c(..1, ..2, ..3, ..4), na.rm = T)
      )
    )
  )

new_variables[["parent"]] <- c(new_variables[["parent"]], "mother_median_height")


# Exclude extreme outliers for aam

values$mother_age_at_menarche[values$mother_age_at_menarche < 9 | values$mother_age_at_menarche > 17] <- NA


# Check for duplicates

if (length(unique(values$child_id)) != nrow(values)) {
  
  stop("Duplicates introduced during cleaning.")
  
}


# Save the variables in different tables

print(paste(Sys.time(), " Saving the phenotypes to tables"))

missing_default_columns <- which(! default_columns %in% names(values))

if (length(missing_default_columns) > 0) {
  
  stop(psate0("Missing default column in phenotype table: ", paste(default_columns[missing_default_columns], sep = ", "), "."))
  
}

table_names <- unique(c(variable_mapping$project_table, names(new_variables)))

for (project_table_name in table_names) {
  
  if (!is.na(project_table_name)) {
    
    print(paste0(Sys.time(), "     ", project_table_name))
    
    variables <- variable_mapping %>% 
      filter(
        project_table == project_table_name
      ) %>% 
      pull(
        project_variable
      )
    
    if (project_table_name == "child_anthropometrics_raw") {
      
      print("- Variables from config")
      print(variables, sep = ", ")
      
      print("- New variables")
      print(new_variables[[project_table_name]], sep = ", ")
      
    }
    
    if (project_table_name %in% names(new_variables)) {
      
      variables <- c(variables, new_variables[[project_table_name]])
      
    }
    
    if (project_table_name == "child_anthropometrics_raw") {
      
      print("- Missing variables")
      print(variables[!variables %in% names(values)], sep = ", ")
      
    }
    
    variables <- unique(variables)
    variables <- variables[!variables %in% default_columns]
    variables <- c(default_columns, variables)
    variables <- variables[variables %in% names(values)]
    
    table_content <- values %>% 
      select(
        all_of(variables)
      )
    
    write.table(
      x = table_content,
      file = gzfile(file.path(tablesFolder, paste0(project_table_name, ".gz"))),
      row.names = F,
      col.names = T,
      sep = "\t",
      quote = F
    )
    
  }
}


# Save qc

print(paste(Sys.time(), "Saving number of samples."))

for (name in names(number_values_tables)) {
  
  n_df <- number_values_tables[[name]]
  
  write.table(
    x = n_df,
    file = file.path(qcFolder, paste0("n_", name)),
    row.names = F,
    col.names = T,
    sep = "\t",
    quote = F
  )
  
}


# Done

print(paste(Sys.time(), "Process complete."))
