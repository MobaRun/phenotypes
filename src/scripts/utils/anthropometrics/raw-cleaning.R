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

child_id_linkage_raw_table_path <- args[1]
mother_id_linkage_raw_table_path <- args[2]
father_id_linkage_raw_table_path <- args[3]
genomics_fam_file_path <- args[4]
unrelated_children_id_path <- args[5]
mfr_raw_table_path <- args[6]
q1m_raw_table_path <- args[7]
q1f_raw_table_path <- args[8]
q2_raw_table_path <- args[9]
q3_raw_table_path <- args[10]
q4_raw_table_path <- args[11]
q5_raw_table_path <- args[12]
q6_raw_table_path <- args[13]
q7_raw_table_path <- args[14]
q8_raw_table_path <- args[15]
q9_raw_table_path <- args[16]
kostUngdom_raw_table_path <- args[17]
qcFolder <- args[18]
project_number <- args[19] 


##
#
# Debug Marc - do not uncomment
# args to run standalone
# 
child_id_linkage_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/linkage/20220516_MoBaGeneticsTot_Child_PDB2824.gz"
mother_id_linkage_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/linkage/20220516_MoBaGeneticsTot_Mother_PDB2824.gz"
father_id_linkage_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/linkage/20220516_MoBaGeneticsTot_Father_PDB2824.gz"
genomics_fam_file_path <- "/mnt/archive/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam"
unrelated_children_id_path <- "/mnt/work/marc/unrelated_samples/children_id_unrelated"
mfr_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_MFR_541_v12.gz"
q1m_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema1_v12.gz"
q1f_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_SkjemaFar_v12.gz"
q2_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema2_v12.gz"
q3_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema3_v12.gz"
q4_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema4_6mnd_v12.gz"
q5_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema5_18mnd_v12.gz"
q6_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema6_3aar_v12.gz"
q7_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema5aar_v12.gz"
q8_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema7aar_v12.gz"
q9_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_Skjema8aar_v12.gz"
kostUngdom_raw_table_path <- "/mnt/work/marc/pheno_22-09-19/raw/moba_ques/PDB2824_SkjemaKostUngdom_v12.gz"
qcFolder <- "/mnt/work/marc/pheno_22-09-19/qc_tmp"
project_number <- 2824
#
##


# Libraries

libFolder <- "~/R/R_4.1"

library(stringr, lib = libFolder)
library(crayon, lib = libFolder)
library(dplyr, lib = libFolder)
library(janitor, lib = libFolder)
library(purrr, lib = libFolder)
library(glue, lib = libFolder)


# Housekeeping

if (!dir.exists(qcFolder)) {
  
  dir.create(
    path = qcFolder,
    showWarnings = T,
    recursive = T
  )
  
}


## Parameters

# The variable mapping
source("src/scripts/utils/anthropometrics/variables_mapping.R")

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

print(paste0(Sys.time(), "    Loading triads from fam file"))

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


# Load the values from the raw phenotype tables

print(paste0(Sys.time(), " - Loading raw phenotypes"))

mfr_raw_table <- read.table(
  file = mfr_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q1m_raw_table <- read.table(
  file = q1m_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q1f_raw_table <- read.table(
  file = q1f_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

# q2_raw_table <- read.table(
#   file = q2_raw_table_path,
#   header = T,
#   sep = "\t",
#   stringsAsFactors = F
# )

q3_raw_table <- read.table(
  file = q3_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q4_raw_table <- read.table(
  file = q4_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q5_raw_table <- read.table(
  file = q5_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q6_raw_table <- read.table(
  file = q6_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q7_raw_table <- read.table(
  file = q7_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q8_raw_table <- read.table(
  file = q8_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

q9_raw_table <- read.table(
  file = q9_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

kostUngdom_raw_table <- read.table(
  file = kostUngdom_raw_table_path,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)



# Extract variables

mfrVariablesMapping <- mfrVariablesMapping[mfrVariablesMapping %in% names(mfr_raw_table)]

mfr_table <- mfr_raw_table %>%
  select(
    all_of(
      mfrVariablesMapping
    )
  )

q1mVariablesMapping <- q1mVariablesMapping[q1mVariablesMapping %in% names(q1m_raw_table)]

q1m_table <- q1m_raw_table %>%
  select(
    all_of(
      q1mVariablesMapping
    )
  )

q1fVariablesMapping <- q1fVariablesMapping[q1fVariablesMapping %in% names(q1f_raw_table)]

q1f_table <- q1f_raw_table %>%
  select(
    all_of(
      q1fVariablesMapping
    )
  )

# q2VariablesMapping <- q2VariablesMapping[q2VariablesMapping %in% names(q2_raw_table)]
# 
# q2_table <- q2_raw_table %>%
#   select(
#     all_of(
#       q2VariablesMapping
#     )
#   )

q4VariablesMapping <- q4VariablesMapping[q4VariablesMapping %in% names(q4_raw_table)]

q4_table <- q4_raw_table %>%
  select(
    all_of(
      q4VariablesMapping
    )
  )

q5VariablesMapping <- q5VariablesMapping[q5VariablesMapping %in% names(q5_raw_table)]

q5_table <- q5_raw_table %>%
  select(
    all_of(
      q5VariablesMapping
    )
  )

q6VariablesMapping <- q6VariablesMapping[q6VariablesMapping %in% names(q6_raw_table)]

q6_table <- q6_raw_table %>%
  select(
    all_of(
      q6VariablesMapping
    )
  )

q7VariablesMapping <- q7VariablesMapping[q7VariablesMapping %in% names(q7_raw_table)]

q7_table <- q7_raw_table %>%
  select(
    all_of(
      q7VariablesMapping
    )
  )

q8VariablesMapping <- q8VariablesMapping[q8VariablesMapping %in% names(q8_raw_table)]

q8_table <- q8_raw_table %>%
  select(
    all_of(
      q8VariablesMapping
    )
  )

q9VariablesMapping <- q9VariablesMapping[q9VariablesMapping %in% names(q9_raw_table)]

q9_table <- q9_raw_table %>%
  select(
    all_of(
      q9VariablesMapping
    )
  )

kostUngdomVariablesMapping <- kostUngdomVariablesMapping[kostUngdomVariablesMapping %in% names(kostUngdom_raw_table)]

kostUngdom_table <- kostUngdom_raw_table %>%
  select(
    all_of(
      kostUngdomVariablesMapping
    )
  )


# Combination of variables from the same questionnaire

if ("breastmilk_freq_18m" %in% names(q5_table) && "breastmilk_freq_18m_1" %in% names(q5_table) && "breastmilk_freq_18m_2" %in% names(q5_table)) {
  
  q5_table$breastmilk_freq_18m <- as.character(q5_table$breastmilk_freq_18m_1)
  q5_table$breastmilk_freq_18m[is.na(q5_table$breastmilk_freq_18m)] <- as.character(q5_table$breastmilk_freq_18m_2)
  q5_table$formula_freq_18m <- as.character(q5_table$formula_freq_18m_1)
  q5_table$formula_freq_18m[is.na(q5_table$formula_freq_18m)] <- as.character(q5_table$formula_freq_18m_2)
  
} else {
  
  q5_table$breastmilk_freq_18m <- NA
  
}

q6_table$diabetes_3y <- NA

if ("GG49" %in% names(q6_raw_table) && "GG51" %in% names(q6_raw_table)) {
  
  q6_table$diabetes_3y[q6_raw_table$GG49 == 1] <- 0
  q6_table$diabetes_3y[q6_raw_table$GG50 == 1 | q6_raw_table$GG51 == 1] <- 1
  
}

q6_table$underweight_3y <- NA

if ("GG53" %in% names(q6_raw_table) && "GG54" %in% names(q6_raw_table) && "GG55" %in% names(q6_raw_table)) {
  
  q6_table$underweight_3y[q6_raw_table$GG53 == 1] <- 0
  q6_table$underweight_3y[q6_raw_table$GG54 == 1 | q6_raw_table$GG55 == 1] <- 1
  
}

q6_table$overweight_3y <- NA

if ("GG57" %in% names(q6_raw_table) && "GG58" %in% names(q6_raw_table) && "GG59" %in% names(q6_raw_table)) {
  
  q6_table$overweight_3y[q6_raw_table$GG57 == 1] <- 0
  q6_table$overweight_3y[q6_raw_table$GG58 == 1 | q6_raw_table$GG59 == 1] <- 1
  
}

q8_table$length_7y <- ifelse(!is.na(q8_raw_table$JJ408), q8_raw_table$JJ408, q8_raw_table$JJ324*100)


# Make a single data frame

nrow_mfr <- nrow(mfr_table)

rawPheno <- mfr_table %>%
  left_join(q1m_table, by = "preg_id") %>%
  left_join(q1f_table, by = "preg_id") %>%
  # left_join(q2_table, by = "preg_id") %>%
  left_join(q4_table, by = c("preg_id", "rank_siblings")) %>%
  left_join(q5_table, by = c("preg_id", "rank_siblings")) %>%
  left_join(q6_table, by = c("preg_id", "rank_siblings")) %>%
  left_join(q7_table, by = c("preg_id", "rank_siblings")) %>%
  left_join(q8_table, by = c("preg_id", "rank_siblings")) %>%
  left_join(q9_table, by = c("preg_id", "rank_siblings")) %>%
  left_join(kostUngdom_table, by = c("preg_id", "rank_siblings"))

# Add sentrix ids

rawPheno$mother_id[rawPheno$mother_id == ""] <- NA
rawPheno$father_id[rawPheno$father_id == ""] <- NA

rawPheno <- rawPheno %>% 
  left_join(
    childIdDF,
    by = c("preg_id", "rank_siblings")
  ) %>% 
  left_join(
    famDF,
    by = "sentrix_id"
  ) %>% 
  mutate(
    child_id = ifelse(is.na(child_id), paste0(preg_id, "_", rank_siblings), child_id)
  )

nrow_pheno <- nrow(rawPheno)

if (nrow_mfr != nrow_pheno) {
  
  stop("Duplictes were introduced when merging data frames.")
  
}

print(glue("Phenotypes loaded:"))
print(glue("- Children in birth registry: {nrow(rawPheno)}"))
print(glue("- Children genotyped: {sum(!is.na(rawPheno$sentrix_id))}"))
print(glue("- Mothers genotyped linked to a child: {length(unique(rawPheno$mother_sentrix_id))}"))
print(glue("- Fathers genotyped linked to a child: {length(unique(rawPheno$father_sentrix_id))}"))

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


# Convert to number

print(paste(Sys.time(), " Converting string input to number"))

rawPheno$age_birth <- 0

rawPheno$age_14y <- rawPheno$age_kost

rawPheno$weight_14y <- str_remove_all(rawPheno$weight_kost, " KG")
rawPheno$weight_14y <- str_remove_all(rawPheno$weight_14y, "MINDRE ENN ")
rawPheno$weight_14y <- str_remove_all(rawPheno$weight_14y, "MER ENN ")
rawPheno$weight_14y <- as.numeric(rawPheno$weight_14y)

rawPheno$length_14y <- str_remove_all(rawPheno$length_kost, " CM")
rawPheno$length_14y <- str_remove_all(rawPheno$length_14y, "LAVERE ENN ")
rawPheno$length_14y <- str_remove_all(rawPheno$length_14y, "HØYERE ENN ")
rawPheno$length_14y <- as.numeric(rawPheno$length_14y)


# Merge columns

print(paste(Sys.time(), " Combining phenotypes"))

if ("mother_height_self" %in% names(rawPheno) & "mother_height" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      mother_height = ifelse(is.na(mother_height), mother_height_self, mother_height)
    )
  
}

if ("mother_weight_beginning_self" %in% names(rawPheno) & "mother_weight_beginning" %in% names(rawPheno)) {
  
  rawPheno <- rawPheno %>%
    mutate(
      mother_weight_beginning = ifelse(is.na(mother_weight_beginning), mother_weight_beginning_self, mother_weight_beginning)
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


# Get a list of unrelated kids

print(paste(Sys.time(), " Loading identifiers for unrelated children"))

unrelatedDF <- read.table(
  file = unrelated_children_id_path,
  header = F,
  stringsAsFactors = F
)
unrelatedIds <- unrelatedDF[, 2]

rawPheno <- rawPheno %>%
  mutate(
    unrelated_children = ifelse(sentrix_id %in% unrelatedIds, 1, 0)
  )

print(glue("Unrelated children:"))
print(glue("- Children genotyped: {sum(!is.na(rawPheno$sentrix_id) & rawPheno$unrelated_children == 1)}"))
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
    sex = ifelse(registry_sex == "Mann", 1, sex),
    sex = ifelse(registry_sex == "Kvinne", 2, sex),
    sex = ifelse(!is.na(genetic_sex) & genetic_sex != 0, genetic_sex, sex)
  )


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

term_pregnancies <- sum(!is.na(rawPheno$sentrix_id) & rawPheno$pregnancy_duration_term == 1)
n_genotyped <- sum(!is.na(rawPheno$sentrix_id))

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
    age_14y = age_14y / 12 * 365.25,
    
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
    n_genotyped = sum(!is.na(rawPheno[[column]]) & !is.na(rawPheno$sentrix_id)),
    n_genotyped_unrelated = sum(!is.na(rawPheno[[column]]) & !is.na(rawPheno$sentrix_id) & rawPheno$unrelated == 1),
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
  
  mean_value <- mean(values[[column]][!is.na(values[[column]]) & !is.na(values$sentrix_id)])
  sd_value <- sd(values[[column]][!is.na(values[[column]]) & !is.na(values$sentrix_id)])
  
  toExclude <- !is.na(values[[column]]) & (values[[column]] < mean_value - 5 * sd_value | values[[column]] > mean_value + 5 * sd_value)
  values[[column]][toExclude] <- NA
  
  n_df <- data.frame(
    phenotype = column,
    n_outliers = sum(toExclude),
    n_outliers_genotyped = sum(toExclude & !is.na(values$sentrix_id)),
    n_outliers_genotyped_unrelated = sum(toExclude & !is.na(values$sentrix_id) & values$unrelated == 1),
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
      n_breast_feeding_genotyped = sum(!is.na(values[[column]]) & !is.na(values$sentrix_id)),
      n_breast_feeding_genotyped_unrelated = sum(!is.na(values[[column]]) & !is.na(values$sentrix_id) & values$unrelated == 1),
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
      n_breast_feeding_genotyped = sum(!is.na(values[[column]]) & !is.na(values$sentrix_id)),
      n_breast_feeding_genotyped_unrelated = sum(!is.na(values[[column]]) & !is.na(values$sentrix_id) & values$unrelated == 1),
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
      n_diabetes_genotyped = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$sentrix_id)),
      n_diabetes_unrelated = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$sentrix_id) & values$unrelated == 1),
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
      n_diabetes_genotyped = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$sentrix_id)),
      n_diabetes_unrelated = sum(!is.na(values[[column]]) & values[[column]] == 1 & !is.na(values$sentrix_id) & values$unrelated == 1),
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
    ),
    mother_median_weight = as.numeric(
      pmap(
        list(mother_weight_beginning, mother_weight_3y, mother_weight_5y, mother_weight_8y),
        ~median(c(..1, ..2, ..3, ..4), na.rm = T)
      )
    )
  )


# Save the variables in different tables

print(paste(Sys.time(), " Saving the phenotypes to tables"))

id_columns <- id_columns[id_columns %in% names(values)]

idValues <- values %>% 
  select(
    all_of(id_columns)
  )

write.table(
  x = idValues,
  file = gzfile(file.path(qcFolder, "ids.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

pregnancy_columns <- pregnancy_columns[pregnancy_columns %in% names(values)]
pregnancyValues <- values %>% 
  select(
    all_of(c(default_columns, pregnancy_columns))
  )

write.table(
  x = pregnancyValues,
  file = gzfile(file.path(qcFolder, "pregnancy.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

delivery_columns <- delivery_columns[delivery_columns %in% names(values)]
deliveryValues <- values %>% 
  select(
    all_of(c(default_columns, delivery_columns))
  )

write.table(
  x = deliveryValues,
  file = gzfile(file.path(qcFolder, "delivery.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

pregnancy_nutrition_columns <- pregnancy_nutrition_columns[pregnancy_nutrition_columns %in% names(values)]
pregnancyNutritionValues <- values %>% 
  select(
    all_of(c(default_columns, pregnancy_nutrition_columns))
  )

write.table(
  x = pregnancyNutritionValues,
  file = gzfile(file.path(qcFolder, "pregnancy_nutrition.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

mother_nutrition_columns <- mother_nutrition_columns[mother_nutrition_columns %in% names(values)]
motherNutritionValues <- values %>% 
  select(
    all_of(c(default_columns, mother_nutrition_columns))
  )

write.table(
  x = motherNutritionValues,
  file = gzfile(file.path(qcFolder, "mother_nutrition.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

child_nutrition_columns <- child_nutrition_columns[child_nutrition_columns %in% names(values)]
childNutritionValues <- values %>% 
  select(
    all_of(c(default_columns, child_nutrition_columns))
  )

write.table(
  x = childNutritionValues,
  file = gzfile(file.path(qcFolder, "child_nutrition.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

child_columns <- child_columns[child_columns %in% names(values)]
childValues <- values %>% 
  select(
    all_of(c(default_columns, child_columns))
  )

write.table(
  x = childValues,
  file = gzfile(file.path(qcFolder, "child.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

child_health_columns <- child_health_columns[child_health_columns %in% names(values)]
childHealth <- values %>% 
  select(
    all_of(c(default_columns, child_health_columns))
  )

write.table(
  x = childHealth,
  file = gzfile(file.path(qcFolder, "child_health.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

parent_values_columns <- parent_values_columns[parent_values_columns %in% names(values)]
parentValues <- values %>% 
  select(
    all_of(c(default_columns, parent_values_columns))
  )

write.table(
  x = parentValues,
  file = gzfile(file.path(qcFolder, "parents.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

mother_health_columns <- mother_health_columns[mother_health_columns %in% names(values)]
motherHealthValues <- values %>% 
  select(
    all_of(c(default_columns, mother_health_columns))
  )

write.table(
  x = motherHealthValues,
  file = gzfile(file.path(qcFolder, "mother_health.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

father_health_columns <- father_health_columns[father_health_columns %in% names(values)]
fatherHealthValues <- values %>% 
  select(
    all_of(c(default_columns, father_health_columns))
  )

write.table(
  x = fatherHealthValues,
  file = gzfile(file.path(qcFolder, "father_health.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

age_columns <- age_columns[age_columns %in% names(values)]
weight_columns <- weight_columns[weight_columns %in% names(values)]
length_columns <- length_columns[length_columns %in% names(values)]
head_circumference_columns <- head_circumference_columns[head_circumference_columns %in% names(values)]
childAnthropometricValues <- values %>% 
  select(
    all_of(c(default_columns, age_columns, weight_columns, length_columns, head_circumference_columns))
  )

write.table(
  x = childAnthropometricValues,
  file = gzfile(file.path(qcFolder, "child_anthropometrics_raw.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

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
