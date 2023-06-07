##
#
# This script processes phenotypes and exports files for genetic analyses and associated documentation.
# 
##

print(paste0(Sys.time(), " - Organize phenotypes in a clean table and write documentation"))



# Libararies, need to be available via the conda environment

library(conflicted)
library(foreign)
library(stringr)
library(glue)
library(crayon)
library(tidyr)
library(dplyr)
library(janitor)


# Avoid namespace conflicts

conflicts_prefer(dplyr::filter)



# Paths from command line
args <- commandArgs(TRUE)

rawTablesFolder <- args[1]
child_id_linkage_raw_table_path <- args[2]
mother_id_linkage_raw_table_path <- args[3]
father_id_linkage_raw_table_path <- args[4]
mfr_raw_table <- args[5]
msis_raw_table_path <- args[6]
child_msis_id_mapping_raw_table_path <- args[7]
mother_msis_id_mapping_raw_table_path <- args[8]
father_msis_id_mapping_raw_table_path <- args[9]
sysvak_raw_table_path <- args[10]
child_sysvak_id_mapping_raw_table_path <- args[11]
mother_sysvak_id_mapping_raw_table_path <- args[12]
father_sysvak_id_mapping_raw_table_path <- args[13]
covidTable <- args[14]
docsFolder <- args[15]
mobaProjectNumber <- args[16]

### DEBUG

# rawTablesFolder <- "/mnt/work/marc/pheno_covid_23-03-03/raw"
# child_id_linkage_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/linkage/20220516_MoBaGeneticsTot_Child_PDB2824.gz"
# mother_id_linkage_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/linkage/20220516_MoBaGeneticsTot_Mother_PDB2824.gz"
# father_id_linkage_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/linkage/20220516_MoBaGeneticsTot_Father_PDB2824.gz"
# mfr_raw_table <- "/mnt/work/marc/pheno_covid_23-03-03/raw/moba_ques/PDB2824_MFR_541_v12.gz"
# msis_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/msis/PDB2824_MSIS-data_MoBa.gz"
# child_msis_id_mapping_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/msis/Barn_ID_2824_2021_11_17sav.gz"
# mother_msis_id_mapping_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/msis/Mor_ID_2824_2021_11_17sav.gz"
# father_msis_id_mapping_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/msis/Far_ID_2824_2021_11_17sav.gz"
# sysvak_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/sysvak/SYSVAK210043_KOBLET_MOBA_01022022.gz"
# child_sysvak_id_mapping_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/sysvak/2022_02_01_Barn_koblingsbro_2824.gz"
# mother_sysvak_id_mapping_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/sysvak/2022_02_01_Mor_koblingsbro_2824.gz"
# father_sysvak_id_mapping_raw_table_path <- "/mnt/work/marc/pheno_covid_23-03-03/raw/sysvak/2022_02_01_Far_koblingsbro_2824_.gz"
# covidTable <- "/mnt/work/marc/pheno_covid_23-03-03/covid/moba_covid_phenotypes.gz"
# docsFolder <- "docs/covid/23-03-03/covid"
# mobaProjectNumber <- 2824

###

# Functions

source("src/covid/scripts/utils/long_covid_docs.R")
source("src/covid/scripts/utils/variables_mapping.R")


# Parameters

batchOrder <- c("NORMENT_JAN21", "NORMENT_MAI21", "NORMENT_SEP20_R996R1029", "NORMENT-JAN20", "NORMENT-FEB18", "NORMENT-JUN15", "NORMENT-MAY16", "NORMENT-JAN15", "ROTTERDAM1", "ROTTERDAM2", "HARVEST", "TED", "PDB1382_R875_R876")

theme_set(theme_bw(base_size = 16))

updateDocs <- T


# Paths

docsFile <- file.path(docsFolder, "phenos.md")

longCovidDocsFolder <- file.path(docsFolder, "long_covid")
longCovidDocsFile <- file.path(longCovidDocsFolder, "long_covid.md")

quesFolder <- file.path(rawTablesFolder, "covid_ques")


# Housekeeping

if (!dir.exists(longCovidDocsFolder)) {
  
  dir.create(
    path = longCovidDocsFolder,
    showWarnings = T,
    recursive = T
  )
  
}

# Project specific variables

preg_id_column <- paste0("preg_id_", mobaProjectNumber)
mother_id_column <- paste0("m_id_", mobaProjectNumber)
father_id_column <- paste0("f_id_", mobaProjectNumber)
parent_id_column <- paste0("p_id_", mobaProjectNumber)
kIdColumn <- paste0("pid_k_", mobaProjectNumber)
msisIdColumn <- paste0("pid_msis_", mobaProjectNumber)


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
  mutate(
    child_id = paste0(!!sym(preg_id_column), "_", barn_nr)
  )

motherIdDF <- read.table(
  file = mother_id_linkage_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()


fatherIdDF <- read.table(
  file = father_id_linkage_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()

idDF <- rbind(
  childIdDF %>% 
    select(
      id = child_id, sentrix_id, role, batch
    ),
  motherIdDF %>% 
    select(
      id = !!sym(mother_id_column), sentrix_id, role, batch
    ),
  fatherIdDF %>% 
    select(
      id = !!sym(father_id_column), sentrix_id, role, batch
    )
)

if (sum(is.na(idDF$sentrix_id)) > 0) {
  
  stop("Missing sentrix id")
  
}


# Check that all batches are supported

if (sum(!idDF$batch %in% batchOrder) > 0) {
  
  missingBatches <- unique(idDF$batch[!idDF$batch %in% batchOrder])
  
  stop(paste0("Batch not supported: ", paste(missingBatches, collapse = ", ")))
  
}


# Remove duplicates

print(glue("{Sys.time()} - Removing duplicates"))

nBefore <- nrow(idDF)

sampleOccurrenceDF <- as.data.frame(
  table(idDF$id),
  stringsAsFactors = F
) %>% 
  clean_names() %>% 
  filter(
    freq > 1
  )

for (duplicate in sampleOccurrenceDF$var1) {
  
  duplicateIds <- idDF %>% 
    filter(
      id == duplicate
    ) %>% 
    mutate(
      batchFactor = factor(batch, levels = batchOrder)
    ) %>% 
    arrange(
      batchFactor
    )
  
  idDF <- idDF %>% 
    filter(
      id != duplicate | sentrix_id == duplicateIds$sentrix_id[1]
    )
  
}

if (length(unique(idDF$id)) != length(unique(idDF$sentrix_id))) {
  
  stop("Duplicate sentrix_ids found")
  
}

if (nrow(idDF) != length(unique(idDF$id))) {
  
  stop("Duplicate ids found")
  
}

nAfter <- nrow(idDF)

print(paste0("Duplicates removed: ", nBefore - nAfter, " (", round(100 * (nBefore - nAfter)/nBefore, 1), "%)."))


# Load vaccine information

sysvak_table <- read.table(
  file = sysvak_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()

child_sysvak_id_mapping <- read.table(
  file = child_sysvak_id_mapping_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()

mother_sysvak_id_mapping <- read.table(
  file = mother_sysvak_id_mapping_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()

father_sysvak_id_mapping <- read.table(
  file = father_sysvak_id_mapping_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()

child_sysvak <- sysvak_table %>% 
  left_join(
    child_sysvak_id_mapping %>%
      mutate(
        id = paste0(!!sym(preg_id_column), "_", barn_nr)
      ) %>% 
      select(
        id, !!sym(kIdColumn)
      ),
    by = kIdColumn
  ) %>% 
  filter(
    !is.na(id)
  ) %>% 
  select(
    id, 
    vaccination_date = konsultasjonsdato, 
    vaccine_code = vaksinekode_vasket
  )

mother_sysvak <- sysvak_table %>% 
  left_join(
    mother_sysvak_id_mapping %>%
      select(
        id = !!sym(mother_id_column), !!sym(kIdColumn)
      ),
    by = kIdColumn
  ) %>% 
  filter(
    !is.na(id)
  ) %>% 
  select(
    id, 
    vaccination_date = konsultasjonsdato, 
    vaccine_code = vaksinekode_vasket
  )

father_sysvak <- sysvak_table %>% 
  left_join(
    father_sysvak_id_mapping %>%
      select(
        id = !!sym(father_id_column), !!sym(kIdColumn)
      ),
    by = kIdColumn
  ) %>% 
  filter(
    !is.na(id)
  ) %>% 
  select(
    id, 
    vaccination_date = konsultasjonsdato, 
    vaccine_code = vaksinekode_vasket
  )

vaccination_table <- rbind(
  child_sysvak, mother_sysvak, father_sysvak
) %>% 
  group_by(
    id, vaccination_date
  ) %>% 
  summarize(
    vaccine_code = paste(sort(vaccine_code), collapse = ","),
    .groups = "drop"
  ) %>% 
  mutate(
    vaccination_date = as.Date(vaccination_date, "%d.%m.%Y")
  ) %>% 
  arrange(
    vaccination_date
  ) %>% 
  group_by(
    id
  ) %>% 
  mutate(
    dose_number = row_number()
  ) %>% 
  ungroup()


# Set up a data frame for phenotypes

print(glue("{Sys.time()} - Setting up pheno table"))

pheno_variable_to_question <- list()

phenoDF <- idDF %>%
  mutate(
    sick_past_14_days = NA,
    suspected_or_confirmed_covid_doctor_past_14_days = 0,
    suspected_or_confirmed_covid_doctor_past_14_days_last_reported = NA,
    tested_positive = 0,
    tested_positive_last_reported = NA,
    tested_positive_pcr = NA,
    tested_positive_pcr_last_reported = NA,
    tested_positive_ab = NA,
    tested_positive_ab_last_reported = NA,
    menstruating = NA,
    pill = NA,
    finished_menstruating = NA
  )


# Effect of vaccines

for (i in 1:length(influenza_variables)) {
  
  variable <- influenza_variables[i]
  question <- influenza_questions[i]
  
  phenoDF[[variable]] <- NA
  pheno_variable_to_question[[variable]] <- c(question)
  
}
for (i in 1:length(covid_vaccination_variables)) {
  
  variable <- covid_vaccination_variables[i]
  question <- covid_vaccination_questions[i]
  
  phenoDF[[variable]] <- NA
  pheno_variable_to_question[[variable]] <- c(question)
  
}
for (i in 1:length(covid_vaccination_menstruation_variables)) {
  
  variable <- covid_vaccination_menstruation_variables[i]
  question <- covid_vaccination_menstruation_questions[i]
  
  phenoDF[[variable]] <- NA
  pheno_variable_to_question[[variable]] <- c(question)
  
}

for (original_variable in covid_vaccination_variables) {
  
  phenoDF[[original_variable]] <- NA
  
  if (endsWith(original_variable, "_first_dose")) {
    
    variable <- substr(original_variable, 1, nchar(original_variable) - nchar("_first_dose"))
    
  } else if (endsWith(original_variable, "_last_dose")) {
    
    variable <- substr(original_variable, 1, nchar(original_variable) - nchar("_last_dose"))
    
  }
  
  pheno_variable_to_question[[variable]] <- unique(c(pheno_variable_to_question[[variable]], pheno_variable_to_question[[original_variable]]))
  
  new_variable <- paste0(variable, "_after_bnt")
  phenoDF[[new_variable]] <- NA
  pheno_variable_to_question[[new_variable]] <- unique(c(pheno_variable_to_question[[new_variable]], pheno_variable_to_question[[variable]]))
  
  new_variable <- paste0(variable, "_after_mod")
  phenoDF[[new_variable]] <- NA
  pheno_variable_to_question[[new_variable]] <- unique(c(pheno_variable_to_question[[new_variable]], pheno_variable_to_question[[variable]]))
  
  new_variable <- paste0(variable, "_after_asz")
  phenoDF[[new_variable]] <- NA
  pheno_variable_to_question[[new_variable]] <- unique(c(pheno_variable_to_question[[new_variable]], pheno_variable_to_question[[variable]]))
  
  new_variable <- paste0(variable, "_after_sin")
  phenoDF[[new_variable]] <- NA
  pheno_variable_to_question[[new_variable]] <- unique(c(pheno_variable_to_question[[new_variable]], pheno_variable_to_question[[variable]]))
  
  new_variable <- paste0(variable, "_after_jan")
  phenoDF[[new_variable]] <- NA
  pheno_variable_to_question[[new_variable]] <- unique(c(pheno_variable_to_question[[new_variable]], pheno_variable_to_question[[variable]]))
  
}


# Long covid

print(glue("{Sys.time()} - Setting up long covid function"))

longCovidPhenoImport <- function(
    quesDF,
    phenoDF,
    quesNumber,
    phenoName,
    folder
) {
  
  if (!phenoName %in% names(phenoDF)) {
    
    phenoDF[[phenoName]] <- NA
    
  }
  
  phenoLast <- paste0(phenoName, "_last_reported")
  
  if (!phenoLast %in% names(phenoDF)) {
    
    phenoDF[[phenoLast]] <- NA
    pheno_variable_to_question[[phenoLast]] <- pheno_variable_to_question[[phenoName]]
    
  }
  
  tempDF <- quesDF %>% 
    filter(
      !is.na(!!sym(quesNumber)) & !!sym(quesNumber) == 1
    )
  
  if (nrow(tempDF) > 0) {
    
    tempDF <- tempDF %>% 
      select(
        id, !!sym(quesNumber), fill_in_date
      ) %>% 
      group_by(
        id
      ) %>% 
      arrange(
        desc(!!sym(quesNumber))
      ) %>% 
      filter(
        row_number() == 1
      )
    
    if (nrow(tempDF) != length(unique(tempDF$id))) {
      
      stop(glue("Duplicate ids found in {folder}"))
      
    }
    
    phenoDF <- phenoDF %>% 
      left_join(
        tempDF,
        by = "id"
      ) %>% 
      mutate(
        !!phenoName := ifelse(
          is.na(!!sym(phenoName)) | !is.na(!!sym(phenoName)) & !is.na(!!sym(quesNumber)) & !!sym(phenoName) == 0, 
          !!sym(quesNumber), 
          !!sym(phenoName)
        ),
        !!phenoLast := ifelse(
          !is.na(!!sym(phenoName)) & !!sym(phenoName) == 1 & is.na(!!sym(phenoLast)) | !is.na(!!sym(phenoLast)) & !is.na(!!sym(quesNumber)) & !!sym(quesNumber) == 1 & fill_in_date > !!sym(phenoLast), 
          fill_in_date, 
          !!sym(phenoLast)
        )
      ) %>% 
      select(
        -!!sym(quesNumber), -fill_in_date
      )
    
  }
  
  if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
    
    stop(glue("Duplicate sentrix id introduced during processing of {phenoName}."))
    
  }
  
  return(phenoDF)
  
}


# Load covid questionnaire data

print(glue("{Sys.time()} - Loading long covid questionnaires"))

quesNames <- c("Foreldre_duplikater_{suffix}.gz", "Foreldre_{suffix}.gz", "Ungdom_duplikater_{suffix}.gz", "Ungdom_{suffix}.gz")

for (folder in list.files(quesFolder)) {
  
  if (folder != "Dokumentasjon") {
    
    print(glue("{Sys.time()} - Loading phenotypes from {folder}"))
    
    if (folder == "2020_06_07_Runde_1_4_MedDup") {
      
      suffixes <- "runde1til4"
      
    } else if (folder == "Runde13") {
      
      suffixes <- paste0(folder, "_korrigert17032021")
      
    } else if (folder == "Runde14") {
      
      suffixes <- paste0(folder, "_korrigert17032021")
      
    } else if (folder == "Runde16") {
      
      suffixes <- c(paste0(folder, "_komplett"), paste0(folder, "_09112020"))
      
    } else if (folder == "Runde22") {
      
      suffixes <- paste0(folder, "_komplett_korrigert17032021")
      
    } else if (folder == "Runde17" || folder == "Runde18" || folder == "Runde19" || folder == "Runde20" || folder == "Runde21" || folder == "Runde23" || folder == "Runde24" || folder == "Runde25" || folder == "Runde26" || folder == "Runde27" || folder == "Runde28" || folder == "Runde29" || folder == "Runde30" || folder == "Runde31" || folder == "Runde32" || folder == "Runde33" || folder == "Runde34" || folder == "Runde35" || folder == "Runde36" || folder == "Runde37" || folder == "Runde38" || folder == "Runde39" || folder == "Runde40" || folder == "Runde41" || folder == "Runde42" || folder == "Runde43" || folder == "Runde44") {
      
      suffixes <- paste0(folder, "_komplett")
      
    } else {
      
      suffixes <- folder
      
    }
    
    for (quesName in quesNames) {
      
      for (suffix in suffixes) {
        
        questionnaireFile <- file.path(quesFolder, folder, glue(quesName))
        
        # print(paste0(questionnaireFile, " - Found: ", file.exists(questionnaireFile)))
        
        if (file.exists(questionnaireFile)) {
          
          # print(glue("{Sys.time()} - Loading phenotypes from {questionnaireFile}"))
          
          quesDF <- read.table(
            file = questionnaireFile,
            sep = "\t",
            header = T,
            quote = "",
            stringsAsFactors = F,
            comment.char = ""
          ) %>% 
            clean_names() %>% 
            mutate(
              role = str_replace_all(
                string = role, 
                pattern = " ",
                repl = ""
              ),
              id = ifelse(
                test = role == "Child",
                yes = paste0(!!sym(preg_id_column), "_", barn_nr),
                no = !!sym(parent_id_column)
              )
            )
          
          if (!"fill_in_date" %in% names(quesDF)) {
            
            stop("'fill_in_date' not found in '{questionnaireFile}': {names(quesDF)}")
            
          }
          
          if ("kf10" %in% names(quesDF)) {
            
            tempDF <- quesDF %>% 
              mutate(
                kf10 = as.numeric(factor(kf10, levels = c("NEI", "JA"))) - 1
              ) %>% 
              select(
                id, kf10
              ) %>% 
              group_by(
                id
              ) %>% 
              summarize(
                kf10 = max(kf10, na.rm = T), 
                .groups = 'drop'
              ) %>% 
              mutate(
                kf10 = ifelse(is.infinite(kf10), NA, kf10)
              )
            
            if (sum(!is.na(tempDF$kf10)) > 0) {
              
              if (nrow(tempDF) != length(unique(tempDF$id))) {
                
                stop(glue("Duplicate ids found in {folder}"))
                
              }
              
              phenoDF <- phenoDF %>% 
                left_join(
                  tempDF,
                  by = "id"
                ) %>% 
                mutate(
                  sick_past_14_days = case_when(
                    is.na(sick_past_14_days) ~ kf10,
                    !is.na(sick_past_14_days) & !is.na(kf10) & sick_past_14_days == 0 ~ kf10,
                    T ~ as.numeric(sick_past_14_days)
                  )
                ) %>% 
                select(
                  -kf10
                )
              pheno_variable_to_question[["sick_past_14_days"]] <- unique(c(pheno_variable_to_question[["sick_past_14_days"]], "kf10"))
              
            }
            
            if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
              
              stop("Duplicate sentrix id introduced during processing of 'kf10'.")
              
            }
          }
          if ("kf30" %in% names(quesDF)) {
            
            caseDF <- quesDF %>% 
              filter(
                !is.na(kf30)
              ) %>% 
              select(
                id, kf30, fill_in_date
              ) %>% 
              group_by(
                id
              ) %>% 
              summarize(
                kf30 = max(kf30, na.rm = T), 
                fill_in_date = max(fill_in_date, na.rm = T), 
                .groups = 'drop'
              ) %>%
              filter(
                !is.na(kf30) & !is.infinite(kf30)
              ) %>% 
              select(
                id, fill_in_date
              )
            
            phenoDF <- phenoDF %>% 
              left_join(
                caseDF,
                by = "id"
              ) %>% 
              mutate(
                suspected_or_confirmed_covid_doctor_past_14_days = ifelse(id %in% caseDF$id, 1, suspected_or_confirmed_covid_doctor_past_14_days),
                suspected_or_confirmed_covid_doctor_past_14_days_last_reported = ifelse(is.na(suspected_or_confirmed_covid_doctor_past_14_days_last_reported) | !is.na(fill_in_date) & fill_in_date > suspected_or_confirmed_covid_doctor_past_14_days_last_reported, fill_in_date, suspected_or_confirmed_covid_doctor_past_14_days_last_reported),
              ) %>% 
              select(
                -fill_in_date
              )
            pheno_variable_to_question[["suspected_or_confirmed_covid_doctor_past_14_days"]] <- unique(c(pheno_variable_to_question[["suspected_or_confirmed_covid_doctor_past_14_days"]], "kf30"))
            pheno_variable_to_question[["suspected_or_confirmed_covid_doctor_past_14_days_last_reported"]] <- unique(c(pheno_variable_to_question[["suspected_or_confirmed_covid_doctor_past_14_days_last_reported"]], "kf30"))
            
            if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
              
              stop("Duplicate sentrix id introduced during processing of 'kf30'.")
              
            }
          }
          if ("kf41" %in% names(quesDF)) {
            
            caseDF <- quesDF %>% 
              filter(
                !is.na(kf41) & kf41 == "JA"
              ) %>% 
              select(
                id, fill_in_date
              ) %>% 
              group_by(
                id
              ) %>% 
              summarize(
                fill_in_date = max(fill_in_date, na.rm = T), 
                .groups = 'drop'
              )
            
            phenoDF <- phenoDF %>% 
              left_join(
                caseDF,
                by = "id"
              ) %>% 
              mutate(
                tested_positive = ifelse(id %in% quesDF$id[!is.na(quesDF$kf41) & quesDF$kf41 == "JA"], 1, tested_positive),
                tested_positive_last_reported = ifelse(is.na(tested_positive_last_reported) | !is.na(fill_in_date) & fill_in_date > tested_positive_last_reported, fill_in_date, tested_positive_last_reported),
              ) %>% 
              select(
                -fill_in_date
              )
            pheno_variable_to_question[["tested_positive"]] <- unique(c(pheno_variable_to_question[["tested_positive"]], "kf41"))
            pheno_variable_to_question[["tested_positive_last_reported"]] <- unique(c(pheno_variable_to_question[["tested_positive_last_reported"]], "kf41"))
            
            if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
              
              stop("Duplicate sentrix id introduced during processing of 'kf41'.")
              
            }
          }
          if ("kf165" %in% names(quesDF)) {
            
            tempDF <- quesDF %>% 
              mutate(
                kf165 = as.numeric(factor(kf165, levels = c("NEI", "JA"))) - 1
              ) %>% 
              mutate(
                kf165 = ifelse(is.infinite(kf165), NA, kf165)
              ) %>% 
              filter(
                !is.na(kf165)
              ) 
            
            if (nrow(tempDF) > 0) {
              
              tempDF <- tempDF %>% 
                select(
                  id, kf165, fill_in_date
                ) %>% 
                group_by(
                  id
                ) %>% 
                arrange(
                  desc(kf165)
                ) %>% 
                filter(
                  row_number() == 1
                )
              
              if (nrow(tempDF) != length(unique(tempDF$id))) {
                
                stop(glue("Duplicate ids found in {folder}"))
                
              }
              
              phenoDF <- phenoDF %>% 
                left_join(
                  tempDF,
                  by = "id"
                ) %>% 
                mutate(
                  tested_positive_pcr = case_when(
                    is.na(tested_positive_pcr) ~ kf165,
                    !is.na(tested_positive_pcr) & !is.na(kf165) & tested_positive_pcr == 0 ~ kf165,
                    T ~ as.numeric(tested_positive_pcr)
                  ),
                  tested_positive_pcr_last_reported = ifelse(
                    !is.na(tested_positive_pcr) & tested_positive_pcr == 1 & is.na(tested_positive_pcr_last_reported) | !is.na(tested_positive_pcr_last_reported) & !is.na(kf165) & kf165 == 1 & fill_in_date > tested_positive_pcr_last_reported, 
                    fill_in_date, 
                    tested_positive_pcr_last_reported
                  )
                ) %>% 
                select(
                  -kf165, -fill_in_date
                )
              pheno_variable_to_question[["tested_positive_pcr"]] <- unique(c(pheno_variable_to_question[["tested_positive_pcr"]], "kf165"))
              pheno_variable_to_question[["tested_positive_pcr_last_reported"]] <- unique(c(pheno_variable_to_question[["tested_positive_pcr_last_reported"]], "kf165"))
            }
            
            if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
              
              stop("Duplicate sentrix id introduced during processing of 'kf165'.")
              
            }
          }
          if ("kf305" %in% names(quesDF)) {
            
            tempDF <- quesDF %>% 
              mutate(
                kf305 = as.numeric(factor(kf305, levels = c("NEI", "JA"))) - 1
              ) %>% 
              mutate(
                kf305 = ifelse(is.infinite(kf305), NA, kf305)
              ) %>% 
              filter(
                !is.na(kf305)
              )
            
            if (nrow(tempDF) > 0) {
              
              tempDF <- tempDF %>% 
                select(
                  id, kf305, fill_in_date
                ) %>% 
                group_by(
                  id
                ) %>% 
                arrange(
                  desc(kf305)
                ) %>% 
                filter(
                  row_number() == 1
                )
              
              if (nrow(tempDF) != length(unique(tempDF$id))) {
                
                stop(glue("Duplicate ids found in {folder}"))
                
              }
              
              phenoDF <- phenoDF %>% 
                left_join(
                  tempDF,
                  by = "id"
                ) %>% 
                mutate(
                  tested_positive_pcr = case_when(
                    is.na(tested_positive_pcr) ~ kf305,
                    !is.na(tested_positive_pcr) & !is.na(kf305) & tested_positive_pcr == 0 ~ kf305,
                    T ~ as.numeric(tested_positive_pcr)
                  ),
                  tested_positive_pcr_last_reported = ifelse(
                    !is.na(tested_positive_pcr) & tested_positive_pcr == 1 & is.na(tested_positive_pcr_last_reported) | !is.na(tested_positive_pcr_last_reported) & !is.na(kf305) & kf305 == 1 & fill_in_date > tested_positive_pcr_last_reported, 
                    fill_in_date, 
                    tested_positive_pcr_last_reported
                  )
                ) %>% 
                select(
                  -kf305, -fill_in_date
                )
              pheno_variable_to_question[["tested_positive_pcr"]] <- unique(c(pheno_variable_to_question[["tested_positive_pcr"]], "kf305"))
              pheno_variable_to_question[["tested_positive_pcr_last_reported"]] <- unique(c(pheno_variable_to_question[["tested_positive_pcr_last_reported"]], "kf305"))
              
            }
            
            if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
              
              stop("Duplicate sentrix id introduced during processing of 'kf305'.")
              
            }
          }
          if ("kf166" %in% names(quesDF)) {
            
            tempDF <- quesDF %>% 
              mutate(
                kf166 = as.numeric(factor(kf166, levels = c("NEI", "JA"))) - 1
              ) %>% 
              mutate(
                kf166 = ifelse(is.infinite(kf166), NA, kf166)
              ) %>% 
              filter(
                !is.na(kf166)
              )
            
            if (nrow(tempDF) > 0) {
              
              tempDF <- tempDF %>% 
                select(
                  id, kf166, fill_in_date
                ) %>% 
                group_by(
                  id
                ) %>% 
                arrange(
                  desc(kf166)
                ) %>% 
                filter(
                  row_number() == 1
                )
              
              if (nrow(tempDF) != length(unique(tempDF$id))) {
                
                stop(glue("Duplicate ids found in {folder}"))
                
              }
              
              phenoDF <- phenoDF %>% 
                left_join(
                  tempDF,
                  by = "id"
                ) %>% 
                mutate(
                  tested_positive_ab = case_when(
                    is.na(tested_positive_ab) ~ kf166,
                    !is.na(tested_positive_ab) & !is.na(kf166) & tested_positive_ab == 0 ~ kf166,
                    T ~ as.numeric(tested_positive_ab)
                  ),
                  tested_positive_ab_last_reported = ifelse(
                    !is.na(tested_positive_ab) & tested_positive_ab == 1 & is.na(tested_positive_ab_last_reported) | !is.na(tested_positive_ab_last_reported) & !is.na(kf166) & kf166 == 1 & fill_in_date > tested_positive_ab_last_reported, 
                    fill_in_date, 
                    tested_positive_ab_last_reported
                  )
                ) %>% 
                select(
                  -kf166, -fill_in_date
                )
              pheno_variable_to_question[["tested_positive_ab"]] <- unique(c(pheno_variable_to_question[["tested_positive_ab"]], "kf166"))
              pheno_variable_to_question[["tested_positive_ab_last_reported"]] <- unique(c(pheno_variable_to_question[["tested_positive_ab_last_reported"]], "kf166"))
              
            }
            
            if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
              
              stop("Duplicate sentrix id introduced during processing of 'kf166'.")
              
            }
          }
          
          for (question_i in 1:length(long_covid_questions)) {
            
            question <- long_covid_questions[question_i]
            variable <- long_covid_variables[question_i]
            
            if (question %in% names(quesDF)) {
              
              if (question == "kf120") {
                
                quesDF$kf120 <- as.numeric(factor(quesDF$kf120, levels = c("NEI", "JA"))) - 1
                
              }
              
              phenoDF <- longCovidPhenoImport(
                quesDF = quesDF,
                phenoDF = phenoDF,
                quesNumber = question,
                phenoName = variable,
                folder = folder
              )
              
            }
          }
          
          
          # Influenza vaccination
          
          if ("kf2051" %in% names(quesDF)) {
            
            vaccinated <- quesDF$id[!is.na(quesDF$kf2051) & quesDF$kf2051 == "JA"]
            
            for (variable in influenza_variables) {
              
              phenoDF[[variable]][is.na(phenoDF[[variable]]) & phenoDF$id %in% vaccinated] <- 0
              
            }
          }
          
          for (question_i in 1:length(influenza_questions)) {
            
            question <- influenza_questions[question_i]
            variable <- influenza_variables[question_i]
            
            if (question %in% names(quesDF)) {
              
              ids <- quesDF$id[!is.na(quesDF[[question]]) & quesDF[[question]] != "NEI"]
              
              phenoDF[[variable]][phenoDF$id %in% ids] <- 1
              
            }
          }
          
          
          # Corona vaccination
          
          for (question_i in 1:length(covid_vaccination_variables)) {
            
            question <- covid_vaccination_questions[question_i]
            variable <- covid_vaccination_variables[question_i]
            
            if (question %in% names(quesDF)) {
              
              ids_nei <- quesDF$id[!is.na(quesDF[[question]]) & quesDF[[question]] == "NEI" & quesDF$id %in% phenoDF$id]
              ids_ja <- quesDF$id[!is.na(quesDF[[question]]) & quesDF[[question]] == "JA" & quesDF$id %in% phenoDF$id]
              
              phenoDF[[variable]][is.na(phenoDF[[variable]]) & phenoDF$id %in% ids_nei] <- 0
              phenoDF[[variable]][phenoDF$id %in% ids_ja] <- 1
              
              if (endsWith(variable, "_first_dose")) {
                
                variable_name <- substr(variable, 1, nchar(variable) - nchar("_first_dose"))
                
                new_variable <- paste0(variable_name, "_after_bnt")
                ids_vac <- vaccination_table$id[vaccination_table$vaccine_code == "BNT03" & vaccination_table$dose_number == 1]
                phenoDF[[new_variable]][is.na(phenoDF[[new_variable]]) & phenoDF$id %in% ids_nei & phenoDF$id %in% ids_vac] <- 0
                phenoDF[[new_variable]][phenoDF$id %in% ids_ja & phenoDF$id %in% ids_vac] <- 1
                
                new_variable <- paste0(variable_name, "_after_mod")
                ids_vac <- vaccination_table$id[vaccination_table$vaccine_code == "MOD03" & vaccination_table$dose_number == 1]
                phenoDF[[new_variable]][is.na(phenoDF[[new_variable]]) & phenoDF$id %in% ids_nei & phenoDF$id %in% ids_vac] <- 0
                phenoDF[[new_variable]][phenoDF$id %in% ids_ja & phenoDF$id %in% ids_vac] <- 1
                
                new_variable <- paste0(variable_name, "_after_asz")
                ids_vac <- vaccination_table$id[vaccination_table$vaccine_code == "ASZ03" & vaccination_table$dose_number == 1]
                phenoDF[[new_variable]][is.na(phenoDF[[new_variable]]) & phenoDF$id %in% ids_nei & phenoDF$id %in% ids_vac] <- 0
                phenoDF[[new_variable]][phenoDF$id %in% ids_ja & phenoDF$id %in% ids_vac] <- 1
                
                new_variable <- paste0(variable_name, "_after_sin")
                ids_vac <- vaccination_table$id[vaccination_table$vaccine_code == "SIN03" & vaccination_table$dose_number == 1]
                phenoDF[[new_variable]][is.na(phenoDF[[new_variable]]) & phenoDF$id %in% ids_nei & phenoDF$id %in% ids_vac] <- 0
                phenoDF[[new_variable]][phenoDF$id %in% ids_ja & phenoDF$id %in% ids_vac] <- 1
                
                new_variable <- paste0(variable_name, "_after_jan")
                ids_vac <- vaccination_table$id[vaccination_table$vaccine_code == "JAN03" & vaccination_table$dose_number == 1]
                phenoDF[[new_variable]][is.na(phenoDF[[new_variable]]) & phenoDF$id %in% ids_nei & phenoDF$id %in% ids_vac] <- 0
                phenoDF[[new_variable]][phenoDF$id %in% ids_ja & phenoDF$id %in% ids_vac] <- 1
                
              } else if (endsWith(variable, "_last_dose")) {
                
                variable_name <- substr(variable, 1, nchar(variable) - nchar("_last_dose"))
                
                for (id_nei in ids_nei) {
                  
                  ques_i <- which(quesDF$id == id_nei)
                  ref_date <- as.Date(quesDF$fill_in_date[ques_i]/86400, origin = "1582-10-14")
                  
                  id_vaccines <- vaccination_table %>% 
                    filter(
                      id == id_nei & vaccination_date <= ref_date
                    ) %>% 
                    arrange(
                      vaccination_date
                    )
                  
                  if (nrow(id_vaccines) > 0) {
                    
                    pheno_id_i <- which(phenoDF$id == id_nei)
                    vaccine_code <- id_vaccines$vaccine_code[nrow(id_vaccines)]
                    
                    if (vaccine_code == "BNT03") {
                      
                      new_variable <- paste0(variable_name, "_after_bnt")
                      
                      if (is.na(phenoDF[[new_variable]][pheno_id_i])) {
                      
                      phenoDF[[new_variable]][pheno_id_i] <- 0
                      
                      }
                      
                    } else if (vaccine_code == "MOD03") {
                      
                      new_variable <- paste0(variable_name, "_after_mod")
                      
                      if (is.na(phenoDF[[new_variable]][pheno_id_i])) {
                        
                        phenoDF[[new_variable]][pheno_id_i] <- 0
                        
                      }
                      
                    } else if (vaccine_code == "ASZ03") {
                      
                      new_variable <- paste0(variable_name, "_after_asz")
                      
                      if (is.na(phenoDF[[new_variable]][pheno_id_i])) {
                        
                        phenoDF[[new_variable]][pheno_id_i] <- 0
                        
                      }
                      
                    } else if (vaccine_code == "SIN03") {
                      
                      new_variable <- paste0(variable_name, "_after_sin")
                      
                      if (is.na(phenoDF[[new_variable]][pheno_id_i])) {
                        
                        phenoDF[[new_variable]][pheno_id_i] <- 0
                        
                      }
                      
                    } else if (vaccine_code == "JAN03") {
                      
                      new_variable <- paste0(variable_name, "_after_jan")
                      
                      if (is.na(phenoDF[[new_variable]][pheno_id_i])) {
                        
                        phenoDF[[new_variable]][pheno_id_i] <- 0
                        
                      }
                    }
                  }
                }
                
                for (id_ja in ids_ja) {
                  
                  ques_i <- which(quesDF$id == id_ja)
                  ref_date <- as.Date(quesDF$fill_in_date[ques_i]/86400, origin = "1582-10-14")
                  
                  id_vaccines <- vaccination_table %>% 
                    filter(
                      id == id_ja & vaccination_date <= ref_date
                    ) %>% 
                    arrange(
                      vaccination_date
                    )
                  
                  if (nrow(id_vaccines) > 0) {
                    
                    pheno_id_i <- which(phenoDF$id == id_ja)
                    vaccine_code <- id_vaccines$vaccine_code[nrow(id_vaccines)]
                    
                    if (vaccine_code == "BNT03") {
                      
                      new_variable <- paste0(variable_name, "_after_bnt")
                      phenoDF[[new_variable]][pheno_id_i] <- 1
                      
                    } else if (vaccine_code == "MOD03") {
                      
                      new_variable <- paste0(variable_name, "_after_mod")
                      phenoDF[[new_variable]][pheno_id_i] <- 1
                      
                    } else if (vaccine_code == "ASZ03") {
                      
                      new_variable <- paste0(variable_name, "_after_asz")
                      phenoDF[[new_variable]][pheno_id_i] <- 1
                      
                    } else if (vaccine_code == "SIN03") {
                      
                      new_variable <- paste0(variable_name, "_after_sin")
                      phenoDF[[new_variable]][pheno_id_i] <- 1
                      
                    } else if (vaccine_code == "JAN03") {
                      
                      new_variable <- paste0(variable_name, "_after_jan")
                      phenoDF[[new_variable]][pheno_id_i] <- 1
                      
                    }
                  }
                }
              }
            }
          }
          
          
          # Corona menstruation
          
          for (question_i in 1:length(covid_vaccination_menstruation_questions)) {
            
            question <- covid_vaccination_menstruation_questions[question_i]
            variable <- covid_vaccination_menstruation_variables[question_i]
            
            if (question %in% names(quesDF)) {
              
              ids <- quesDF$id[!is.na(quesDF[[question]]) & quesDF[[question]] != "JA"]
              
              phenoDF[[variable]][phenoDF$id %in% ids] <- 1
              
              ids <- quesDF$id[!is.na(quesDF[[question]]) & quesDF[[question]] != "NEI"]
              
              phenoDF[[variable]][is.na(phenoDF[[variable]]) & phenoDF$id %in% ids] <- 0
              
            }
          }
        }
      }
    }
  }
}

if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
  
  stop("Duplicate sentrix id introduced during processing of covid questionnaires.")
  
}

# Load medical birth registry

print(glue("{Sys.time()} - Loading MBR data"))

mbrDF <- read.table(
  file = mfr_raw_table,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  mutate(
    child_id = paste0(!!sym(preg_id_column), "_", barn_nr)
  )

child_faarDF <- mbrDF %>% 
  select(
    id = child_id, birth_year = faar
  )

mother_faarDF <- mbrDF %>% 
  select(
    id = !!sym(mother_id_column), birth_year = mor_faar
  ) %>% 
  mutate(
    birth_year = ifelse(birth_year == "<=1958", 1958, birth_year),
    birth_year = ifelse(birth_year == ">=1990", 1990, birth_year)
  )

father_faarDF <- mbrDF %>% 
  select(
    id = !!sym(father_id_column), birth_year = far_faar
  ) %>% 
  mutate(
    birth_year = ifelse(birth_year == "<=1944", 1944, birth_year),
    birth_year = ifelse(birth_year == "1945-1946", 1945.5, birth_year),
    birth_year = ifelse(birth_year == ">=1990", 1990, birth_year)
  )

faarDF <- rbind(child_faarDF, mother_faarDF, father_faarDF) %>% 
  filter(
    !is.na(birth_year) | !is.numeric(birth_year)
  ) %>% 
  mutate(
    birth_year = as.numeric(birth_year)
  ) %>% 
  group_by(
    id
  ) %>% 
  summarize(
    birth_year = mean(birth_year, na.rm = T),
    .groups = 'drop'
  )


if (nrow(faarDF) != length(unique(faarDF$id))) {
  
  stop("Duplicate sentrix id introduced during merging with mbr.")
  
}

phenoDF <- phenoDF %>% 
  left_join(
    faarDF,
    by = "id"
  )

# Sex sex using MBR

print(glue("{Sys.time()} - Assigning sex"))

sexDF <- mbrDF %>% 
  select(
    id = child_id, sex_mbr = kjonn
  ) %>% 
  mutate(
    sex = case_when(
      sex_mbr == "Mann" ~ 1,
      sex_mbr == "Kvinne" ~ 2
    )
  ) %>% 
  select(
    id, sex
  )

phenoDF <- phenoDF %>% 
  left_join(
    sexDF,
    by = "id"
  )

phenoDF$sex[phenoDF$role == "Mother"] <- 2
phenoDF$sex[phenoDF$role == "Father"] <- 1


if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
  
  stop("Duplicate sentrix id introduced during merging with mbr.")
  
}


# Load infection registry and match ids

print(glue("{Sys.time()} - Merging with infection registry"))

msisDF <- read.table(
  file = msis_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()

idMappingDF <- read.table(
  file = mother_msis_id_mapping_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
)

names(idMappingDF) <- c("id", kIdColumn)

msisDF <- msisDF %>% 
  left_join(
    idMappingDF,
    by = kIdColumn
  )

idMappingDF <- read.table(
  file = father_msis_id_mapping_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
)

names(idMappingDF) <- c("temp_id", kIdColumn)

msisDF <- msisDF %>% 
  left_join(
    idMappingDF,
    by = kIdColumn
  ) %>% 
  mutate(
    id = ifelse(!is.na(temp_id), temp_id, id)
  ) %>% 
  select(
    -temp_id
  )

idMappingDF <- read.table(
  file = child_msis_id_mapping_raw_table_path,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>%
  mutate(
    temp_id = paste0(!!sym(preg_id_column), "_", barn_nr)
  ) %>% 
  select(
    temp_id,
    !!kIdColumn := all_of(msisIdColumn)
  )

msisDF <- msisDF %>% 
  left_join(
    idMappingDF,
    by = kIdColumn
  ) %>% 
  mutate(
    id = ifelse(!is.na(temp_id), temp_id, id)
  ) %>% 
  select(
    -temp_id
  )

msisDF <- msisDF %>% 
  mutate(
    registered_date = as.Date(registrert_dato, "%d.%m.%Y"),
    hospitalized = as.numeric(factor(er_innlagt_sykehus, levels = c("Nei", "Ja"))) - 1
  ) %>% 
  select(
    id, registered_date, hospitalized
  ) %>% 
  filter(
    !is.na(id)
  ) %>% 
  group_by(
    id
  ) %>% 
  summarize(
    msis_first_registered = min(registered_date, na.rm = T), 
    msis_last_registered = max(registered_date, na.rm = T), 
    msis_hospitalized = max(hospitalized, na.rm = T), 
    .groups = 'drop'
  ) %>% 
  mutate(
    msis_hospitalized = ifelse(is.infinite(msis_hospitalized), NA, msis_hospitalized),
    hospitalized = ifelse(!is.na(msis_hospitalized) & msis_hospitalized == 1, 1, 0),
    infected = ifelse(!is.na(msis_first_registered), 1, 0)
  )

phenoDF <- phenoDF %>% 
  left_join(
    msisDF,
    by = "id"
  )


if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
  
  stop("Duplicate sentrix id introduced during processing of infection registry.")
  
}


# Get age at diagnosis

print(glue("{Sys.time()} - Getting age at diagnosis"))

phenoDF <- phenoDF %>% 
  mutate(
    age = 2021 - birth_year
  )

for (role in c("Child", "Mother", "Father")) {
  
  meanAge <- mean(phenoDF$age[phenoDF$role == role], na.rm = T)
  phenoDF$age[phenoDF$role == role & is.na(phenoDF$age)] <- meanAge
  
}


# Long covid phenotypes

print(glue("{Sys.time()} - Computing long covid phenotypes"))

longCovidPhenos <- list(
  reduced_smell_taste = "Reduced smell taste (kf120)", 
  brain_fog = "Brain Fog (kf480)", 
  poor_memory = "Poor Memory (kf481)", 
  dizziness = "Dizziness (kf479)", 
  heart_palpitation = "Heart Palpitation (kf476)", 
  fatigue = "Fatigue (kf468)", 
  headache = "Headache (kf484)", 
  skin_rash = "Skin rash (kf487)", 
  anxiety = "Anxiety (kf486)", 
  altered_smell_taste = "Altered smell taste (kf489)", 
  chest_pain = "Chest pain (kf475)", 
  shortness_breath = "Shortness breath (kf470)", 
  lung_function_reduced = "Lung function reduced (kf472)", 
  cough = "Cough (kf471)"
)
pheno_variable_to_question[["reduced_smell_taste"]] <- unique(c(pheno_variable_to_question[["reduced_smell_taste"]], "kf120"))
pheno_variable_to_question[["brain_fog"]] <- unique(c(pheno_variable_to_question[["brain_fog"]], "kf480"))
pheno_variable_to_question[["poor_memory"]] <- unique(c(pheno_variable_to_question[["poor_memory"]], "kf481"))
pheno_variable_to_question[["dizziness"]] <- unique(c(pheno_variable_to_question[["dizziness"]], "kf479"))
pheno_variable_to_question[["heart_palpitation"]] <- unique(c(pheno_variable_to_question[["heart_palpitation"]], "kf476"))
pheno_variable_to_question[["fatigue"]] <- unique(c(pheno_variable_to_question[["fatigue"]], "kf468"))
pheno_variable_to_question[["headache"]] <- unique(c(pheno_variable_to_question[["headache"]], "kf484"))
pheno_variable_to_question[["skin_rash"]] <- unique(c(pheno_variable_to_question[["skin_rash"]], "kf487"))
pheno_variable_to_question[["anxiety"]] <- unique(c(pheno_variable_to_question[["anxiety"]], "kf486"))
pheno_variable_to_question[["altered_smell_taste"]] <- unique(c(pheno_variable_to_question[["altered_smell_taste"]], "kf489"))
pheno_variable_to_question[["chest_pain"]] <- unique(c(pheno_variable_to_question[["chest_pain"]], "kf475"))
pheno_variable_to_question[["shortness_breath"]] <- unique(c(pheno_variable_to_question[["shortness_breath"]], "kf470"))
pheno_variable_to_question[["lung_function_reduced"]] <- unique(c(pheno_variable_to_question[["cough"]], "kf472"))
pheno_variable_to_question[["cough"]] <- unique(c(pheno_variable_to_question[["reduced_smell_taste"]], "kf471"))

longCovidFactorWeights <- data.frame(
  brain_fog = c(0.952392278245466, -0.108327867721912), 
  poor_memory = c(0.815151856835235, 0.0959893863922384), 
  dizziness = c(0.807334011072063, -0.105211767181289), 
  heart_palpitation = c(0.739879182820943, 0.0507881302224967), 
  fatigue = c(0.622492944472992, 0.301909054576685), 
  headache = c(0.572262858913896, 0.181154869329663), 
  skin_rash = c(0.449545853704445, 0.0248297029917953), 
  anxiety = c(0.412869730153714, 0.152068325497297), 
  altered_smell_taste = c(0.335692445248068, 0.0693279720080806), 
  chest_pain = c(0.208312587760769, 0.451259084957272), 
  shortness_breath = c(0.0189268663688903, 0.964408601577604), 
  lung_function_reduced = c(0.0158479449038898, 0.576540782363923), 
  cough = c(0.0457152256335479, 0.60115307666968)
)

phenoDF$anySymptomShort <- 0 
phenoDF$anySymptomLong <- 0 
phenoDF$long_covid_factor1 <- 0
phenoDF$long_covid_factor2 <- 0

long_covid_questions <- c()

for (longCovidPheno in names(longCovidPhenos)) {
  
  long_covid_questions <- c(long_covid_questions, pheno_variable_to_question[[longCovidPheno]])
  
  lastPheno <- paste0(longCovidPheno, "_last_reported")
  
  phenoDF <- phenoDF %>% 
    mutate(
      last_pheno_date = as.Date(!!sym(lastPheno)/86400, origin = "1582-10-14"),
      symptome_duration = as.numeric(last_pheno_date - msis_last_registered)
    )
  
  phenoShort <- paste0(longCovidPheno, "_short")
  phenoLong <- paste0(longCovidPheno, "_long")
  phenoDF[[phenoShort]] <- ifelse(phenoDF$symptome_duration < 90, 1, 0)
  phenoDF[[phenoLong]] <- ifelse(phenoDF$symptome_duration >= 90, 1, 0)
  
  pheno_variable_to_question[[phenoShort]] <- unique(c(pheno_variable_to_question[[phenoShort]], pheno_variable_to_question[[longCovidPheno]]))
  pheno_variable_to_question[[phenoLong]] <- unique(c(pheno_variable_to_question[[phenoLong]], pheno_variable_to_question[[longCovidPheno]]))
  
  phenoDF$anySymptomShort[!is.na(phenoDF$msis_last_registered) & !is.na(phenoDF[[longCovidPheno]]) & phenoDF[[longCovidPheno]] == 1 & phenoDF$symptome_duration < 90] <- 1
  
  symptom_long <- !is.na(phenoDF$msis_last_registered) & !is.na(phenoDF[[longCovidPheno]]) & phenoDF[[longCovidPheno]] == 1 & phenoDF$symptome_duration >= 90
  phenoDF$anySymptomLong[symptom_long] <- 1
  
  if (longCovidPheno %in% names(longCovidFactorWeights)) {
    
    phenoDF$long_covid_factor1[symptom_long] <- phenoDF$long_covid_factor1[symptom_long] + longCovidFactorWeights[[longCovidPheno]][1]
    phenoDF$long_covid_factor2[symptom_long] <- phenoDF$long_covid_factor2[symptom_long] + longCovidFactorWeights[[longCovidPheno]][2]
    
  }
  
  phenoDF <- phenoDF %>% 
    select(
      -last_pheno_date, -symptome_duration
    )
}

long_covid_questions <- unique(long_covid_questions)

pheno_variable_to_question[["anySymptomShort"]] <- long_covid_questions
pheno_variable_to_question[["anySymptomLong"]] <- long_covid_questions
pheno_variable_to_question[["long_covid_factor1"]] <- long_covid_questions
pheno_variable_to_question[["long_covid_factor2"]] <- long_covid_questions

phenoDF <- phenoDF %>% 
  mutate(
    long_covid_vs_recovered = ifelse(!is.na(msis_last_registered), 0, NA),
    long_covid_vs_recovered = ifelse(anySymptomLong == 1, 1, long_covid_vs_recovered),
    long_covid_vs_population = ifelse(anySymptomLong == 1, 1, 0),
    long_covid_factor1_vs_recovered = ifelse(!is.na(long_covid_vs_recovered) & long_covid_vs_recovered == 1, long_covid_factor1, long_covid_vs_recovered),
    long_covid_factor2_vs_recovered = ifelse(!is.na(long_covid_vs_recovered) & long_covid_vs_recovered == 1, long_covid_factor2, long_covid_vs_recovered)
  )
pheno_variable_to_question[["long_covid_vs_recovered"]] <- long_covid_questions
pheno_variable_to_question[["long_covid_vs_population"]] <- long_covid_questions
pheno_variable_to_question[["long_covid_factor1_vs_recovered"]] <- long_covid_questions
pheno_variable_to_question[["long_covid_factor2_vs_recovered"]] <- long_covid_questions


# Variables mapping

variablesDF <- data.frame(
  name = character(length(pheno_variable_to_question)),
  question = character(length(pheno_variable_to_question)),
  stringsAsFactors = F
)

for (i in 1:length(pheno_variable_to_question)) {
  
  variable <- names(pheno_variable_to_question)[i]
  questions <- pheno_variable_to_question[[variable]]
  questionsString <- paste(questions, collapse = ",")
  
  variablesDF$name[i] <- variable
  variablesDF$question[i] <- questionsString
  
}


# Save to file

print(glue("{Sys.time()} - Saving to file"))

write.table(
  x = phenoDF,
  file = gzfile(covidTable),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

write.table(
  x = variablesDF,
  file = gzfile(file.path(docsFolder, "variable_to_question.gz")),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)


# Write long covid docs

print(glue("{Sys.time()} - Writing long covid documentation"))

writeLongCovidDocs(
  longCovidDocsFolder = longCovidDocsFolder,
  longCovidDocsFile = longCovidDocsFile,
  phenoDF = phenoDF
)

