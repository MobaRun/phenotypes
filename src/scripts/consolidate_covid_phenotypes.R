##
#
# This script processes phenotypes and exports files for genetic analyses and associated documentation.
# 
##


print(paste0(Sys.time(), " - Organize phenotypes in a clean table and write documentation"))

# Libraries

libFolder <- "~/R/R_4.1"

# This randomly fails the first time but not the second, seems to be an error in the latest R versions

loadLibraries <- function() {
    
    library(isoband, lib = libFolder)
    library(farver, lib = libFolder)
    library(ellipsis, lib = libFolder)
    library(scales, lib = libFolder)
    library(backports, lib = libFolder)
    library(vctrs, lib = libFolder)
    library(crayon, lib = libFolder)
    library(tidyr, lib = libFolder)
    library(dplyr, lib = libFolder)
    library(MASS, lib = libFolder)
    library(gamlss.data, lib = libFolder)
    library(gamlss.dist, lib = libFolder)
    library(nlme, lib = libFolder)
    library(gamlss, lib = libFolder)
    library(withr, lib = libFolder)
    library(labeling, lib = libFolder)
    library(digest, lib = libFolder)
    library(reshape2, lib = libFolder)
    library(ggplot2, lib = libFolder)
    library(grid, lib = libFolder)
    library(scico, lib = libFolder)
    library(gtable, lib = libFolder)
    #library(conflicted, lib = libFolder)
    library(stringr, lib = libFolder)
    library(jsonlite, lib = libFolder)
    library(glue, lib = libFolder)
    library(janitor, lib = libFolder)
    
}

tryCatch(
    {
        loadLibraries()
        
    }, error = function(error_condition) {
        
        loadLibraries()
        
    }
)


# Solve namespace conflicts

# conflict_prefer("select", "dplyr")
# conflict_prefer("filter", "dplyr")
select <- dplyr::select
filter <- dplyr::filter


# Paths from command line
args <- commandArgs(TRUE)

linkageFolder <- args[1]
quesFolder <- args[2]
rawTablesFolder <- args[3]
covidTable <- args[4]
docsFolder <- args[5]


# Functions

source("src/scripts/utils/covid/long_covid_docs.R")


# Parameters

batchOrder <- c("NORMENT_JAN21", "NORMENT_MAI21", "NORMENT_SEP20_R996R1029", "NORMENT-JAN20", "NORMENT-FEB18", "NORMENT-JUN15", "NORMENT-MAY16", "NORMENT-JAN15", "ROTTERDAM1", "ROTTERDAM2", "HARVEST", "TED", "PDB1382_R875_R876")

theme_set(theme_bw(base_size = 16))

updateDocs <- T


# Paths

docsFile <- file.path(docsFolder, "phenos.md")

longCovidDocsFolder <- file.path(docsFolder, "long_covid")
longCovidDocsFile <- file.path(longCovidDocsFolder, "long_covid.md")


# Housekeeping

if (!dir.exists(longCovidDocsFolder)) {
  
  dir.create(
    path = longCovidDocsFolder,
    showWarnings = T,
    recursive = T
  )
  
}


# Load identifiers

print(paste0(Sys.time(), "    Loading identifiers"))

childIdDF <- read.table(
  file = file.path(rawTablesFolder, "20220516_MoBaGeneticsTot_Child_PDB2824.gz"),
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  mutate(
    child_id = paste0(preg_id_2824, "_", barn_nr)
  )

for (colname in c("sentrix_id", "role", "batch", "sampletype")) {
  
  childIdDF[[colname]] <- str_replace_all(
    string = childIdDF[[colname]], 
    pattern = " ",
    repl = ""
  )
  
}

motherIdDF <- read.table(
  file = file.path(rawTablesFolder, "20220516_MoBaGeneticsTot_Mother_PDB2824.gz"),
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()


for (colname in c("m_id_2824", "sentrix_id", "role", "batch", "sampletype")) {
  
  motherIdDF[[colname]] <- str_replace_all(
    string = motherIdDF[[colname]], 
    pattern = " ",
    repl = ""
  )
  
}


fatherIdDF <- read.table(
  file = file.path(rawTablesFolder, "20220516_MoBaGeneticsTot_Father_PDB2824.gz"),
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names()


for (colname in c("f_id_2824", "sentrix_id", "role", "batch", "sampletype")) {
  
  fatherIdDF[[colname]] <- str_replace_all(
    string = fatherIdDF[[colname]], 
    pattern = " ",
    repl = ""
  )
  
}

idDF <- rbind(
  childIdDF %>% 
    select(
      id = child_id, sentrix_id, role, batch
    ),
  motherIdDF %>% 
    select(
      id = m_id_2824, sentrix_id, role, batch
    ),
  fatherIdDF %>% 
    select(
      id = f_id_2824, sentrix_id, role, batch
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


# Set up a data frame for phenotypes

print(glue("{Sys.time()} - Setting up pheno table"))

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
        vaccine_arm_pain = 0,
        vaccine_feber = 0,
        vaccine_freezing = 0,
        vaccine_feeling_unwell = 0,
        vaccine_bad_appetite = 0,
        vaccine_headache = 0,
        vaccine_pain_other_place_than_injection = 0,
        vaccine_skin_bleeding_ecchymosis = 0,
        vaccine_nose_bleeding = 0,
        vaccine_gum_bleeding = 0,
        vaccine_mouth_ulcer = 0,
        vaccine_thrombus = 0,
        vaccine_unusually_strong_menstruation = 0,
        vaccine_nausea = 0,
        vaccine_belly_pain = 0,
        vaccine_diarrhea = 0,
        vaccine_dizziness = 0,
        vaccine_faintness = 0
    )


# long covid pheno

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
                
                if (file.exists(questionnaireFile)) {
                    
                    print(glue("{Sys.time()} - Loading phenotypes from {questionnaireFile}"))
                    
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
                                yes = paste0(preg_id_2824, "_", barn_nr),
                                no = p_id_2824
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
                            
                        }
                        
                        if (nrow(phenoDF) != length(unique(phenoDF$sentrix_id))) {
                            
                            stop("Duplicate sentrix id introduced during processing of 'kf166'.")
                            
                        }
                    }
                    if ("kf120" %in% names(quesDF)) {
                      
                      print(summary(quesDF$kf120))
                      
                      stop("DEBUG")
                        
                        quesDF$kf120 <- as.numeric(factor(quesDF$kf120, levels = c("NEI", "JA"))) - 1
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf120",
                            phenoName = "reduced_smell_taste",
                            folder = folder
                        )
                    }
                    if ("kf480" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf480",
                            phenoName = "brain_fog",
                            folder = folder
                        )
                    }
                    if ("kf481" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf481",
                            phenoName = "poor_memory",
                            folder = folder
                        )
                    }
                    if ("kf479" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf479",
                            phenoName = "dizziness",
                            folder = folder
                        )
                    }
                    if ("kf476" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf476",
                            phenoName = "heart_palpitation",
                            folder = folder
                        )
                    }
                    if ("kf468" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf468",
                            phenoName = "fatigue",
                            folder = folder
                        )
                    }
                    if ("kf484" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf484",
                            phenoName = "headache",
                            folder = folder
                        )
                    }
                    if ("kf487" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf487",
                            phenoName = "skin_rash",
                            folder = folder
                        )
                    }
                    if ("kf486" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf486",
                            phenoName = "anxiety",
                            folder = folder
                        )
                    }
                    if ("kf489" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf489",
                            phenoName = "altered_smell_taste",
                            folder = folder
                        )
                    }
                    if ("kf475" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf475",
                            phenoName = "chest_pain",
                            folder = folder
                        )
                    }
                    if ("kf470" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf470",
                            phenoName = "shortness_breath",
                            folder = folder
                        )
                    }
                    if ("kf472" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf472",
                            phenoName = "lung_function_reduced",
                            folder = folder
                        )
                    }
                    if ("kf471" %in% names(quesDF)) {
                        
                        phenoDF <- longCovidPhenoImport(
                            quesDF = quesDF,
                            phenoDF = phenoDF,
                            quesNumber = "kf471",
                            phenoName = "cough",
                            folder = folder
                        )
                    }
                    if ("kf1081" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_arm_pain = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1081) & quesDF$kf1081 == "JA"], 1, vaccine_arm_pain)
                            )
                    }
                    if ("kf1126" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_arm_pain = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1126) & quesDF$kf1126 == "JA"], 1, vaccine_arm_pain)
                            )
                    }
                    if ("kf1082" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_feber = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1082) & quesDF$kf1082 == "JA"], 1, vaccine_feber)
                            )
                    }
                    if ("kf1127" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_feber = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1127) & quesDF$kf1127 == "JA"], 1, vaccine_feber)
                            )
                    }
                    if ("kf1083" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_freezing = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1083) & quesDF$kf1083 == "JA"], 1, vaccine_freezing)
                            )
                    }
                    if ("kf1128" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_freezing = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1128) & quesDF$kf1128 == "JA"], 1, vaccine_freezing)
                            )
                    }
                    if ("kf1084" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_feeling_unwell = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1084) & quesDF$kf1084 == "JA"], 1, vaccine_feeling_unwell)
                            )
                    }
                    if ("kf1129" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_feeling_unwell = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1129) & quesDF$kf1129 == "JA"], 1, vaccine_feeling_unwell)
                            )
                    }
                    if ("kf1085" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_bad_appetite = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1085) & quesDF$kf1085 == "JA"], 1, vaccine_bad_appetite)
                            )
                    }
                    if ("kf1130" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_bad_appetite = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1130) & quesDF$kf1130 == "JA"], 1, vaccine_bad_appetite)
                            )
                    }
                    if ("kf1086" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_headache = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1086) & quesDF$kf1086 == "JA"], 1, vaccine_headache)
                            )
                    }
                    if ("kf1131" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_headache = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1131) & quesDF$kf1131 == "JA"], 1, vaccine_headache)
                            )
                    }
                    if ("kf1087" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_pain_other_place_than_injection = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1087) & quesDF$kf1087 == "JA"], 1, vaccine_pain_other_place_than_injection)
                            )
                    }
                    if ("kf1132" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_pain_other_place_than_injection = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1132) & quesDF$kf1132 == "JA"], 1, vaccine_pain_other_place_than_injection)
                            )
                    }
                    if ("kf1088" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_skin_bleeding_ecchymosis = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1088) & quesDF$kf1088 == "JA"], 1, vaccine_skin_bleeding_ecchymosis)
                            )
                    }
                    if ("kf1133" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_skin_bleeding_ecchymosis = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1133) & quesDF$kf1133 == "JA"], 1, vaccine_skin_bleeding_ecchymosis)
                            )
                    }
                    if ("kf1089" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_nose_bleeding = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1089) & quesDF$kf1089 == "JA"], 1, vaccine_nose_bleeding)
                            )
                    }
                    if ("kf1134" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_nose_bleeding = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1134) & quesDF$kf1134 == "JA"], 1, vaccine_nose_bleeding)
                            )
                    }
                    if ("kf1090" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_gum_bleeding = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1090) & quesDF$kf1090 == "JA"], 1, vaccine_nose_bleeding)
                            )
                    }
                    if ("kf1135" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_gum_bleeding = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1135) & quesDF$kf1135 == "JA"], 1, vaccine_nose_bleeding)
                            )
                    }
                    if ("kf1091" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_mouth_ulcer = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1091) & quesDF$kf1091 == "JA"], 1, vaccine_mouth_ulcer)
                            )
                    }
                    if ("kf1136" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_mouth_ulcer = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1136) & quesDF$kf1136 == "JA"], 1, vaccine_mouth_ulcer)
                            )
                    }
                    if ("kf1092" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_thrombus = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1092) & quesDF$kf1092 == "JA"], 1, vaccine_thrombus)
                            )
                    }
                    if ("kf1137" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_thrombus = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1137) & quesDF$kf1137 == "JA"], 1, vaccine_thrombus)
                            )
                    }
                    if ("kf1093" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_unusually_strong_menstruation = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1093) & quesDF$kf1093 == "JA"], 1, vaccine_unusually_strong_menstruation)
                            )
                    }
                    if ("kf1138" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_unusually_strong_menstruation = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1138) & quesDF$kf1138 == "JA"], 1, vaccine_unusually_strong_menstruation)
                            )
                    }
                    if ("kf1095" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_nausea = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1095) & quesDF$kf1095 == "JA"], 1, vaccine_nausea)
                            )
                    }
                    if ("kf1140" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_nausea = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1140) & quesDF$kf1140 == "JA"], 1, vaccine_nausea)
                            )
                    }
                    if ("kf1096" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_belly_pain = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1096) & quesDF$kf1096 == "JA"], 1, vaccine_belly_pain)
                            )
                    }
                    if ("kf1141" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_belly_pain = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1141) & quesDF$kf1141 == "JA"], 1, vaccine_belly_pain)
                            )
                    }
                    if ("kf1097" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_diarrhea = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1097) & quesDF$kf1097 == "JA"], 1, vaccine_diarrhea)
                            )
                    }
                    if ("kf1142" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_diarrhea = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1142) & quesDF$kf1142 == "JA"], 1, vaccine_diarrhea)
                            )
                    }
                    if ("kf1098" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_dizziness = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1098) & quesDF$kf1098 == "JA"], 1, vaccine_dizziness)
                            )
                    }
                    if ("kf1143" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_dizziness = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1143) & quesDF$kf1143 == "JA"], 1, vaccine_dizziness)
                            )
                    }
                    if ("kf1099" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_faintness = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1099) & quesDF$kf1099 == "JA"], 1, vaccine_faintness)
                            )
                    }
                    if ("kf1144" %in% names(quesDF)) {
                        
                        phenoDF <- phenoDF %>% 
                            mutate(
                                vaccine_faintness = ifelse(id %in% quesDF$id[!is.na(quesDF$kf1144) & quesDF$kf1144 == "JA"], 1, vaccine_faintness)
                            )
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
    file = file.path(rawTablesFolder, "PDB2824_MFR_541_v12.gz"),
    sep = "\t",
    header = T,
    quote = "",
    stringsAsFactors = F,
    comment.char = ""
) %>% 
    clean_names() %>% 
    mutate(
        child_id = paste0(preg_id_2824, "_", barn_nr)
    )

child_faarDF <- mbrDF %>% 
    select(
        id = child_id, birth_year = faar
    )

mother_faarDF <- mbrDF %>% 
    select(
        id = m_id_2824, birth_year = mor_faar
    ) %>% 
    mutate(
        birth_year = ifelse(birth_year == "<=1958", 1958, birth_year),
        birth_year = ifelse(birth_year == ">=1990", 1990, birth_year)
    )

father_faarDF <- mbrDF %>% 
    select(
        id = f_id_2824, birth_year = far_faar
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
    file = file.path(rawTablesFolder, "PDB2824_MSIS-data_MoBa.gz"),
    sep = "\t",
    header = T,
    quote = "",
    stringsAsFactors = F,
    comment.char = ""
) %>% 
    clean_names()

idMappingFile <- file.path(rawTablesFolder, glue("Mor_ID_2824_2021_11_17sav.gz"))

idMappingDF <- read.table(
    file = idMappingFile,
    sep = "\t",
    header = T,
    quote = "",
    stringsAsFactors = F,
    comment.char = ""
)

names(idMappingDF) <- c("id", "pid_k_2824")

msisDF <- msisDF %>% 
    left_join(
        idMappingDF,
        by = "pid_k_2824"
    )

idMappingFile <- file.path(rawTablesFolder, glue("Far_ID_2824_2021_11_17sav.gz"))

idMappingDF <- read.table(
    file = idMappingFile,
    sep = "\t",
    header = T,
    quote = "",
    stringsAsFactors = F,
    comment.char = ""
)

names(idMappingDF) <- c("temp_id", "pid_k_2824")

msisDF <- msisDF %>% 
    left_join(
        idMappingDF,
        by = "pid_k_2824"
    ) %>% 
    mutate(
        id = ifelse(!is.na(temp_id), temp_id, id)
    ) %>% 
    select(
        -temp_id
    )

idMappingFile <- file.path(rawTablesFolder, glue("Barn_ID_2824_2021_11_17sav.gz"))

idMappingDF <- read.table(
    file = idMappingFile,
    sep = "\t",
    header = T,
    quote = "",
    stringsAsFactors = F,
    comment.char = ""
) %>% 
    clean_names() %>%
    mutate(
        temp_id = paste0(preg_id_2824, "_", barn_nr)
    ) %>% 
    select(
        temp_id,
        pid_k_2824 = pid_msis_2824
    )

msisDF <- msisDF %>% 
    left_join(
        idMappingDF,
        by = "pid_k_2824"
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


for (longCovidPheno in names(longCovidPhenos)) {
  
  lastPheno <- paste0(longCovidPheno, "_last_reported")
  
  phenoDF <- phenoDF %>% 
    mutate(
      last_pheno_date = as.Date(!!sym(lastPheno)/86400, origin = "1582-10-14"),
      symptome_duration = as.numeric(last_pheno_date - msis_last_registered)
    )
  
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

phenoDF <- phenoDF %>% 
    mutate(
        long_covid_vs_recovered = ifelse(!is.na(msis_last_registered), 0, NA),
        long_covid_vs_recovered = ifelse(anySymptomLong == 1, 1, long_covid_vs_recovered),
        long_covid_vs_population = ifelse(anySymptomLong == 1, 1, 0),
        long_covid_factor1_vs_recovered = ifelse(!is.na(long_covid_vs_recovered) & long_covid_vs_recovered == 1, long_covid_factor1, long_covid_vs_recovered),
        long_covid_factor2_vs_recovered = ifelse(!is.na(long_covid_vs_recovered) & long_covid_vs_recovered == 1, long_covid_factor2, long_covid_vs_recovered)
    )


# Write long covid docs

print(glue("{Sys.time()} - Saving long covid docs"))

writeLongCovidDocs(
    longCovidDocsFolder = longCovidDocsFolder,
    longCovidDocsFile = longCovidDocsFile,
    phenoDF = phenoDF
)


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
