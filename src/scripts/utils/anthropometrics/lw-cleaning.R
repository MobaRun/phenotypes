#!/usr/bin/env Rscript

##
#
# This script processes the length and weight values for the children. 
#
##


## Command line input

# Command line arguments
args <- commandArgs(TRUE)

qcFolder <- args[1]
project_number <- args[2] 


##
#
# Debug Marc - do not uncomment
# args to run locally on Hunt 
# 
qcFolder <- "/mnt/work/marc/pheno_22-09-19/qc_tmp"
project_number <- 2824
#
##


## Libraries import

libFolder <- "~/R/R_4.1"

library(stringr, lib = libFolder)
library(crayon, lib = libFolder)
library(dplyr, lib = libFolder)
library(janitor, lib = libFolder)
library(purrr, lib = libFolder)
library(glue, lib = libFolder)



## Parameters


# Minimal number of missing values allowed
minValidValuesAfter2 <- 2
minValidValuesBefore2 <- 3
minValidValuesBeforeFirstImputed <- 2
minValidValuesAfterLastImputed <- 2


# Cache for the standard values of the normal distribution
anchorUp <- pnorm(1)
anchorCenter <- pnorm(0)
anchorDown <- pnorm(-1)
margin99 <- qnorm(0.99)


# Maximal number of iterations
maxIt <- 100

# If true, modified growth profiles will be plotted
exportProfiles <- F


# External functions and variable mapping

source("src/scripts/utils/anthropometrics/lw-cleaning-functions.R")

source("src/scripts/utils/anthropometrics/variables_mapping.R")


# Housekeeping

lwQcLogFolder <- file.path(qcFolder, "lw_log")

if (!file.exists(lwQcLogFolder)) {
  
  dir.create(
    path = lwQcLogFolder,
    showWarnings = T,
    recursive = T
  )
  
}

lwQcPlotsFolder <- file.path(lwQcLogFolder, "plots")

if (!file.exists(lwQcPlotsFolder)) {
  
  dir.create(
    path = lwQcPlotsFolder,
    showWarnings = T,
    recursive = T
  )
  
}

if (exportProfiles) {
  
  curvesFolder <- file.path(lwQcPlotsFolder, "curves")
  
  if (!dir.exists(curvesFolder)) {
    
    dir.create(
      path = curvesFolder,
      showWarnings = T,
      recursive = T
    )
    
  }
}

# Load values from the table

print(paste(Sys.time(), " Loading length and weight values"))

values <- read.table(
  file = file.path(qcFolder, "child_anthropometrics_raw.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

originalValues <- values
values$log <- ""

# Run cleaning

values <- iterativeCleaning(values)


# Store the log

print(paste(Sys.time(), " Log of the modifications"))

dummyIds <- sample(nrow(values), nrow(values))
logBridgeDF <- values %>%
  select(
    child_id, log
  ) %>%
  mutate(
    dummyId = dummyIds
  )
bridgeDF <- logBridgeDF %>%
  select(
    child_id, dummyId
  )
logDF <- logBridgeDF %>%
  select(
    dummyId, log
  )

write.table(
  x = logDF, 
  file = file.path(lwQcLogFolder, "cleaning_log"), 
  col.names = T, 
  row.names = F, 
  quote = F, 
  sep = "\t"
)
write.table(
  x = bridgeDF, 
  file = file.path(lwQcLogFolder, "log_bridge"), 
  col.names = T, 
  row.names = F, 
  quote = F, 
  sep = "\t"
)


# Plot profiles

print(paste(Sys.time(), " Plotting profiles"))

if (exportProfiles) {
  
  exportCurves(
    originalValues = originalValues,
    values = values,
    bridgeDF = bridgeDF,
    qcFolderLocal = curvesFolder
  )
  
}


# Plot length-weight

print(paste(Sys.time(), " Plotting length-weight"))

exportLenghtWeight(
  originalValues = originalValues,
  values = values,
  qcFolderLocal = lwQcPlotsFolder
)


# Store the number of corrected values

print(paste(Sys.time(), " Storing the number of modified values"))

qcStepV <- character()
correctionV <- character()
phenoV <- character()
nV <- numeric()
k <- 1

for (variables in list(length_columns, weight_columns)) {
  
  for (phenoName in variables) {
    
    diff <- rep("none", nrow(values))
    
    is <- is.na(originalValues[[phenoName]]) & is.na(values[[phenoName]])
    diff[is] <- "NA"
    
    is <- !is.na(originalValues[[phenoName]]) & is.na(values[[phenoName]])
    diff[is] <- "outlier"
    
    is <- is.na(originalValues[[phenoName]]) & !is.na(values[[phenoName]])
    diff[is] <- "imputed"
    
    is <- !is.na(originalValues[[phenoName]]) & !is.na(values[[phenoName]]) & originalValues[[phenoName]] != values[[phenoName]]
    diff[is] <- "corrected"
    
    table <- as.data.frame(table(diff))
    
    print(
      paste0(phenoName, ":")
    )
    print(table)
    
    for (i in 1:nrow(table)) {
      
      correction <- as.character(table[i, "diff"])
      n <- table[i, "Freq"]
      
      qcStepV[k] <- "imputative-cleaning"
      correctionV[k] <- correction
      phenoV[k] <- phenoName
      nV[k] <- n
      
      k <- k + 1
      
    }
  }
}

imputativeQcDF <- data.frame(
  qcStep = qcStepV,
  correction = correctionV,
  pheno = phenoV,
  n = nV,
  stringsAsFactors = F
)

qcFile <- file.path(qcFolder, "qcDF")
qcDF <- read.table(qcFile, header = T, stringsAsFactors = F)
qcDF <- rbind(qcDF, imputativeQcDF)
write.table(qcDF, file = qcFile, col.names = T, row.names = F, quote = F)


# Save phenotypes to a new table

print(paste(Sys.time(), "Saving the results."))

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
  file = gzfile(file.path(qcFolder, "child_anthropometrics.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)
