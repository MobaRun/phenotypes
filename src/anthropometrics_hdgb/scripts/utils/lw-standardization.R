#!/usr/bin/env Rscript

##
#
# This script processes the length and weight values for the children. 
#
##


## Command line input

# Command line arguments
args <- commandArgs(TRUE)

tablesFolder <- args[1]
project_number <- args[2]
docs_folder <- args[3]


##
#
# Debug Marc - do not uncomment
# args to run locally on Hunt 
# extern
# tablesFolder <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-05-07"
# project_number <- 2824
# docs_folder <- "docs/anthropometrics/24-05-07/standardization"
#
##


## Libraries import

library(cli)
library(janitor)
library(stringr)
library(dplyr)
library(gamlss)
library(glue)
library(ggplot2)
library(gtable)
library(grid)



## Parameters


# External functions

source("src/anthropometrics_hdgb/scripts/utils/standardization_functions.R")
source("src/anthropometrics_hdgb/scripts/utils/standardization_docs_functions.R")
source("src/anthropometrics_hdgb/scripts/utils/variables_mapping.R")


# Currently exclude values after 8y

length_columns <- length_columns[1:length(length_columns)]
weight_columns <- weight_columns[1:length(weight_columns)]
bmi_columns <- bmi_columns[1:length(bmi_columns)]


# Housekeeping

docs_file <- file.path(docs_folder, "standardization.md")



# Load values from the table

print(paste(Sys.time(), " Loading length and weight values"))

values <- read.table(
  file = file.path(tablesFolder, "child_anthropometrics.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

pregnancy_values <- read.table(
  file = file.path(tablesFolder, "pregnancy.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  select(
    child_id, pregnancy_duration
  )

values <- values %>% 
  left_join(
    pregnancy_values,
    by = "child_id"
  )

# Compute BMI

print(paste(Sys.time(), "Computing BMI."))

for (i in 1:length(length_columns)) {
  
  length_name <- length_columns[i]
  weight_name <- weight_columns[i]
  bmi_name <- bmi_columns[i]
  
  values <- values %>% 
    mutate(
      !!bmi_name := ifelse(!is.na(!!sym(length_name)) & !is.na(!!sym(weight_name)), 10000 * !!sym(weight_name) / (!!sym(length_name) * !!sym(length_name)), NA)
    )
  
}


# Set up docs

print(paste0(Sys.time(), "    Setting up docs"))

n_children <- length(unique(values$child_sentrix_id[!is.na(values$child_sentrix_id)]))
n_mother <- length(unique(values$mother_sentrix_id[!is.na(values$mother_sentrix_id)]))
n_father <- length(unique(values$father_sentrix_id[!is.na(values$father_sentrix_id)]))
n_pregnancies <- sum(!is.na(values$child_sentrix_id) | !is.na(values$mother_sentrix_id) | !is.na(values$father_sentrix_id))

write(x = "# Phenotypes\n", file = docs_file, append = F)
write(x = glue("{nrow(values)} pregnancies, {n_pregnancies} with genotypes: {n_children} children, {n_mother} mothers, {n_father} fathers\n\n"), file = docs_file, append = T)


# Covariates

median_duration <- median(values$pregnancy_duration, na.rm = T)
values$pregnancy_duration_1 <- ifelse(is.na(values$pregnancy_duration), median_duration, values$pregnancy_duration) 


# Length

for (phenoI in 1:length(length_columns)) {
  
  phenoName <- length_columns[phenoI]
  
  print(paste0(Sys.time(), "    Standardizing ", phenoName))
  
  timePoint <- strsplit(phenoName, "_")[[1]][2]
  
  zPhenoName <- paste0("z_", phenoName)
  phenoLabel <- paste0("Length/height at ", timePoint)
  phenoUnit <- "cm"
  familyName <- "NO"
  family <- NO
  
  sigmaFormula <- paste0(" ~ pregnancy_duration_1")
  
  formula <- paste0(phenoName, " ~ fp(pregnancy_duration_1)")
  formulaLinear <- paste0(phenoName, " ~ pregnancy_duration_1")
  
  refvalues <- values %>% 
    filter(
      unrelated_children == 1
    )
  
  zDF <- standardizeBySex(
    trainingvalues = refvalues,
    values = values,
    id = "child_id",
    x = "pregnancy_duration_1",
    y = phenoName,
    zY = zPhenoName,
    formula = as.formula(formula),
    sigmaFormula = as.formula(sigmaFormula),
    formula2 = as.formula(formulaLinear),
    sigmaFormula2 = as.formula(sigmaFormula),
    family = family
  )
  
  values <- values %>%
    left_join(
      zDF,
      by = "child_id"
    )
  
  writeDocs(
    phenoName = phenoName,
    zPhenoName = zPhenoName,
    phenoLabel = phenoLabel,
    phenoUnit = phenoUnit,
    formula = formula,
    sigmaFormula = sigmaFormula,
    family = familyName,
    values = values,
    docsFolder = docs_folder,
    docsFile = docs_file
  )
}


# Weight

for (phenoI in 1:length(weight_columns)) {
  
  phenoName <- weight_columns[phenoI]
  
  print(paste0(Sys.time(), "    Standardizing ", phenoName))
  
  timePoint <- strsplit(phenoName, "_")[[1]][2]
  
  zPhenoName <- paste0("z_", phenoName)
  phenoLabel <- paste0("Weight at ", timePoint)
  phenoUnit <- "kg"
  familyName <- "NO"
  family <- NO
  
  sigmaFormula <- paste0(" ~ pregnancy_duration_1")
  
  formula <- paste0(phenoName, " ~ fp(pregnancy_duration_1)")
  formulaLinear <- paste0(phenoName, " ~ pregnancy_duration_1")
  
  refvalues <- values %>% 
    filter(
      unrelated_children == 1
    )
  
  zDF <- standardizeBySex(
    trainingvalues = refvalues,
    values = values,
    id = "child_id",
    x = "pregnancy_duration_1",
    y = phenoName,
    zY = zPhenoName,
    formula = as.formula(formula),
    sigmaFormula = as.formula(sigmaFormula),
    formula2 = as.formula(formulaLinear),
    sigmaFormula2 = as.formula(sigmaFormula),
    family = family
  )
  
  values <- values %>%
    left_join(
      zDF,
      by = "child_id"
    )
  
  writeDocs(
    phenoName = phenoName,
    zPhenoName = zPhenoName,
    phenoLabel = phenoLabel,
    phenoUnit = phenoUnit,
    formula = formula,
    sigmaFormula = sigmaFormula,
    family = familyName,
    values = values,
    docsFolder = docs_folder,
    docsFile = docs_file
  )
}


# BMI

for (phenoI in 1:length(bmi_columns)) {
  
  phenoName <- bmi_columns[phenoI]
  
  print(paste0(Sys.time(), "    Standardizing ", phenoName))
  
  timePoint <- strsplit(phenoName, "_")[[1]][2]
  
  zPhenoName <- paste0("z_", phenoName)
  phenoLabel <- paste0("BMI at ", timePoint)
  phenoUnit <- "kg/m²"
  familyName <- "LOGNO"
  family <- LOGNO
  
  sigmaFormula <- paste0(" ~ pregnancy_duration_1")
  
  formula <- paste0(phenoName, " ~ fp(pregnancy_duration_1)")
  formulaLinear <- paste0(phenoName, " ~ pregnancy_duration_1")
  
  refvalues <- values %>% 
    filter(
      unrelated_children == 1
    )
  
  zDF <- standardizeBySex(
    trainingvalues = refvalues,
    values = values,
    id = "child_id",
    x = "pregnancy_duration_1",
    y = phenoName,
    zY = zPhenoName,
    formula = as.formula(formula),
    sigmaFormula = as.formula(sigmaFormula),
    formula2 = as.formula(formulaLinear),
    sigmaFormula2 = as.formula(sigmaFormula),
    family = family
  )
  
  values <- values %>%
    left_join(
      zDF,
      by = "child_id"
    )
  
  writeDocs(
    phenoName = phenoName,
    zPhenoName = zPhenoName,
    phenoLabel = phenoLabel,
    phenoUnit = phenoUnit,
    formula = formula,
    sigmaFormula = sigmaFormula,
    family = familyName,
    values = values,
    docsFolder = docs_folder,
    docsFile = docs_file
  )
}


# Save phenotypes to a new table

print(paste(Sys.time(), "Saving the results."))

values_to_export <- values %>% 
  select(
    all_of(c(default_columns, "pregnancy_duration", age_columns, length_columns, weight_columns, bmi_columns, paste0("z_", length_columns), paste0("z_", weight_columns), paste0("z_", bmi_columns)))
  )


write.table(
  x = values_to_export,
  file = gzfile(file.path(tablesFolder, "child_anthropometrics_standardized.gz")),
  row.names = F,
  col.names = T,
  sep = "\t",
  quote = F
)

