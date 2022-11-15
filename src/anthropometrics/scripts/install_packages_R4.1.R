
##
#
# This script installs the packages needed for the conda R_4.1 environment.
#
# Run in the same conda environment as the scripts that will run using R 4.1 environment.
# Execute from the root of the repo.
#
##

# Parameters

library <- "~/R/R_4.1"
repo <- "https://cran.uib.no/"


# Housekeeping

if (!file.exists(library)) {
  
  dir.create(
    path = library,
    showWarnings = T,
    recursive = T
  )
  
}


# Install all packages required by the pipeline

packages <- c(
  "memoise", 
  "crayon", 
  "conflicted", 
  "stringr", 
  "glue", 
  "purrr",
  "dplyr", 
  "janitor", 
  "foreign", 
  "tidyr", 
  "withr",
  "labeling",
  "digest",
  "farver",
  "ggplot2"
)

for (package in packages) {
  
  install.packages(package, lib = library, repos = repo)
  
  library(package, lib.loc = library, character.only = T)
  
}

# Log 

args <- commandArgs(TRUE)

log_file <- args[1]

write(x = glue("{Sys.time()} - Done"), file = log_file, append = F)

