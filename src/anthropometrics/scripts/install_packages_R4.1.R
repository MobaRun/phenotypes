
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

# Install all packages required by the pipeline

packages <- c(
  "janitor", 
  "conflicted", 
  "foreign", 
  "stringr", 
  "glue", 
  "crayon", 
  "tidyr", 
  "dplyr", 
  "purrr",
  "withr",
  "labeling",
  "digest",
  "farver",
  "ggplot2",
  "grid"
)

for (package in packages) {
  
  install.packages(package, lib = library, repos = repo)
  
  library(package, lib.loc = library)
  
}

# Log 

args <- commandArgs(TRUE)

log_file <- args[1]

write(x = glue("{Sys.time()} - Done"), file = log_file, append = F)

