
##
#
# Writes the main documentation.
#
##


# Libraries

library(here)
library(conflicted)
library(yaml)
library(glue)
library(dplyr)

# Solve name space conflicts
conflicts_prefer(dplyr::filter)


# Files

args <- commandArgs(TRUE)

config_file <- args[1]

if (!file.exists(config_file)) {
  
  stop(glue("Config file {config_file} not found."))
  
}

moba_delivery_folder <- args[2]
raw_tables_folder <- args[3]
tables_folder <- args[4]

# Import config file

config <- read_yaml(config_file)


# Write readme for release

readme_file <- file.path(here(), "docs", "anthropometrics", "README.md")

write(
  x = glue("# Anthropometric traits\n"),
  file = readme_file,
  append = F
)
write(
  x = glue("Phenotype analysis pipeline for anthropometric traits in the [Norwegian Mother, Father and Child Cohort Study (MoBa)](fhi.no/en/studies/moba).\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("### Current release\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("- [{config$release}]({config$release}/README.md)\n"),
  file = readme_file,
  append = T
)


write(
  x = glue("#### Previous releases\n"),
  file = readme_file,
  append = T
)

for (file in list.files(file.path(here(), "docs", "anthropometrics"))) {
  
  if (dir.exists(file.path(here(), "docs", "anthropometrics", file)) && file != config$release) {
    
    write(
      x = glue("- [{file}](file)\n"),
      file = readme_file,
      append = T
    )
    
  }
}


# Write readme for release

readme_file <- file.path(here(), "docs", "anthropometrics", config$release, "README.md")

write(
  x = glue("# Anthropometric traits\n"),
  file = readme_file,
  append = F
)
write(
  x = glue("Phenotype analysis pipeline for anthropometric traits in the [Norwegian Mother, Father and Child Cohort Study (MoBa)](fhi.no/en/studies/moba).\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("The current version of the pipeline and associated files is {config$release}. The phenotypes obtained from MoBa correspond to the version {config$moba_version} released for project number {config$moba_project_number} and linkage {config$moba_linkage_relase_date}. The individual-level files are available on HUNT Cloud in the _{config$hunt_lab}_ digital lab, see below for individual paths.\n"),
  file = readme_file,
  append = T
)

write(
  x = glue("### Raw tables\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("Raw tables refer to the original tables obtained from MoBa. These tables contain the raw information from registries and questionnaires. Information on the different tables and variables can be found on the [MoBa website](https://www.fhi.no/en/studies/moba/for-forskere-artikler/questionnaires-from-moba/).\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("- Path to files from MoBa: `{moba_delivery_folder}`.\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("- [Raw tables documentation](raw/data.md).\n"),
  file = readme_file,
  append = T
)

write(
  x = glue("### Curated tables\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("Curated tables refer to tables curated and consolidated by this pipeline.\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("- Path to files from MoBa: `{tables_folder}`.\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("- [Curated tables documentation](tables/phenotypes.md).\n"),
  file = readme_file,
  append = T
)