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
release_name <- args[3] 


##
#
# Debug Marc - do not uncomment
# args to run locally on Hunt 
# 
# qcFolder <- "/mnt/work/marc/pheno_22-09-19/qc_tmp"
# project_number <- 2824
# release_name <- "22-09-19
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
library(withr, lib = libFolder)
library(labeling, lib = libFolder)
library(digest, lib = libFolder)
library(ggplot2, lib = libFolder)
library(grid, lib = libFolder)



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
source("src/scripts/utils/anthropometrics/lw-docs-functions.R")

source("src/scripts/utils/anthropometrics/variables_mapping.R")


# Currently exclude values after 8y

length_columns <- length_columns[1:(length(length_columns) - 1)]
weight_columns <- weight_columns[1:(length(weight_columns) - 1)]


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

if (!dir.exists(qcFolder)) {
  
  dir.create(
    path = qcFolder,
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

if (file.exists(qcFile)) {
  
  qcDF <- read.table(qcFile, header = T, stringsAsFactors = F)
  qcDF <- rbind(qcDF, imputativeQcDF)
  
} else {
  
  qcDF <- imputativeQcDF
  
}

write.table(qcDF, file = qcFile, col.names = T, row.names = F, quote = F)


# Write documentation

print(paste(Sys.time(), " Writing documentation"))

docs_folder <- glue("docs/{release_name}/anthropometrics")
md_file <- file.path(docs_folder, "readme.md")
docs_plots_folder <- file.path(docs_folder, "plots")

if (!dir.exists(docs_plots_folder)) {
  
  dir.create(
    path = docs_plots_folder,
    showWarnings = T,
    recursive = T
  )
  
}

write(
  x = "# Phenotypes", 
  file = md_file, 
  append = F
)

write(
  x = "### Number of values", 
  file = md_file, 
  append = T
)


write(
  x = "![](plots/n.png)", 
  file = md_file, 
  append = T
)

write(
  x = "### Length vs weight", 
  file = md_file, 
  append = T
)

for (age_i in 1:length(length_columns)) {
  
  write(
    x = glue("![](plots/length_weight_{age_i}.png)"), 
    file = md_file, 
    append = T
  )
  
}

n_values <- get_n_values(
  originalValues = originalValues,
  values = values
)

write.table(
  x = n_values, 
  file = file.path(docs_folder, "n_values"), 
  col.names = T, 
  row.names = F, 
  quote = F, 
  sep = "\t"
)

n_values <- n_values %>% 
  mutate(
    qc_factor = factor(qc_class, levels = c("Raw", "QC")),
    phenotype_factor = factor(pheno, levels = c("Length", "Weight", "BMI")),
    age_factor = factor(age_i, levels = 1:length(length_columns))
  )

levels(n_values$age_factor) <- c("Birth", "6 w", "3 m", "6 m", "8 m", "1 y", "1.5 y", "2 y", "3 y", "5 y", "7 y", "8 y", "14 y")

n_plot <- ggplot() +
  theme_bw(
    base_size = 24
  ) +
  geom_col(
    data = n_values,
    mapping = aes(
      x = age_factor,
      y = n,
      fill = qc_factor
    ),
    position = "dodge"
  ) +
  facet_grid(
    phenotype_factor ~ .
  ) +
  scale_y_continuous(
    name = "# values",
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  scale_fill_manual(
    values = c("darkblue", "darkgreen")
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

png(
  filename = file.path(docs_plots_folder, "n.png"),
  height = 300 * 3,
  width = 900
)
grid.draw(n_plot)
dummy <- dev.off()


exportLenghtWeight(
    originalValues = originalValues,
    values = values,
    plot_folder = docs_plots_folder
)

write(
  x = "### Imputation", 
  file = md_file, 
  append = T
)

write(
  x = glue("- Children with no data point altered: {sum(values$log == '')}"), 
  file = md_file, 
  append = T
)

write(
  x = glue("- Children with at least one data point altered: {sum(values$log != '')}"), 
  file = md_file, 
  append = T
)

log_df <- values %>% 
  select(
    child_id, log
  ) %>% 
  filter(
    log != ""
  ) %>% 
  mutate(
    log_length = nchar(log)
  ) %>% 
  arrange(
    desc(log_length)
  )

sampled_rows <- which(values$child_id %in% sample(log_df$child_id, 20, replace = F))

for (row in sampled_rows) {
  
  dummyIdI <- bridgeDF$dummyId[row]
  
  write(
    x = glue("#### Random example: {dummyIdI}"), 
    file = md_file, 
    append = T
  )
  
  write(
    x = glue("> {values$log[i]}"), 
    file = md_file, 
    append = T
  )
  
  write(
    x = glue("![](plots/{dummyIdI}_length.png)"), 
    file = md_file, 
    append = T
  )
  
  write(
    x = glue("![](plots/{dummyIdI}_weight.png)"), 
    file = md_file, 
    append = T
  )
  
  plots <- get_annotated_curves(
    originalValues = originalValues,
    values = values,
    i = row
  )
  
  plotFile <- file.path(docs_plots_folder, paste0(dummyIdI, "_length.png"))
  png(plotFile, width = 900, height = 600)
  grid.draw(plots[[1]])
  dummy <- dev.off()
  
  plotFile <- file.path(docs_plots_folder, paste0(dummyIdI, "_weight.png"))
  png(plotFile, width = 900, height = 600)
  grid.draw(plots[[2]])
  dummy <- dev.off()
  
}


for (i in 1:20) {
  
  child_id <- log_df$child_id[i]
  
  row <- which(values$child_id == child_id)
  
  dummyIdI <- bridgeDF$dummyId[row]
  
  write(
    x = glue("#### Most extreme example ({i}): {dummyIdI}"), 
    file = md_file, 
    append = T
  )
  
  write(
    x = glue("> {values$log[i]}"), 
    file = md_file, 
    append = T
  )
  
  write(
    x = glue("![](plots/{dummyIdI}_length.png)"), 
    file = md_file, 
    append = T
  )
  
  write(
    x = glue("![](plots/{dummyIdI}_weight.png)"), 
    file = md_file, 
    append = T
  )
  
  plots <- get_annotated_curves(
    originalValues = originalValues,
    values = values,
    i = row
  )
  
  plotFile <- file.path(docs_plots_folder, paste0(dummyIdI, "_length.png"))
  png(plotFile, width = 900, height = 600)
  grid.draw(plots[[1]])
  dummy <- dev.off()
  
  plotFile <- file.path(docs_plots_folder, paste0(dummyIdI, "_weight.png"))
  png(plotFile, width = 900, height = 600)
  grid.draw(plots[[2]])
  dummy <- dev.off()
  
}

