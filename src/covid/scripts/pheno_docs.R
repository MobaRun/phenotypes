##
#
# This script builds basic documentation on the phenotypes.
# 
##

library(conflicted)
library(janitor)
library(glue)
library(dplyr)
library(ggplot2)
library(grid)

conflicts_prefer(dplyr::filter)

# General parameters
theme_set(theme_bw(base_size = 16))

# Command line arguments
args <- commandArgs(TRUE)

moba_version <- args[1]
release_version <- args[2]
project_number <- args[3]
covid_table <- args[4]
docs_folder <- args[5]

# moba_version <- "V12"
# release_version <- "23-07-12"
# project_number <- 315
# 
# covid_table <- "/mnt/work/marc/phenotypes/pheno_covid_23-05-28/covid/moba_covid_phenotypes.gz"

pheno_file <- file.path(docs_folder, "phenotypes.md")

write(
  x = glue("# Phenotypes\n"), 
  file = pheno_file, 
  append = F
)
write(
  x = glue("Descriptives of the phenotypic files for version '{moba_version}' of the phenotypes, release '{release_version}' of the QC pipeline.\n"), 
  file = pheno_file, 
  append = T
)

table_name <- "moba_covid_phenotypes"

print(paste(Sys.time(), " Processing", table_name))

dir.create(file.path(docs_folder, table_name))

table <- read.table(
  file = covid_table,
  header = T,
  sep = "\t"
)

write(
  x = glue("## {table_name}\n\n"), 
  file = pheno_file, 
  append = T
)

excluded_columns <- c("id", "sentrix_id", "role", "batch")

for (column in names(table)) {
  
  if (!column %in% excluded_columns) {
    
    column_file <- file.path(docs_folder, table_name, glue("{column}.md"))
    
    write(
      x = glue("- [{column}]({table_name}/{column}.md)\n"), 
      file = pheno_file, 
      append = T
    )
    
    write(
      x = glue("# {column}\n"), 
      file = column_file, 
      append = F
    )
    
    write(
      x = glue("- Number of values:\n\n"), 
      file = column_file, 
      append = T
    )
    
    all_values <- table[[column]]
    child_genotyped_values <- table[[column]][!is.na(table$sentrix_id) & table$role == "Child"]
    mother_genotyped_values <- table[[column]][!is.na(table$sentrix_id) & table$role == "Mother"]
    father_genotyped_values <- table[[column]][!is.na(table$sentrix_id) & table$role == "Father"]
    parent_genotyped_values <- table[[column]][!is.na(table$sentrix_id) & (table$role == "Mother" | table$role == "Father")]
    
    write(
      x = glue("| Value | Total | Child genotyped | Mother genotyped | Father genotyped | Parents genotyped |"), 
      file = column_file, 
      append = T
    )
    
    write(
      x = glue("| ----- | ----- | --------------- | ---------------- | ---------------- |---------------- |"), 
      file = column_file, 
      append = T
    )
    
    write(
      x = glue("| Missing | {sum(is.na(all_values))} | {sum(is.na(child_genotyped_values))} | {sum(is.na(mother_genotyped_values))} | {sum(is.na(father_genotyped_values))} | {sum(is.na(parent_genotyped_values))} |"), 
      file = column_file, 
      append = T
    )
    write(
      x = glue("| Non-missing | {sum(!is.na(all_values))} | {sum(!is.na(child_genotyped_values))} | {sum(!is.na(mother_genotyped_values))} | {sum(!is.na(father_genotyped_values))} | {sum(!is.na(parent_genotyped_values))} |\n\n"), 
      file = column_file, 
      append = T
    )
    
    
    write(
      x = glue("| Value | Total | Child genotyped | Mother genotyped | Father genotyped | Parents genotyped |"), 
      file = column_file, 
      append = T
    )
    
    write(
      x = glue("| ----- | ----- | --------------- | ---------------- | ---------------- |---------------- |"), 
      file = column_file, 
      append = T
    )
    
    all_values_numeric <- as.numeric(all_values)
    all_character <- all_values[!is.na(all_values) & is.na(all_values_numeric)]
    child_values_numeric <- all_values_numeric[!is.na(table$sentrix_id) & table$role == "Child"]
    mother_values_numeric <- all_values_numeric[!is.na(table$sentrix_id) & table$role == "Mother"]
    father_values_numeric <- all_values_numeric[!is.na(table$sentrix_id) & table$role == "Father"]
    parent_values_numeric <- all_values_numeric[!is.na(table$sentrix_id) & (table$role == "Mother" | table$role == "Father")]
    all_values_numeric <- all_values_numeric[!is.na(all_values_numeric)]
    child_values_numeric <- child_values_numeric[!is.na(child_values_numeric)]
    mother_values_numeric <- mother_values_numeric[!is.na(mother_values_numeric)]
    father_values_numeric <- father_values_numeric[!is.na(father_values_numeric)]
    parent_values_numeric <- parent_values_numeric[!is.na(parent_values_numeric)]
    
    if (length(all_character) > 0) {
      
      unique_character <- sort(unique(all_character))
      
      for (value in unique_character) {
        
        write(
          x = glue("| {value} | {sum(!is.na(all_values) & all_values == value)} | {sum(!is.na(child_genotyped_values) & child_genotyped_values == value)} | {sum(!is.na(mother_genotyped_values) & mother_genotyped_values == value)} | {sum(!is.na(father_genotyped_values) & father_genotyped_values == value)} | {sum(!is.na(father_genotyped_values) & father_genotyped_values == value)} | {sum(!is.na(parent_genotyped_values) & parent_genotyped_values == value)} |"),  
          file = column_file, 
          append = T
        )
        
      }
    }
    
    if (length(all_values_numeric) > 0) {
      
      unique_numeric <- sort(unique(all_values_numeric))
      
      if (length(unique_numeric) < 20) {
        
        for (value in unique_numeric) {
          
          write(
            x = glue("| {value} | {sum(all_values_numeric == value)} | {sum(child_values_numeric == value)} | {sum(mother_values_numeric == value)} | {sum(father_values_numeric == value)} | {sum(parent_values_numeric == value)} |"), 
            file = column_file, 
            append = T
          )
          
        }
      } else {
        
        write(
          x = glue("| 25th percentile | {quantile(all_values_numeric, 0.25)} | {quantile(child_values_numeric, 0.25)} | {quantile(mother_values_numeric, 0.25)} | {quantile(father_values_numeric, 0.25)} | {quantile(parent_values_numeric, 0.25)} |"), 
          file = column_file, 
          append = T
        )
        write(
          x = glue("| 50th percentile | {quantile(all_values_numeric, 0.5)} | {quantile(child_values_numeric, 0.5)} | {quantile(mother_values_numeric, 0.5)} | {quantile(father_values_numeric, 0.5)} | {quantile(parent_values_numeric, 0.5)} |"), 
          file = column_file, 
          append = T
        )
        write(
          x = glue("| 75th percentile | {quantile(all_values_numeric, 0.75)} | {quantile(child_values_numeric, 0.75)} | {quantile(mother_values_numeric, 0.75)} | {quantile(father_values_numeric, 0.75)} | {quantile(parent_values_numeric, 0.75)} |"), 
          file = column_file, 
          append = T
        )
        
      }
      
      write(
        x = "\n\n",  
        file = column_file, 
        append = T
      )
      
      plot_data <- data.frame(
        value = c(all_values, child_genotyped_values, mother_genotyped_values, father_genotyped_values, parent_genotyped_values),
        category = c(rep("Total", length(all_values)), rep("Child genotyped", length(child_genotyped_values)), rep("Mother genotyped", length(mother_genotyped_values)), rep("Father genotyped", length(father_genotyped_values)), rep("Parent genotyped", length(parent_genotyped_values)))
      ) %>% 
        mutate(
          value = as.numeric(value)
        ) %>% 
        filter(
          !is.na(value)
        )
      
      plot <-  ggplot() +
        geom_histogram(
          data = plot_data,
          mapping = aes(
            x = value,
            fill = category
          ),
          position = "dodge",
          bins = 30
        ) +
        scale_fill_manual(
          values = c("black", "darkgreen", "darkblue", "darkred", "purple")
        ) +
        scale_y_continuous(
          name = "Number of values",
          expand = expansion(mult = c(0, .1))
        ) +
        theme(
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
      
      png(
        filename = file.path(docs_folder, table_name, glue("{column}_n.png")),
        width = 900,
        height = 600
      )
      grid.draw(plot)
      device <- dev.off()
      
      write(
        x = glue("![]({column}_n.png)"), 
        file = column_file, 
        append = T
      )
      
    }
    
    write(
      x = "\n\n",  
      file = column_file, 
      append = T
    )
  }
}
