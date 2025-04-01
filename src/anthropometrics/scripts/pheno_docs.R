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
child_id_file <- args[4]
mother_id_file <- args[5]
father_id_file <- args[6]
variables_mapping_file <- args[7]
tables_folder <- args[8]
docs_folder <- args[9]


# The variable mapping
source("src/anthropometrics/scripts/utils/variables_mapping.R")

variable_mapping <- read.table(
  file = variables_mapping_file,
  header = T
)

tables <- unique(variable_mapping$project_table)

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

for (table_name in tables) {
  
  print(paste(Sys.time(), " Processing", table_name))
  
  dir.create(file.path(docs_folder, table_name))
  
  table_file <- file.path(tablesFolder, paste0(project_table_name, ".gz"))
  
  table <- read.table(
    file = table_file,
    header = T,
    sep = "\t"
  )
  
  write(
    x = glue("## {table_name}\n\n"), 
    file = pheno_file, 
    append = T
  )
  
  for (column in names(table)) {
    
    if (!column %in% default_columns) {
      
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
      
      if (column %in% variable_mapping$project_variable) {
        
        i <- which(variable_mapping$project_variable)
        
        write(
          x = glue("Variable mapping to `{variable_mapping$moba_variable[i]}` in `{variable_mapping$moba_table[i]}`."), 
          file = column_file, 
          append = T
        )
        
      } else {
        
        write(
          x = glue("Variable created during phenotype curation."), 
          file = column_file, 
          append = T
        )
        
      }
      
      write(
        x = glue("- Number of values:\n\n"), 
        file = column_file, 
        append = T
      )
      
      write(
        x = glue("| Value | Total | Child genotyped | Mother genotyped | Father genotyped |"), 
        file = column_file, 
        append = T
      )
      
      write(
        x = glue("| ----- | ----- | --------------- | ---------------- | ---------------- |"), 
        file = column_file, 
        append = T
      )
      
      all_values <- table[[column]]
      child_genotyped_values <- table[[column]][!is.na(table$child_sentrix_id)]
      mother_genotyped_values <- table[[column]][!is.na(table$mother_sentrix_id)]
      father_genotyped_values <- table[[column]][!is.na(table$father_sentrix_id)]
      
      write(
        x = glue("| Missing | {sum(is.na(all_values))} | {sum(is.na(child_genotyped_values))} | {sum(is.na(mother_genotyped_values))} | {sum(is.na(father_genotyped_values))} |"), 
        file = column_file, 
        append = T
      )
      write(
        x = glue("| Non-missing | {sum(!is.na(all_values))} | {sum(!is.na(child_genotyped_values))} | {sum(!is.na(mother_genotyped_values))} | {sum(!is.na(father_genotyped_values))} |"), 
        file = column_file, 
        append = T
      )
      
      all_values_numeric <- as.numeric(all_values)
      all_character <- all_values[!is.na(all_values) & is.na(all_values_numeric)]
      child_values_numeric <- all_values_numeric[!is.na(table$child_sentrix_id)]
      mother_values_numeric <- all_values_numeric[!is.na(table$mother_sentrix_id)]
      father_values_numeric <- all_values_numeric[!is.na(table$father_sentrix_id)]
      all_values_numeric <- all_values_numeric[!is.na(all_values_numeric)]
      child_values_numeric <- child_values_numeric[!is.na(child_values_numeric)]
      mother_values_numeric <- mother_values_numeric[!is.na(mother_values_numeric)]
      father_values_numeric <- father_values_numeric[!is.na(father_values_numeric)]
      
      if (length(all_character) > 0) {
        
        unique_character <- sort(unique(all_character))
        
        for (value in unique_character) {
          
          write(
            x = glue("| {value} | {sum(!is.na(all_values) & all_values == value)} | {sum(!is.na(child_genotyped_values) & child_genotyped_values == value)} | {sum(!is.na(mother_genotyped_values) & mother_genotyped_values == value)} |{sum(!is.na(father_genotyped_values) & father_genotyped_values == value)} |"),  
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
              x = glue("| {value} | {sum(all_values_numeric == value)} | {sum(child_values_numeric == value)} | {sum(mother_values_numeric == value)} | {sum(father_values_numeric == value)} |"), 
              file = column_file, 
              append = T
            )
            
          }
        } else {
          
          write(
            x = glue("| 25th percentile | {quantile(all_values_numeric, 0.25)} | {quantile(child_values_numeric, 0.25)} | {quantile(mother_values_numeric, 0.25)} | {quantile(father_values_numeric, 0.25)} |"), 
            file = column_file, 
            append = T
          )
          write(
            x = glue("| 50th percentile | {quantile(all_values_numeric, 0.5)} | {quantile(child_values_numeric, 0.5)} | {quantile(mother_values_numeric, 0.5)} | {quantile(father_values_numeric, 0.5)} |"), 
            file = column_file, 
            append = T
          )
          write(
            x = glue("| 75th percentile | {quantile(all_values_numeric, 0.75)} | {quantile(child_values_numeric, 0.75)} | {quantile(mother_values_numeric, 0.75)} | {quantile(father_values_numeric, 0.75)} |"), 
            file = column_file, 
            append = T
          )
          write(
            x = glue("| Mean | {mean(all_values_numeric, na.rm = T)} | {mean(child_values_numeric, na.rm = T)} | {mean(mother_values_numeric, na.rm = T)} | {mean(father_values_numeric, na.rm = T)} |"), 
            file = column_file, 
            append = T
          )
          write(
            x = glue("| Standard deviation | {sd(all_values_numeric, na.rm = T)} | {sd(child_values_numeric, na.rm = T)} | {sd(mother_values_numeric, na.rm = T)} | {sd(father_values_numeric, na.rm = T)} |"), 
            file = column_file, 
            append = T
          )
          write(
            x = glue("| N | {sum(!is.na(all_values_numeric))} | {sum(!is.na(child_values_numeric))} | {sum(!is.na(mother_values_numeric))} | {sum(!is.na(father_values_numeric))} |"), 
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
          value = c(all_values, child_genotyped_values, mother_genotyped_values, father_genotyped_values),
          category = c(rep("Total", length(all_values)), rep("Child genotyped", length(child_genotyped_values)), rep("Mother genotyped", length(mother_genotyped_values)), rep("Father genotyped", length(father_genotyped_values)))
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
            values = c("black", "darkgreen", "darkblue", "darkred")
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
}

