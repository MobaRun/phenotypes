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
theme_set(theme_bw(base_size = 14))

pheno_file <- glue("docs/{release_version}/phenotypes.md")
moba_version <- "V12"
release_version <- "23-05-28"
tables_folder <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_23-05-28"

children_id_file <- file.path(tables_folder, "id", "children_id")
mother_id_file <- file.path(tables_folder, "id", "mothers_id")
father_id_file <- file.path(tables_folder, "id", "fathers_id")

pregnancy_table <- file.path(tables_folder, "pregnancy.gz")
delivery_table <- file.path(tables_folder, "delivery.gz")
pregnancy_nutrition_table <- file.path(tables_folder, "pregnancy_nutrition.gz")
mother_nutrition_table <- file.path(tables_folder, "mother_nutrition.gz")
child_nutrition_table <- file.path(tables_folder, "child_nutrition.gz")
child_table <- file.path(tables_folder, "child.gz")
child_health_table <- file.path(tables_folder, "child_health.gz")
parents_table <- file.path(tables_folder, "parents.gz")
mother_health_table <- file.path(tables_folder, "mother_health.gz")
father_health_table <- file.path(tables_folder, "father_health.gz")
child_anthropometrics_table <- file.path(tables_folder, "child_anthropometrics.gz")

project_number <- 315

# The variable mapping
source("src/anthropometrics/scripts/utils/variables_mapping.R")

tables <- list(
  pregnancy_table = pregnancy_table,
  delivery_table = delivery_table,
  pregnancy_nutrition_table = pregnancy_nutrition_table,
  mother_nutrition_table = mother_nutrition_table,
  child_nutrition_table = child_nutrition_table,
  child_table = child_table,
  child_health_table = child_health_table,
  parents_table = parents_table,
  mother_health_table = mother_health_table,
  father_health_table = father_health_table,
  child_anthropometrics_table = child_anthropometrics_table
)
columns <- list(
  pregnancy_table = pregnancy_columns,
  delivery_table = delivery_columns,
  pregnancy_nutrition_table = pregnancy_nutrition_columns,
  mother_nutrition_table = mother_nutrition_columns,
  child_nutrition_table = child_nutrition_columns,
  child_table = child_columns,
  child_health_table = child_health_columns,
  parents_table = parent_values_columns,
  mother_health_table = mother_health_columns,
  father_health_table = father_health_columns,
  child_anthropometrics_table = c(age_columns, weight_columns, length_columns, head_circumference_columns)
)
pheno_to_question <- list(
  mfr = mfrVariablesMapping,
  q1m = q1mVariablesMapping,
  q1f = q1fVariablesMapping,
  q2 = q2VariablesMapping,
  q3 = q3VariablesMapping,
  q4 = q4VariablesMapping,
  q5 = q5VariablesMapping,
  q6 = q6VariablesMapping,
  q7 = q7VariablesMapping,
  q8 = q8VariablesMapping,
  q9 = q9VariablesMapping,
  kost_ungdom = kostUngdomVariablesMapping,
  ungdomsskjema_barn = ungdomsskjemaBarnVariablesMapping
)

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

for (table_name in names(tables)) {
  
  dir.create(glue("docs/{release_version}/{table_name}"))
  
  table <- read.table(
    file = tables[[table_name]],
    header = T,
    sep = "\t"
  )
  
  table_columns <- columns[[table_name]]
  
  write(
    x = glue("## {table_name}\n"), 
    file = pheno_file, 
    append = T
  )
  
  for (column in table_columns) {
    
    column_file <- glue("docs/{release_version}/{table_name}/{column}.md")
    
    write(
      x = glue("[{column}]({table_name}/{column}.md)\n"), 
      file = pheno_file, 
      append = T
    )
    
    write(
      x = glue("# {column}\n"), 
      file = column_file, 
      append = F
    )
    
    questions <- ""
    
    for (questionnaire in names(pheno_to_question)) {
      
      if (column %in% names(pheno_to_question[[questionnaire]])) {
        
        if (nchar(questions) > 0) {
          
          questions <- paste0(questions, "; ")
          
        }
        
        questions <- paste0(questions, "questionnaire: ", questionnaire, ", question ", paste(pheno_to_question[[questionnaire]][[column]], collapse = ","))
        
      }
    }
    
    if (nchar(questions) > 0) {
      
      questions <- paste0(questions, ".\n")
      
      write(
        x = glue("Variable mapping to {questions}."), 
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
    child_genotyped_values <- table[[column]][!is.na(table$sentrix_id)]
    mother_genotyped_values <- table[[column]][!is.na(table$mother_sentrix_id)]
    father_genotyped_values <- table[[column]][!is.na(table$father_sentrix_id)]
    
    write(
      x = glue("| Missing | {sum(is.na(all_values))} | {sum(is.na(child_genotyped_values))} | {sum(is.na(mother_genotyped_values))} | {sum(is.na(father_genotyped_values))} |"), 
      file = column_file, 
      append = T
    )
    
    all_values_numeric <- as.numeric(all_values)
    all_character <- all_values[!is.na(all_values) & is.na(all_values_numeric)]
    all_values_numeric <- all_values_numeric[!is.na(all_values_numeric)]
    
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
      
      if (length(unique_numeric < 20)) {
        
        for (value in unique_numeric) {
          
          write(
            x = glue("| {value} | {sum(!is.na(all_values) & all_values == value)} | {sum(!is.na(child_genotyped_values) & child_genotyped_values == value)} | {sum(!is.na(mother_genotyped_values) & mother_genotyped_values == value)} |{sum(!is.na(father_genotyped_values) & father_genotyped_values == value)} |"), 
            file = column_file, 
            append = T
          )
          
        }
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
          position = "dodge"
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
        filename = glue("docs/{release_version}/{table_name}/{column}_n.png"),
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
