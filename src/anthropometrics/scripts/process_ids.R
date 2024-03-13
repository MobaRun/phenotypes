##
#
# This script maps genotypic to phenotypic ids and extracts tables of id mapping.
# 
##

set.seed(123456)


# Libraries

library(dplyr)


# Command line parameters

args <- commandArgs(TRUE)

kinship_file <- args[1]
linkage_child <- args[2]
linkage_mother <- args[3]
linkage_father <- args[4]
id_folder <- args[5]
project_number <- args[6]


# Parameters

preg_id_column <- paste0("PREG_ID_", project_number)
mother_id_column <- paste0("M_ID_", project_number)
father_id_column <- paste0("F_ID_", project_number)

# Function

##
# This function takes a data frame of identifiers and exports lists of ids and related ids.
##
process_ids <- function(
    linkage_table,
    relatedness_df,
    fam_id_df,
    export_folder,
    file_name
) {
  
  print(paste0(Sys.time(), " - Processing ", file_name))
  
  identifiers <- linkage_table %>% 
    select(
      iid = SENTRIX_ID, id
    ) %>% 
    inner_join(
      fam_id_df,
      by = "iid"
    ) %>% 
    rename(
      sentrix_id = iid,
      pheno_id = id,
      family_id = fid
    )
  
  # Extract unrelated individuals
  
  unrelated_individuals <- identifiers$sentrix_id[! identifiers$sentrix_id %in% c(related_ids_table$ID1, related_ids_table$ID2)]
  
  relatedness_graph <- graph_from_data_frame(
    related_ids_table,
    directed = F
  )
  
  relatedness_components <- components(relatedness_graph)
  
  for (component_i in 1:relatedness_components$no) {
    
    component_nodes <- names(relatedness_components$membership)[relatedness_components$membership == component_i]
    
    unrelated_individuals <- c(unrelated_individuals, sample(component_nodes, 1))
    
  }

  identifiers_unrelated <- identifiers[identifiers$sentrix_id %in% unrelated_individuals, ]

  print(paste0(Sys.time(), " - Exporting identifiers for ", file_name))
  
  output_file <- file.path(export_folder, file_name)
  
  write.table(
    x = identifiers,
    file = output_file,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
  )
  
  output_file <- file.path(export_folder, paste0(file_name, "_plink"))
  
  write.table(
    x = identifiers[ , c("family_id", "sentrix_id")],
    file = output_file,
    col.names = F,
    row.names = F,
    quote = F,
    sep = " "
  )
  
  
  output_file <- file.path(export_folder, paste0(file_name, "_unrelated"))
  
  write.table(
    x = identifiers_unrelated,
    file = output_file,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
  )
  
  output_file <- file.path(export_folder, paste0(file_name, "_unrelated_plink"))
  
  write.table(
    x = identifiers_unrelated[ , c("family_id", "sentrix_id")],
    file = output_file,
    col.names = F,
    row.names = F,
    quote = F,
    sep = " "
  )


  sampled_indexes <- sort(sample(x = indexes, size = min(length(indexes), 30000)))
  identifiers_unrelated <- identifiers[sampled_indexes, ]
  
  output_file <- file.path(export_folder, paste0(file_name, "_unrelated_30k"))
  
  write.table(
    x = identifiers_unrelated,
    file = output_file,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
  )
  
  output_file <- file.path(export_folder, paste0(file_name, "_unrelated_30k_plink"))

  write.table(
    x = identifiers_unrelated[ , c("family_id", "sentrix_id")],
    file = output_file,
    col.names = F,
    row.names = F,
    quote = F,
    sep = " "
  )
  
}

print(paste0(Sys.time(), " - Loading kinship"))

kinship_df <- read.table(
  file = kinship_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

# Make list of related individuals

relatedness_threshold <- min(kinship_df$PropIBD[kinship_df$InfType == "3rd"])

relatedness_df <- kinship_df[kinship_df$PropIBD > relatedness_threshold, c("ID1", "ID2")]


# Get family ids

fam_id_df_1 <- kinship_df[, c("FID", "ID1")]
names(fam_id_df_1) <- c("fid", "iid")
fam_id_df_2 <- kinship_df[, c("FID", "ID2")]
names(fam_id_df_2) <- c("fid", "iid")
fam_id_df <- rbind(fam_id_df_1, fam_id_df_2) %>% 
  distinct()


# Sanity check that we have only one family id per individual

fid_n <- fam_id_df %>% 
  group_by(
    iid
  ) %>% 
  summarize(
    n = n()
  )

if (max(fid_n$n) > 1) {
  
  stop("More than one family id per individual")
  
}


# Children

print(paste0(Sys.time(), " - Loading identifiers"))

child_linkage_table <- read.table(
  file = linkage_child,
  sep = "\t", 
  header = T
) %>% 
  mutate(
    id = paste0(!!sym(preg_id_column), "_", BARN_NR)
  ) %>% 
  select(
    SENTRIX_ID, id
  ) %>% 
  distinct()
  
process_ids(
    linkage_table = child_linkage_table,
    relatedness_df = relatedness_df,
    fam_id_df = fam_id_df,
    export_folder = id_folder,
    file_name = "children_id"
)

mother_linkage_table <- read.table(
  file = linkage_mother,
  sep = "\t", 
  header = T
) %>% 
  select(
    SENTRIX_ID, id = !!sym(mother_id_column)
  ) %>% 
  distinct()

process_ids(
  linkage_table = mother_linkage_table,
  relatedness_df = relatedness_df,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "mothers_id"
)

father_linkage_table <- read.table(
  file = linkage_father,
  sep = "\t", 
  header = T
) %>% 
  select(
    SENTRIX_ID, id = !!sym(father_id_column)
  ) %>% 
  distinct()

process_ids(
  linkage_table = father_linkage_table,
  relatedness_df = relatedness_df,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "fathers_id"
)

parents_linkage_table <- rbind(mother_linkage_table, father_linkage_table)
process_ids(
  linkage_table = parents_linkage_table,
  relatedness_df = relatedness_df,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "parents_id"
)

all_individuals_linkage_table <- rbind(parents_linkage_table, child_linkage_table)
process_ids(
  linkage_table = all_individuals_linkage_table,
  relatedness_df = relatedness_df,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "all_individuals_id"
)

