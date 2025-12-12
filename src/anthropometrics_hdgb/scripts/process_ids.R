##
#
# This script maps genotypic to phenotypic ids and extracts tables of id mapping.
# 
##

set.seed(123456)


# Libraries

library(janitor)
library(dplyr)
library(igraph)


# Command line parameters

args <- commandArgs(TRUE)

kinship_file <- args[1]
psam_file <- args[2]
linkage_preg <- args[3]
linkage_child <- args[4]
linkage_mother <- args[5]
linkage_father <- args[6]
id_folder <- args[7]
project_number <- args[8]


# Debug files
# 
# kinship_file <- "/mnt/work/mikko/github/new_geno_release_test/kinship_test/moba_genotypes_2025.01.30_withFIDn_clean.kin0"
# linkage_preg <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-05-07/raw/linkage/PDB315_SV_INFO_V12_20250131.gz"
# linkage_child <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-05-07/raw/linkage/PDB315_MoBaGeneticsTot_Child_20221228.gz"
# linkage_mother <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-05-07/raw/linkage/PDB315_MoBaGeneticsTot_Mother_20221228.gz"
# linkage_father <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-05-07/raw/linkage/PDB315_MoBaGeneticsTot_Father_20221228.gz"
# id_folder <- "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-05-07/id"
# project_number <- 315
#


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
    related_ids_table,
    fam_id_df,
    export_folder,
    file_name,
    population
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
  
  # Keep only relationships within population
  
  population_related_ids_table <- related_ids_table %>% 
    filter(
      ID1 %in% identifiers$sentrix_id & ID2 %in% identifiers$sentrix_id
    )
  
  
  # Extract unrelated individuals
  
  unrelated_individuals <- identifiers$sentrix_id[! identifiers$sentrix_id %in% c(population_related_ids_table$ID1, population_related_ids_table$ID2)]
  
  if (nrow(population_related_ids_table) > 0) {
    
    relatedness_graph <- graph_from_data_frame(
      population_related_ids_table,
      directed = F
    )
    
    relatedness_components <- components(relatedness_graph)
    
    for (component_i in 1:relatedness_components$no) {
      
      component_nodes <- names(relatedness_components$membership)[relatedness_components$membership == component_i]
      
      if (length(component_nodes) > 0) {
        
        unrelated_individuals <- c(unrelated_individuals, sample(component_nodes, 1))
        
      }
    }
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
  
  n_unrelated <- nrow(identifiers_unrelated)
  
  if (n_unrelated <= 30000) {
    
    stop(paste0("Less than 30000 unrelated individuals found (", n_unrelated, ") for ", population, "."))
    
  }
  
  sampled_indexes <- sort(sample(x = 1:nrow(identifiers_unrelated), size = 30000))
  identifiers_unrelated <- identifiers_unrelated[sampled_indexes, ]
  
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

print(paste0(Sys.time(), " - Loading psam file"))

psam_df <- read.table(
  file = psam_file,
  header = F,
  sep = "\t",
  stringsAsFactors = F
)
names(psam_df) <- c("fid", "iid", "pat", "mat", "sex")


# Make list of related individuals

relatedness_threshold <- min(kinship_df$PropIBD[kinship_df$InfType == "3rd"])

related_ids_table <- kinship_df[kinship_df$PropIBD > relatedness_threshold, c("ID1", "ID2")]

# Get family ids

fam_id_df <- psam_df[, c("fid", "iid")]


# Children

print(paste0(Sys.time(), " - Loading identifiers"))

trioIdDF <- read.table(
  file = linkage_preg,
  sep = "\t",
  header = T,
  quote = "",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  select(
    preg_id = !!sym(preg_id_column),
    mother_id = !!sym(mother_id_column),
    father_id = !!sym(father_id_column)
  )

child_linkage_table <- read.table(
  file = linkage_child,
  sep = "\t", 
  header = T
) %>% 
  rename(
    preg_id = !!sym(preg_id_column)
  ) %>% 
  filter(
    preg_id %in% trioIdDF$preg_id
  ) %>% 
  mutate(
    id = paste0(preg_id, "_", BARN_NR)
  ) %>% 
  select(
    SENTRIX_ID, id
  ) %>% 
  distinct()

process_ids(
  linkage_table = child_linkage_table,
  related_ids_table = related_ids_table,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "children_id",
  "children"
)

mother_linkage_table <- read.table(
  file = linkage_mother,
  sep = "\t", 
  header = T
) %>% 
  select(
    SENTRIX_ID, id = !!sym(mother_id_column)
  ) %>% 
  filter(
    id %in% trioIdDF$mother_id
  ) %>% 
  distinct()

process_ids(
  linkage_table = mother_linkage_table,
  related_ids_table = related_ids_table,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "mothers_id",
  "mothers"
)

father_linkage_table <- read.table(
  file = linkage_father,
  sep = "\t", 
  header = T
) %>% 
  select(
    SENTRIX_ID, id = !!sym(father_id_column)
  ) %>% 
  filter(
    id %in% trioIdDF$father_id
  ) %>% 
  distinct()

process_ids(
  linkage_table = father_linkage_table,
  related_ids_table = related_ids_table,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "fathers_id",
  "fathers"
)

parents_linkage_table <- rbind(mother_linkage_table, father_linkage_table)
process_ids(
  linkage_table = parents_linkage_table,
  related_ids_table = related_ids_table,
  fam_id_df = fam_id_df,
  export_folder = id_folder,
  file_name = "parents_id",
  "parents"
)

