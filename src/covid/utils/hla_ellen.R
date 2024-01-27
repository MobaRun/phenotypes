location_table <- read.table("~/mnt/work/marc/phenotypes/pheno_covid_23-09-08/hla_location.gz", header = T, sep = "\t")

linkage_child <- read.table("~/mnt/work/marc/phenotypes/pheno_covid_23-09-08/raw/linkage/20220516_MoBaGeneticsTot_Child_PDB2824.gz", header = T, sep = "\t")
linkage_mother <- read.table("~/mnt/work/marc/phenotypes/pheno_covid_23-09-08/raw/linkage/20220516_MoBaGeneticsTot_Mother_PDB2824.gz", header = T, sep = "\t")
linkage_father <- read.table("~/mnt/work/marc/phenotypes/pheno_covid_23-09-08/raw/linkage/20220516_MoBaGeneticsTot_Father_PDB2824.gz", header = T, sep = "\t")

library(dplyr)

location_matched <- location_table %>% 
  left_join(
    linkage_mother %>% 
      select(
        id = M_ID_2824,
        SENTRIX_ID_m = SENTRIX_ID
      ),
    by = "id",
    relationship = "many-to-many"
  ) %>% 
  left_join(
    linkage_father %>% 
      select(
        id = F_ID_2824,
        SENTRIX_ID_f = SENTRIX_ID
      ),
    by = "id",
    relationship = "many-to-many"
  ) %>% 
  mutate(
    SENTRIX_ID = ifelse(is.na(SENTRIX_ID_m), SENTRIX_ID_f, SENTRIX_ID_m)
  ) %>% 
  filter(
    !is.na(SENTRIX_ID)
  ) %>% 
  select(
    SENTRIX_ID, fylkenummer, fylkenavn
  ) %>% 
  distinct() %>% 
  group_by(
    SENTRIX_ID
    ) %>% 
  summarize(
    # fylkenummer = paste(fylkenummer, sep = ","),
    fylkenavn = paste(fylkenavn, sep = ",")
  )

genes <- c("A", "B", "C", "DPB1", "DQB1", "DRB1")

for (gene in genes) {
  
  hibag_file <- paste0("hibag_", gene, ".csv")
  hibag_table <- read.table(hibag_file, sep = ",", header = T)
  
  hibag_location <- location_matched %>% 
    inner_join(
      hibag_table %>% 
        rename(
          SENTRIX_ID = sample.id
        ),
      by = "SENTRIX_ID"
    ) %>% 
    select(
      -SENTRIX_ID
    )
  
  location_file <- paste0("~/hla/", gene, "_location.gz")
  
  write.table(
    hibag_location,
    gzfile(location_file),
    col.names = T,
    sep = "\t"
  )
  
}

