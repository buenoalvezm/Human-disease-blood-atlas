# Prepare PDAC data
library(tidyverse)

source("scripts/functions/functions_utility.R")

data_b4 <-  read_tsv("data/final_data/data_phase2_batch4_curated_20250425.tsv") 
manifest <- import_df("data/samples_2025-05-12.xlsx")

data_b4 |> 
  left_join(manifest, by = "DAid") |> 
  count(Disease) |> 
  pull(Disease)

pdac_samples <- 
  manifest |> 
  filter(Disease == "Pancreatic ductal adenocarcinoma")

pdac_data <- 
  data_b4 |> 
  filter(DAid %in% pdac_samples$DAid) 


write_tsv(pdac_data, "data/processed/pdac_ht_data.tsv")
