
# Prepare olink to uniprot mapping

library(tidyverse)
resource_data <- read_csv("data/processed/final_data/data_resource_20240604.csv")

resource_data |> 
  distinct(Assay, UniProt) |> 
  write_tsv("data/olink_uniprot_mapping.tsv")

# For HR data
data_bamse <- read_tsv("data/final_data/HPA/v24_2/bamse_data_ht_phase2.tsv")

data_bamse |> 
  distinct(Assay, OlinkID, UniProt) |> 
  write_tsv("data/olink_uniprot_mapping_ht.tsv")