
# Prepare olink to uniprot mapping

library(tidyverse)
resource_data <- read_csv("data/processed/final_data/data_resource_20240604.csv")

resource_data |> 
  distinct(Assay, UniProt) |> 
  write_tsv("data/olink_uniprot_mapping.tsv")
