library(tidyverse)

# Read in data
resource_data <- read_csv("data/processed/final_data/data_resource_20240604.csv")
resource_meta <- read_csv("data/processed/final_data/metadata_resource_20240715.csv")

# Filter to match data & metadata samples
resource_meta <- 
  resource_meta |> 
  filter(DAid %in% resource_data$DAid) |> 
  mutate(Disease = ifelse(Disease == "Scleroderma", "Systemic sclerosis", Disease),
         Age = ifelse(Age == 0 & Class != "Pediatric", NA, Age)) 

resource_data <- 
  resource_data |> 
  filter(DAid %in% resource_meta$DAid)

# Check number of samples
resource_meta |> 
  distinct(DAid)

resource_data |> 
  distinct(DAid)

resource_meta |> 
  count(Disease) |> pull(n) |> sum()
  filter(is.na(Disease))

# Write files
write_tsv(resource_data, "data/final_data/data_hdba_hpa_p1_v24.tsv")
write_tsv(resource_meta, "data/final_data/meta_hdba_hpa_p1_v24.tsv")
