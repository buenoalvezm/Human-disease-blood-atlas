
# Load libraries
library(tidyverse)
library(readxl)
library(OlinkAnalyze)

# Read Olink HT data file and manifest
data_b4 <-  read_NPX("data/VL-3530B4_NPX_2025-03-25.parquet")
manifest <- read_excel("../metadata/samples_2025-05-12.xlsx")

# Extract PCAD data & metadata
pdac_meta <- 
  manifest |> 
  filter(Cohort == "PAC1") 

pdac_mapping <- 
  data_b4 |> 
  distinct(SampleID) |> 
  separate(SampleID, into = c("DAid", "Batch", "Plate"), remove = FALSE) |> 
  filter(DAid %in% pdac_meta$DAid)

pdac_data <- 
  data_b4 |> 
  filter(SampleID %in% pdac_mapping$SampleID)

# Select 16 bridge samples in the PDAC data
set.seed(213)
pdac_bridge_samples <- 
  pdac_data |> 
  olink_bridgeselector(sampleMissingFreq = 0.1,
                     n = 16) 

# Plot PCA highlighting bridging samples
pdac_data |> 
  mutate(Bridge = ifelse(SampleID %in% pdac_bridge_samples$SampleID, "Bridge", "Sample")) |> 
  olink_pca_plot(color_g = "Bridge")

# Save briding samples
pdac_bridge_selection <- 
  pdac_mapping |> 
  filter(SampleID %in% pdac_bridge_samples$SampleID) |> 
  left_join(manifest, by = "DAid") |> 
  select(SampleID, DAid, `Vial barcode`) |> 
  relocate(`Vial barcode`, "DAid")

write_tsv(pdac_bridge_selection, "data/processed/pdac_bridge_samples_20251016.tsv")



