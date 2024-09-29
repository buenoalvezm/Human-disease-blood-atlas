
# Explore UCAN 

library(OlinkAnalyze)
library(tidyverse)

# Read internal functions
source("scripts/functions/functions_utility.R")
source("scripts/functions/themes_palettes.R")
source("scripts/functions/functions_visualization.R")

# Read data
ucan_ht <- read_NPX("data/raw_data/VL-3530B3_NPX_2024-09-25.parquet")
manifest <- import_df("data/samples_2024-09-27.xlsx")
targets <-  import_df("data/olink_targets/overlap_olink_platforms.csv")


replicate_assays <- 
  ucan_ht |> 
  count(Assay) |> 
  arrange(-n) |> 
  head(2)

replicate_exclude <- 
  ucan_ht |> 
  distinct(Assay, OlinkID, Block) |> 
  filter(Assay %in% replicate_assays$Assay,
         Block != "5") |> 
  pull(OlinkID)

filtered_ucan <- 
  ucan_ht |>
  filter(SampleType == "SAMPLE",
         AssayType == "assay",
         !OlinkID %in% replicate_exclude) |> 
  separate(SampleID, into = c("DAid", "Batch", "Plate"), sep = "-") 

# Number of samples
pal_ucan <- c("LUNG" = "#ADC74F",
              "CRC" = "#B89B74", 
              "BRC" = "#E8A29A",
              "OVC" = "#603479",
              "PRC" = "#E7662B")

filtered_ucan |> 
  distinct(DAid) |> 
  left_join(manifest, by = "DAid") |> 
  filter(Cohort == "UCA2") |> 
  count(Diagnose) |> 
  mutate(Diagnose = factor(Diagnose, levels = names(pal_ucan))) |> 
  ggplot(aes(Diagnose, n, fill = Diagnose)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(values = pal_ucan) +
  theme_hpa(angled = T)

ggsave(savepath("ucan_samples.png"), width = 4, height = 5)

# PCA
pca_ucan <- 
  do_pca(data = filtered_ucan, #|> rename(Sample = DAid), 
       meta = manifest,
       wide = F)

# UMAP
umap_ucan <- 
  do_umap(data = filtered_ucan, 
         meta = manifest,
         wide = F)

# DE
filtered_ucan