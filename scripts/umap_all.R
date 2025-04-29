
# UMAP all batches
library(tidyverse)
source("scripts/functions/functions_analyses.R")
source("scripts/functions/functions_visualization.R")
source("scripts/functions/functions_utility.R")
source("scripts/functions/themes_palettes.R")
data_b1 <- read_tsv("data/final_data/data_phase2_batch1_curated_20241217.tsv") 
data_b3 <- read_tsv("data/final_data/data_phase2_batch2_curated_20250127.tsv") 
data_b4 <-  read_tsv("data/final_data/data_phase2_batch4_curated_20250425.tsv") 
manifest <- import_df("data/samples_2025-04-07.xlsx")

# Combine
combined_data <- 
  data_b1 |> 
  select(DAid, Assay, NPX) |> 
  mutate(Batch = "1") |> 
  bind_rows(data_b3 |> 
              select(DAid, Assay, NPX) |> 
              mutate(Batch = "3")) |> 
  bind_rows(data_b4 |> 
              select(DAid, Assay, NPX) |> 
              mutate(Batch = "4")) |> 
  mutate(DAid = paste(DAid, Batch, sep = "_")) |> 
  mutate(DAid_original = str_remove(DAid, "_[0-9]+$")) 
  


umap <-
  do_umap(data = combined_data,
          wide = F,
          plots = F)

saveRDS(umap, savepath_data(folder = "UMAP", savename = "umap.rds"))
umap <- readRDS(savepath_data(folder = "UMAP", savename = "umap.rds"))


meta <- 
  combined_data |> 
  distinct(DAid, Batch)

saveRDS(meta, savepath_data(folder = "UMAP", savename = "umap_meta.rds"))

umap |>
  rename(DAid = Sample) |>
  left_join(meta, by = "DAid") |>
  ggplot(aes(UMAP1, UMAP2, color = Batch)) +
  geom_point(alpha = 0.7, size = 0.8) +
  theme_hpa() +
  ggtitle("UMAP: B1, B3, B4") +
  
  umap |>
  rename(DAid = Sample) |>
  mutate(DAid = str_remove(DAid, "_[0-9]+$")) |> 
  left_join(meta, by = "DAid") |>
  left_join(manifest, by = "DAid") |> #distinct(Class)
  ggplot(aes(UMAP1, UMAP2, color = Cohort)) +
  geom_point(alpha = 0.7, size = 0.8) +
#  scale_color_manual(values = pal_class) +
  theme_hpa() 

ggsave(savepath("umap_b1_b3_b4.png"), h = 5, w = 10)

umap |>
  rename(DAid = Sample) |>
  left_join(manifest, by = "DAid") |>
  mutate(Cohort = paste(Cohort, `PI - sample owner`, sep = "_")) |>
  ggplot(aes(UMAP1, UMAP2, color = Cohort)) +
  geom_point(alpha = 0.7, size = 0.8) +
 # scale_color_manual(values = pal_phase2) +
 # theme_hpa() +
  ggtitle("UMAP: B1, B3, B4")

umap |> 
  filter(UMAP1 > 0, 
         UMAP2 > 5) |> 
  rename(DAid = Sample) |> 
  mutate(DAid = str_remove(DAid, "_[0-9]+$")) |> 
  left_join(manifest, by = "DAid") |> 
#  count(Cohort)
  filter(Cohort == "UCA2")


# Pregancy
meta_preg  <- 
  manifest |> 
  filter(DAid %in% unique(combined_data$DAid_original),
         Sex == "F") |> 
  select(DAid, Disease, Cohort) |> 
  mutate(Disease = ifelse(Cohort == "PREG", "Pregnancy", "Control")) |> 
  select(-Cohort)

data_preg <- 
  combined_data |> 
  select(DAid = DAid_original, Assay, NPX) |> 
  filter(DAid %in% meta_preg$DAid)

de_preg <- 
  do_limma_disease(data_wide = data_preg, 
                 metadata = meta_preg,
                 disease = "Pregnancy",
                 correct = F,
                 cutoff = 0.5)

# Only women
data_preg |> 
  left_join


         