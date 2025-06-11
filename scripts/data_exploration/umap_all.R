
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

data_b5 <-  read_tsv(file = "data/final_data/data_phase2_batch5_raw_20250512.tsv")
data_b2 <-  read_tsv(file = "data/final_data/data_phase2_batch2_raw_20250512.tsv")


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
  select(DAid, Disease, Cohort, Sex, Age, BMI) |> 
  mutate(Disease = ifelse(Cohort == "PREG", "Pregnancy", "Control")) |> 
  select(-Cohort)

data_preg <- 
  combined_data |> 
  select(DAid = DAid_original, Assay, NPX) |> 
  filter(DAid %in% meta_preg$DAid) |> 
  pivot_wider(names_from = Assay, values_from = NPX) 

de_preg <- 
  do_limma_disease(data_wide = data_preg, 
                 metadata = meta_preg,
                 disease = "Pregnancy", color = UCAN_outlier
                 correct = F,
                 cutoff = 0.5)

plot_volcano(de_preg)
ggsave(savepath("volcano_pregnancy.png"), h = 5, w = 5)

# Look at proteins
# top_down <- 
#   de_preg |>
#   arrange(logFC) |> 
#   head(5)
# 
# top_up <- 
#   de_preg |>
#   arrange(-logFC) |> 
#   head(5)


top <- 
  de_preg |> 
  head(12)

# Only women
data_preg |> 
  pivot_longer(cols = -DAid, names_to = "Assay", values_to = "NPX") |> 
  filter(Assay %in% top$Assay) |>
  left_join(manifest) |> 
  mutate(Assay = factor(Assay, levels = top$Assay) ) |>
  ggplot(aes(Cohort, NPX, color = Cohort)) +
  geom_quasirandom() +
  geom_boxplot() +
  facet_wrap(~Assay) +
  theme_hpa(angled = T)

high_csh1 <- 
  data_preg |> 
  filter(CSH1 > 5)
data_preg |> 
  pivot_longer(cols = -DAid, names_to = "Assay", values_to = "NPX") |> 
  filter(Assay %in% top$Assay) |>
  left_join(manifest) |> 
  mutate(Assay = factor(Assay, levels = top$Assay), 
         UCAN_outlier = ifelse(DAid == "DA14976", "Yes", "No")) |>
  ggplot(aes(Cohort, NPX)) +
  geom_quasirandom(aes(color = UCAN_outlier)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(~Assay) +
  theme_hpa(angled = T)

ggsave(savepath("pregnancy_proteins.png"), h = 8, w = 12)

data_preg |> 
  pivot_longer(cols = -DAid, names_to = "Assay", values_to = "NPX") |> 
  filter(Assay %in% top$Assay) |>
  left_join(manifest) |> 
  mutate(Assay = factor(Assay, levels = top$Assay), 
         CSH1_outlier = ifelse(DAid %in% high_csh1$DAid, "Yes", "No")) |>
  ggplot(aes(Cohort, NPX)) +
  geom_quasirandom(aes(color = CSH1_outlier)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(~Assay) +
  theme_hpa(angled = T)

ggsave(savepath("pregnancy_proteins_2.png"), h = 8, w = 12)


# UMAP B4

data_b4_w <- 
  data_b4 |> 
  select(DAid, Assay, NPX) |> 
  pivot_wider(names_from = Assay, values_from = NPX)

umap_final <- 
  do_umap(data = data_b4_w,
          wide = F,
          plots = F)

# umap_final |>
#   rename(DAid = Sample) |>
#   left_join(manifest, by = "DAid") |>
#   mutate(Cohort = paste(Cohort, `PI - sample owner`, sep = "_")) |>
#   ggplot(aes(UMAP1, UMAP2, color = Cohort)) +
#   geom_point(alpha = 0.7, size = 0.8) +
#   scale_color_manual(values = brewer.pal(n = 11, name = "Set3")) +  
#   theme_hpa() +
#   ggtitle("Final UMAP")

umap_final |> filter(UMAP1 < -4.3) |>
  rename(DAid = Sample) |>
  left_join(manifest)

umap_final |>
  rename(DAid = Sample) |>
  left_join(manifest, by = "DAid") |>
  mutate(Disease = ifelse(Cohort == "PSC1", Disease, NA)) |>
  ggplot(aes(UMAP1, UMAP2, color = Disease)) +
  geom_point(alpha = 0.7, size = 0.8) +
  #scale_color_manual(values = brewer.pal(n = 11, name = "Set3")) +
  theme_hpa() +
  ggtitle("Final UMAP")

umap_final |>
  rename(DAid = Sample) |>
  left_join(manifest, by = "DAid") |>
  filter(Cohort != "FTD2") |> 
  mutate(Disease = ifelse(Cohort == "PSC1", Disease, NA)) |>
  ggplot(aes(UMAP1, UMAP2, color = Disease)) +
  geom_point(alpha = 0.7, size = 0.8) +
  #scale_color_manual(values = brewer.pal(n = 11, name = "Set3")) +
  theme_hpa() +
  ggtitle("Final UMAP")

ggsave(savepath("final_umap_no_ftd_filter.png"), h = 6, w = 6)



#Remake
data_b4_no_FTD <- 
  data_b4 |> 
  left_join(manifest |> 
              select(DAid, Cohort), by = "DAid") |> 
  filter(Cohort != "FTD1") |> 
  select(DAid, Assay, NPX) |> 
           pivot_wider(names_from = Assay, values_from = NPX)
  
umap_final_noFTD <- 
  do_umap(data = data_b4_no_FTD,
          plots = F)

umap_final_noFTD <- umap_res

umap_final_noFTD |>
  rename(DAid = Sample) |>
  left_join(manifest, by = "DAid") |>
  mutate(Disease = ifelse(Cohort == "PSC1", Disease, NA)) |>
  ggplot(aes(UMAP1, UMAP2, color = Disease)) +
  geom_point(alpha = 0.7, size = 0.8) +
  #scale_color_manual(values = brewer.pal(n = 11, name = "Set3")) +
  theme_hpa() +
  ggtitle("Final UMAP")

ggsave(savepath("final_umap_no_ftd_recalc.png"), h = 6, w = 8)

         