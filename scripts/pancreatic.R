library(tidyverse)
library(OlinkAnalyze)
library(arrow)
library(readxl)
library(ggsci)

source("scripts/functions/functions_utility.R")
source("scripts/functions/functions_visualization.R")

data_b4 <- read_NPX("data/VL 3530B4 NPX Dec 5 2024.parquet")
new_manifest <- read_excel("data/samples_2024-12-17.xlsx")

daid_patients_2_3 <- c("DA12362", "DA12363")

data <- 
  data_b4 |> 
  filter(SampleType == "SAMPLE",
         AssayType == "assay",
         !Assay %in% c("GBP1", "MAP2K1")) |> 
  mutate(DAid = str_extract(SampleID, "^[^-]+")) |> 
  filter(!DAid %in% daid_patients_2_3)

metadata <- 
  new_manifest |> 
  filter(DAid %in% data$DAid)

# Number of samples
new_manifest |> 
  filter(DAid %in% b4_id) |> 
  count(Cohort) |> 
  ggplot(aes(x = fct_reorder(Cohort, -n), y = n, fill = Cohort)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_simpsons() +
  xlab("Cohort") +
  theme_bw() 

ggsave("results/cohorts_b4.png", width = 6, height = 6)


# PCA
pca_b4 <- 
  do_pca(data = data,
         wide = F)

pca_b4$pca_res |> 
  left_join(metadata, by = c("Sample" = "DAid")) |>
  ggplot(aes(PC1, PC2, color = Cohort)) +
  geom_point() +
  scale_color_simpsons() +
  theme_hpa()

# UMAP
umap_b4 <- 
  do_umap(data = data,
          wide = F)

umap_b4 |> 
  left_join(metadata, by = c("Sample" = "DAid")) |>
  ggplot(aes(UMAP1, UMAP2, color = Cohort)) +
  geom_point() +
  scale_color_simpsons() +
  theme_hpa()

# Differential expression
disease_meta <- 
  metadata |> 
  mutate(Disease = Cohort) |> 
  mutate(Disease = ifelse(Disease == "PAC1", "Pancreatic cancer", "Control")) |> 
  select(DAid = DAid, Disease, Sex, Age, BMI)

disease_data <- 
  data |> 
  filter(DAid %in% disease_meta$DAid) |> 
  select(DAid, Assay, NPX) |> 
  pivot_wider(names_from = "Assay", values_from = "NPX")

library(limma)
de_pancreatic <- 
  do_limma_disease(data_wide = disease_data, 
                 metadata = disease_meta,
                 disease = "Pancreatic cancer",
                 correct = F,
                 cutoff = 0)

plot_volcano(de_pancreatic) +
  ggtitle("HT Phase 2 - Pancreatic cancer against other diseases")

ggsave(savepath("pancreatic_ht.png"))

de_pancreatic |> 
  filter(sig == "significant down") 


# phase 1
de_pancreatic_p1 <- 
  readRDS("../Pan-disease-profiling/data/processed/DE_v6/combined_de.rds") |> 
  filter(Control == "ALl other diseases",
         Disease == "Pancreatic cancer")


de_pancreatic_p1 |> 
  select(Assay, logFC_P1 = logFC, pval_P1 = adj.P.Val)  |> 
  left_join(de_pancreatic |> 
              select(Assay, logFC_P2 = logFC, pval_P2 = adj.P.Val), by = "Assay") |> 
  ggplot(aes(logFC_P1, logFC_P2)) +
  geom_point() +
  geom_text_repel(aes(label = Assay)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_hpa() +


de_pancreatic_p1 |> 
  select(Assay, logFC_P1 = logFC, pval_P1 = adj.P.Val)  |> 
  left_join(de_pancreatic |> 
              select(Assay, logFC_P2 = logFC, pval_P2 = adj.P.Val), by = "Assay") |> 
  ggplot(aes(-log10(pval_P1), -log10(pval_P2))) +
  geom_point() +
  geom_text_repel(aes(label = Assay)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_hpa()
 
ggsave(savepath("pancreatic_de_p1_p2.png"))
