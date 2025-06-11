library(tidyverse)
library(eulerr)
library(ggplotify)
library(OlinkAnalyze)
library(readxl)

raw_1463 <- 
  readr::read_delim("data/UD-2950_B1-B4_NPX_2023-03-29/UD-2950_B1-B4_NPX_2023-03-29.csv", delim = ";")

targets_1463 <- 
  raw_1463 |> 
  distinct(Assay)

# Assays overlap
raw_ht <- read_NPX("data/raw_data/VL-3530B3_NPX_2024-09-25.parquet")

targets_ht <- 
  raw_ht |> 
  filter(AssayType == "assay") |> 
  distinct(Assay)
  

targets_3072 <- 
  read_excel("data/olink_targets/targets_3k.xlsx") |> 
  distinct(Assay = `Gene name`)

venn_list <- list("HT" = targets_ht$Assay, 
                  "3K" = targets_3072$Assay,
                  "1.5K" = targets_1463$Assay)

plot(euler(venn_list, shape = "ellipse"), quantities = TRUE) |> #, fills = pal_platforms
  as.ggplot() +
  ggtitle("Tarets from Olink Explore platforms") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 10)) 

targets_ht |> 
  mutate(Platform = "Olink Explore HT",
         Platform = ifelse(Assay %in% targets_3072$Assay, "Olink Explore 3072", Platform),
         Platform = ifelse(Assay %in% targets_1463$Assay, "Olink Explore 1463", Platform)) |> 
  write_tsv("data/olink_targets/targets_olink_platforms.tsv")

ggsave(savepath("overlap_olink_platforms.png"), h = 5, w = 5)

targets_ht_extra <-  
  tibble(Assay = targets_ht) |> 
  filter(!Assay %in% c(targets_3072, targets_1463)) |> 
  pull()

targets_3072_extra <-  
  tibble(Assay = targets_3072) |> 
  filter(!Assay %in% targets_1463) |> 
  pull()

targets <- 
  tibble(Assay = targets_ht) |> 
  mutate(Platform = case_when(Assay %in% targets_ht_extra ~ "HT",
                              Assay %in% targets_3072_extra ~ "3K",
                              Assay %in% targets_1463 ~ "1.5K")) # |>  
#count(Platform)  

write_csv(targets, "data/olink_targets/overlap_olink_platforms.csv")