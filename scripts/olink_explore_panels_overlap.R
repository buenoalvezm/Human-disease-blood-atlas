library(tidyverse)
library(eulerr)
library(ggplotify)

source("scripts/functions_utility.R")
source("scripts/themes_palettes.R")

data_NPX_phase1 <- import_df(file = "data/final_data/final_data_phase1.csv") 

# Assays overlap
targets_ht <- 
  read_delim("data/olink_targets/targets_ht.csv") |> 
  distinct(Gene) |> 
  pull()

targets_3072 <- 
  import_df("data/olink_targets/targets_3k.xlsx") |> 
  distinct(`Gene name`) |> 
  pull()

targets_1463 <- 
  data_NPX_phase1 |> 
  distinct(Assay) |> 
  pull()



venn_list <- list("HT" = targets_ht, 
          "3K" = targets_3072,
          "1.5K" = targets_1463)

plot(euler(venn_list, shape = "ellipse"), quantities = TRUE, fills = c("#D69DCA","#D6EDDA", "#A7C7E7")) |> 
  as.ggplot() +
  ggtitle("Tarets from Olink Explore platforms") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 10)) 

ggsave(savepath("overlap_olink_platforms.png"), h = 5, w = 5)

targets_ht_extra <-  
  tibble(Assay = targets_ht) |> 
  filter(!Assay %in% c(targets_3072, targets_1463)) |> 
  pull()

targets_3072_extra <-  
  tibble(Assay = targets_3072) |> 
  filter(!Assay %in% targets_1463) |> 
  pull()

tibble(Assay = targets_ht) |> 
  mutate(Platform = case_when(Assay %in% targets_ht_extra ~ "HT",
                              Assay %in% targets_3072_extra ~ "3K",
                              Assay %in% targets_1463 ~ "1.5K")) |> 
  #count(Platform) |> 
  write_csv("data/olink_targets/overlap_olink_platforms.csv")

