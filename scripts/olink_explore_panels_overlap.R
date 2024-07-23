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

plot(euler(venn_list, shape = "ellipse"), quantities = TRUE, fills = pal_platforms) |> 
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

targets <- 
  tibble(Assay = targets_ht) |> 
  mutate(Platform = case_when(Assay %in% targets_ht_extra ~ "HT",
                              Assay %in% targets_3072_extra ~ "3K",
                              Assay %in% targets_1463 ~ "1.5K")) # |>  
  #count(Platform)  

write_csv(targets, "data/olink_targets/overlap_olink_platforms.csv")


hpa <- read_tsv("data/hpa/proteinatlas.tsv")
# HPA secretome
hpa |> colnames()

pal_secreted <- c("Secreted to blood" = "#B30000", 
                  "Secreted in brain" = "#FFDD00", 
                  "Secreted to digestive system" = "#1280C4", 
                  "Secreted in male reproductive system" = "#95D4F5", 
                  "Secreted in female reproductive system" = "#F8BDD7", 
                  "Secreted to extracellular matrix"  = "#7F6A9C", 
                  "Secreted in other tissues" = "#FFD480", 
                  "Secreted - unknown location" = "#A1A8AA", 
                  "Intracellular and membrane" = "#F9A266", 
                  "Unknown" = "grey80")

targets |> 
  left_join(hpa |> 
              select(Gene, `Secretome location`), by = c("Assay" = "Gene")) |> 
  group_by(Platform, `Secretome location`) |>
  count() |> 
  ggplot(aes(Platform, n, fill = `Secretome location`)) +
  geom_col() +
  scale_fill_manual(values = pal_secreted) +
  theme_hpa()

ggsave(savepath("olink_secretome.png"), h = 5, w = 6)

# HPA tissue specificity
targets |> 
  left_join(hpa |> 
              select(Gene, `Secretome location`), by = c("Assay" = "Gene")) |> 
  group_by(Platform, `Secretome location`) |>
  count() |> 
  ggplot(aes(Platform, n, fill = `Secretome location`)) +
  geom_col() +
  scale_fill_manual(values = pal_secreted) +
  theme_hpa()

ggsave(savepath("olink_secretome.png"), h = 5, w = 6)


targets |> 
  left_join(hpa |> 
              select(Gene, `Secretome location`), by = c("Assay" = "Gene")) |> 
  group_by(Platform, `Secretome location`) |>
  count() |> 
  filter(`Secretome location` == "Secreted to blood") 
