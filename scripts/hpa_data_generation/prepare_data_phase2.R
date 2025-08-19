library(tidyverse)

data_b1 <- read_tsv("data/final_data/data_phase2_batch1_curated_20241217.tsv") 
#data_b2 <-  read_tsv(file = "data/final_data/data_phase2_batch2_raw_20250512.tsv")
data_b3 <- read_tsv("data/final_data/data_phase2_batch2_curated_20250127.tsv") 
data_b4 <-  read_tsv("data/final_data/data_phase2_batch4_curated_20250425.tsv") 
data_b5 <-  read_tsv(file = "data/final_data/data_phase2_batch5_curated_20250701.tsv")

manifest <- readxl::read_excel("data/samples_2025-04-07.xlsx")

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
   bind_rows(data_b4 |> 
              select(DAid, Assay, NPX) |> 
              mutate(Batch = "5")) |> 
  mutate(DAid_batch = paste(DAid, Batch, sep = "_"))

write_tsv(combined_data, "../olink_data/phase2_unbridged_data.tsv")