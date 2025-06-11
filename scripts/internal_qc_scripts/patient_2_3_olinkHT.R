
library(tidyverse)
library(OlinkAnalyze)
library(arrow)
library(readxl)
library(ggsci)

data_b4 <- read_NPX("data/VL 3530B4 NPX Dec 5 2024.parquet")
data_b1_b3 <- read_NPX("data/cs-transfer/VL-3530B1_complete_NPX_2024-07-01.parquet")
manifest <- read_excel("data/samples_2024-04-19.xlsx")
new_manifest <- read_excel("data/samples_2024-12-17.xlsx")

daid_patients_2_3 <- c("DA12362", "DA12363")

# Look at samples run in batch 4 (first two kits)
b4_id <- 
  data_b4 |> 
  filter(SampleType == "SAMPLE") |> 
  mutate(DAid = str_extract(SampleID, "^[^-]+")) |>
  distinct(DAid) |> 
  pull()

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


# Save data for Olink
data_b4 |> 
  mutate(DAid = str_extract(SampleID, "^[^-]+")) |> 
  filter(DAid %in% daid_patients_2_3) |> 
  select(-DAid) |> 
  write_parquet("data/ht_b4_patients_2_3.parquet")

data_b1_b3 |> 
  mutate(DAid = str_extract(SampleID, "^[^-]+")) |> 
  filter(DAid %in% daid_patients_2_3) |> 
  select(-DAid) |> 
  write_parquet("data/ht_b1_b3_patients_2_3.parquet")

  
data_b1_b3 |> 
  mutate(Batch = str_extract(SampleID, "[^\\-]+$")) |> 
  distinct(Batch) |> View()

#Check samples in manifest

dat |> count(SampleID)

data_b4 |> 
  mutate(DAid = str_extract(SampleID, "^[^-]+")) |> 
  filter(DAid == "DA17568") |> 
  count(SampleID)

# CHeck if bridging samples
