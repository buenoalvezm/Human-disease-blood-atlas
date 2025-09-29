library(tidyverse)
library(readxl)
setwd( "/Users/mariabueno/Library/CloudStorage/OneDrive-KTH/Repos/Human-disease-blood-atlas/")

data_b1 <- read_tsv("data/final_data/data_phase2_batch1_curated_20241217.tsv") 
data_b2 <-  read_tsv(file = "data/final_data/data_phase2_batch2_curated_20250925.tsv")
data_b3 <- read_tsv("data/final_data/data_phase2_batch2_curated_20250127.tsv") 
data_b4 <-  read_tsv("data/final_data/data_phase2_batch4_curated_20250425.tsv") 
data_b5 <-  read_tsv(file = "data/final_data/data_phase2_batch5_curated_20250701.tsv")

manifest <- read_excel("../metadata/samples_2025-05-12.xlsx")

 
# Combine
combined_data <- 
  data_b1 |> 
  select(DAid, Assay, OlinkID, NPX) |> 
  mutate(Batch = "1") |> 
  bind_rows(data_b2 |> 
              select(DAid, Assay, OlinkID, NPX) |> 
              mutate(Batch = "2")) |> 
  bind_rows(data_b3 |> 
              select(DAid, Assay, OlinkID, NPX) |> 
              mutate(Batch = "3")) |> 
  bind_rows(data_b4 |> 
              select(DAid, Assay, OlinkID, NPX) |> 
              mutate(Batch = "4")) |> 
   bind_rows(data_b5 |> 
              select(DAid, Assay, OlinkID, NPX) |> 
              mutate(Batch = "5")) |> 
  mutate(DAid_batch = paste(DAid, Batch, sep = "_"))

write_tsv(combined_data, "../olink_data/phase2_unbridged_data_HT_20250926.tsv")


# Soma data
soma_data <- read_tsv(file = "../soma_data/data_phase2_somalogic_curated_20250922.tsv") 


# Prepare metadata (overlapping soma)
manifest_common <- 
  manifest |> 
#        filter(!(Disease == "Lung cancer" & Cohort == "ANMU")) |> 

    filter(DAid %in% unique(soma_data$DAid)) |> 
    mutate(Disease = case_when(Disease == "Pancreatic ductal adenocarcinoma" ~ "Pancreatic cancer", 
  Diagnose == "OVC" ~ "Ovarian cancer", 
Diagnose == "BRC" ~ "Breast cancer",
Diagnose == "LUNG" ~ "Lung cancer",
Diagnose == "CRC" ~ "Colorectal cancer",
Diagnose == "PRC" ~ "Prostate cancer",
Diagnose == "4_random" ~ "BAMSE - visit 4",
Diagnose == "8_random" ~ "BAMSE - visit 8",
Diagnose == "16_random" ~ "BAMSE - visit 16",
Diagnose == "24_random" ~ "BAMSE - visit 24",
Diagnose == "pregnancy" ~ "Pregnancy",
T ~ Disease),
Disease = ifelse(Cohort == "EPIL", "Epilepsy", Disease),
Disease = ifelse(Cohort %in% c("FIBR", "PARD"), Diagnose, Disease),
Disease = ifelse(Disease %in% c("Healthy control", "healthy", "Control"), "Healthy", Disease),
Class = case_when(Disease %in% c("Breast cancer", "Colorectal cancer", "Lung cancer", "Ovarian cancer", "Prostate cancer") ~ "Cancer",
Disease == "PD" ~ "Neurologic",
Cohort == "BAMS" ~ "Healthy",
Cohort == "PREG" ~ "Healthy",
Disease == "Healthy" ~ "Healthy",
T ~ Class))  

manifest_common |> 
  filter(is.na(Disease)) |> count(Cohort, Diagnose)

# Filter data 
proteins_block8 <- 
  data_b1 |> 
  distinct(Assay, Block) |> 
  filter(Block == 8)

data_olink_unbridged <- 
  combined_data |> 
  filter(DAid %in% manifest_common$DAid,
        !Assay %in% proteins_block8$Assay)


manifest_overlapping <- 
  manifest_common |> 
  select(DAid, Cohort, Age, Sex, BMI, Disease, Class)

manifest_overlapping |> filter(is.na((Class))) |> count(Cohort, Disease)


write_tsv(manifest_overlapping, "../olink_data/metadata_olink_soma_phase2.tsv")

# Filter data
write_tsv(data_olink_unbridged, "../olink_data/phase2_unbridged_data_overlapping_HT_pandisease_20250926.tsv")

manifest_overlapping |> distinct(DAid) |> nrow()
data_olink_unbridged |> distinct(DAid) |> nrow()
soma_data|> distinct(DAid) |> nrow()

manifest_overlapping |> 
  filter(!DAid %in% data_olink_unbridged$DAid) |> 
  count(Disease)
