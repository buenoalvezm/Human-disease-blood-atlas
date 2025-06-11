
# Prepare data for Konstantinos project 
# 2024-09-27

library(tidyverse)

source("scripts/functions/functions_utility.R")

data_phase1 <- import_df("data/data_phase1_20240604.csv")
manifest <- import_df("data/samples_2024-09-27.xlsx")

data_igt_scapis <- 
  data_phase1 |> 
  left_join(manifest |> 
              distinct(DAid, Cohort), by = "DAid") |> 
  filter(Cohort %in% c("IGTS", "SCAP")) |> 
  select(-Cohort, -SampleID_Original) |> 
  relocate(DAid)


write_tsv(data_igt_scapis, "data/processed/LIMS/data_IGT_SCAPIS_preprocessed.csv")



data_phase1 |> 
  left_join(manifest |> 
              distinct(DAid, Cohort), by = "DAid") |> 
  filter(Cohort %in% c("IGTS", "SCAP")) 
  
  
  select(-Cohort, -SampleID_Original) |> 
  relocate(DAid)
  
  manifest |> 
    filter(Disease == "Healthy") |> 
    count(Cohort) 
  