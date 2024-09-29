
library(tidyverse)

# Read in data
resource_meta <- read_csv("data/processed/final_data/metadata_resource_20240715.csv") 
resource_data <- read_csv(file = "data/processed/final_data/data_resource_20240604.csv") 


# Save metadata
resource_meta |> 
  filter(DAid %in% resource_data$DAid) |> 
  select(DAid, Class, Disease, Subcategory, Age, Sex, BMI) |> 
  mutate(Disease = ifelse(Disease == "Scleroderma", "Systemic sclerosis", Disease)) |>
  write_tsv("data/data_to_IT/da_phase1_metadata.tsv")

# Save data
resource_data |> 
  filter(DAid %in% resource_meta$DAid) |> 
  select(DAid, OlinkID, NPX, LOD) |> 
  write_tsv("data/data_to_IT/da_phase1_data.tsv")

# Save disease class order
class_order <- c("Healthy", "Cardiovascular","Metabolic","Cancer","Psychiatric","Autoimmune","Infection","Pediatric")
disease_class_order <-
  resource_meta |> 
  distinct(Class, Disease, Cohort) |> 
  mutate(Disease = ifelse(Disease == "Scleroderma", "Systemic sclerosis", Disease)) |>
  mutate(Class = factor(Class, levels = class_order)) |> 
  arrange(Class, Cohort) |> 
  pull(Disease) |>
  unique()

#saveRDS(pal_class, "../../DBA_data/Phase1_2024/pal_class.rds")
saveRDS(disease_class_order, "data/data_to_IT/disease_levels.rds")

# #Save differential expression results
combined_de <- readRDS("../Pan-disease-profiling/data/processed/DE_v5/combined_de.rds")

combined_de |> 
  filter(Disease == "Systemic sclerosis")

combined_de |> 
  filter(adj.P.Val == 0)

combined_de |> 
  arrange(-logFC)

combined_de|>
  left_join(resource_data |>
              distinct(Assay, OlinkID), by = "Assay") |>
  select(OlinkID, Assay, Disease, NPX_difference = logFC, p.adjusted = adj.P.Val, Significance = sig, Control) |>
  write_tsv("data/data_to_IT/differential_expression_data.tsv")

# #Save protein importance
# combined_importance <- readRDS(savepath_data("ML_v2", "combined_importance.rds"))
# 
# combined_importance |> 
#   group_by(Disease, Against) |> 
#   top_n(10, Importance)|> 
#   ungroup() |> 
#   left_join(resource_data |> 
#               distinct(Assay, OlinkID), by = "Assay") |> 
#   select(OlinkID, Assay, Disease, Overall_importance = Importance, Sign, Control = Against) |> 
#   mutate(Disease = case_when(Disease == "HIV_baseline" ~ "HIV", 
#                              Disease == "Long COVID" ~ "Pediatric long COVID", 
#                              Disease == "t2d" ~ "Type 2 diabetes",
#                              Disease == "MetS" ~ "Metabolic syndrome",
#                              Disease == "NAFLD" ~ "Nonalcoholic fatty liver disease",
#                              Disease == "Chronic Liver Disease (CLD)" ~ "Chronic liver disease",
#                              Disease == "Breast cancer DCIS" ~ "Breast ductal carcinoma in situ",
#                              Disease == "Pediatric Sarcoma" ~ "Pediatric sarcoma",
#                              Disease == "Pediatric Lymphoma" ~ "Pediatric lymphoma",
#                              Disease == "Pediatric Blastoma" ~ "Pediatric blastoma",
#                              Disease == "Pediatric Glioma" ~ "Pediatric glioma",
#                              Disease == "Pediatric Brain tumor other" ~ "Pediatric brain tumor other",
#                              T ~ Disease)) |>
#   write_tsv("../../DBA_data/Phase1_2024/disease_models_data.tsv")
