library(tidyverse)
library(readxl)

# Starting LIMS file
warehouse_meta_raw <- readxl::read_excel("data/Disease LIMS samples July 15.xlsx")

# Diseases to exclude from resource
exclude_disease <- c("Influenza follow up",
                     "Influenza healthy",
                     "Lean healthy",
                     "Liver healthy",
                     "Malaria healthy",
                     "Teenager healthy",
                     "Wellness",
                     "CAD other", 
                     "Rare blood cancers",
                     "Other tek", 
                     "Healthy Turkey",
                     "EpiHealth",
                     "NAFLD AT",
                     "Colorectal cancer AT",
                     "Obesity Turkey",
                     "Stroke", 
                     "Myocardial infarction", 
                     "Coronary artery interv",
                     "Malaria 10D", 
                     "Pediatric tumor other", 
                     "Langerhans cellhistiocytos",
                     "Occlusion or stenosis of the carotid artery",
                     "HIV_followup",
                     "Rheumatoid arthritis (compromised)",
                     "Previous venous thromboembolism (controls)",
                     "Acute coronary syndrome (controls)",
                     "Acute venous thromboembolism (controls)",
                     "Technical control",
                     "UCAN controls",
                     "NAFLD",
                     "IGST/SCAPIS other")

# Exclude UCAN samples based on UCAN file
ucan_meta <- read_tsv("data/clinical_metadata/20230504metadata_UCAN.tsv")

ucan_samples_exclude <- 
  warehouse_meta_raw |> 
  filter(Cohort == "UCAN",
         Disease != "EpiHealth") |> 
  filter(!DAid %in% ucan_meta$DAid) |> 
  pull(DAid)

ucan_diseases <- 
  warehouse_meta_raw |> 
  distinct(Disease, Cohort) |> 
  filter(Cohort == "UCAN",
         !Disease %in% c("UCAN controls", "EpiHealth")) |> 
  pull(Disease)

# Exclusions CVD/Infectious based on Emil's files
exclusions_cvd <- read_csv("data/exclusions/resource_exclusion_list_CVD.csv")
exclusions_infectious <- read_csv("data/exclusions/resource_exclusion_list_INFECT.csv")


# Regrouping of liver diseases based on Ozlem/Ligqi's files
liver_meta <- read_excel("data/clinical_metadata/_Lingqi Meng_Plasma Proteomics Clinical Data.v2.xlsx")

arld_samples <- 
  warehouse_meta_raw  |>
  filter(Disease == "Fatty Liver Disease") |> 
  select(DAid) |> 
  left_join(liver_meta, by = "DAid")  |> 
  filter(Diagnosis == "ARLD")
  
masld_samples <- 
  warehouse_meta_raw  |>
  filter(Disease == "Fatty Liver Disease") |> 
  select(DAid) |> 
  left_join(liver_meta, by = "DAid")  |> 
  filter(Diagnosis == "MASLD")

# Exclude Elite controllers based on Manos' analyses
elite_controllers <- 
  warehouse_meta_raw |> 
  filter(Disease == "HIV_baseline",
         Subcategory == "Elite Controller") |> 
  pull(DAid)

# Filter data and modify Disease abd Class banes
resource_meta <- 
  warehouse_meta_raw |> 
  mutate(Disease = ifelse(Cohort %in% c("IGTS", "SCAP") & is.na(Disease), "IGST/SCAPIS other", Disease),
         Class = ifelse(Cohort %in% c("IGTS", "SCAP") & is.na(Disease), "Healthy", Class),
         Class = ifelse(Disease %in% ucan_diseases, "Cancer", Class), 
         Disease = case_when(Subcategory == "Neuroblastoma" ~ "Pediatric neuroblastoma",
                             Subcategory == "Retinoblastoma" ~ "Pediatric retinoblastoma",
                             Subcategory == "Wilms tumor" ~ "Pediatric kidney tumor",
                             Subcategory %in% c("Pilocytic astrocytoma", "Diffuse midline glioma (incl DIPG)", "Glioblastoma multiforme") ~ "Pediatric diffuse astrocytic and oligodendro. tumor",
                             Subcategory %in% c("Medulloblastoma", "Ganglioglioma", "Anaplastic Ependymoma", "Meningioma", "Optical pathway glioma in NF-1" ) ~ "Pediatric CNS tumor",
                             Subcategory %in% c("Ewing sarcoma", "Osteosarcoma") ~ "Pediatric bone tumor",
                             Subcategory == "Long Covid" & Class == "Pediatric" ~ "Pediatric long COVID",
                             Disease == "Pediatric Sarcoma" ~ "Pediatric sarcoma",
                             Disease == "Pediatric Lymphoma" ~ "Pediatric lymphoma",
                             T ~ Disease)) |> 
  filter(!is.na(Disease),
         Disease != "NA",
         Cohort != "TECH",
         !Disease %in% c(exclude_disease, "HIV"),
         !DAid %in% elite_controllers,
         !DAid %in% c(ucan_samples_exclude, exclusions_cvd$DAid, exclusions_infectious$DAid),
         !Subcategory %in% c("Mild Covid", "Recovered COVID")#,
         #DAid %in% resource_data$DAid
  ) |> 
  mutate(Class = case_when(Disease == "Coronary artery calcification" ~ "Cardiovascular",
                           Disease == "Previous venous thromboembolism" ~ "Cardiovascular",
                           Disease == "Rheumatoid arthritis" ~ "Autoimmune",
                           Disease == "Acute coronary syndrome" ~ "Cardiovascular",
                           Disease == "Acute venous thromboembolism" ~ "Cardiovascular",
                           T ~ Class),
         
         Disease = ifelse(DAid %in% arld_samples$DAid, "Alcohol-related liver disease", Disease),
         Disease = ifelse(DAid %in% masld_samples$DAid, "MASLD", Disease),
         Disease = case_when(Disease == "HIV_baseline" ~ "HIV", 
                             Disease == "Long COVID" ~ "Pediatric long COVID", 
                             Disease == "t2d" ~ "Type 2 diabetes",
                             Disease == "MetS" ~ "Metabolic syndrome",
                             Disease == "NAFLD" ~ "Nonalcoholic fatty liver disease",
                             Disease == "Chronic Liver Disease (CLD)" ~ "Chronic liver disease",
                             Disease == "Breast cancer DCIS" ~ "Breast ductal carcinoma in situ",
                             Disease == "Systemisk Lupus Erythematosus" ~ "Systemic lupus erythematosus",
                             Disease == "Pneumococcal Pneumonia" ~ "Pneumococcal pneumonia",
                             Disease == "Pediatric Sarcoma" ~ "Pediatric sarcoma",
                             Disease == "Pediatric Lymphoma" ~ "Pediatric lymphoma",
                             Disease == "Pediatric Blastoma" ~ "Pediatric blastoma",
                             Disease == "Pediatric Glioma" ~ "Pediatric glioma",
                             Disease == "Pediatric Brain tumor other" ~ "Pediatric brain tumor other",
                             T ~ Disease),
         Age = ifelse(Disease == "Pediatric long COVID", round(Age/12, 0), Age),
         Class = ifelse(Disease == "Pediatric long COVID", "Pediatric", Class)) |>
  filter(Disease != "Fatty Liver Disease")


# Save final data
write_csv(resource_meta, "data/processed/final_data/metadata_resource_20240715.csv") 
