library(tidyverse)
library(OlinkAnalyze)
library(umap)
library(readxl)


data_curated <-  read_tsv("data/final_data/data_phase2_batch4_curated_20250425.tsv") 
data_raw <-  read_tsv("data/final_data/data_phase2_batch4_raw_20240425.tsv") 
data_ht <- read_NPX("data/VL 3530B4 Extended NPX Mar 25 2025.parquet")

data <- read_NPX("data/final_data/data_phase2_batch4_curated_20250425.tsv") 

manifest <- read_excel("data/samples_2025-04-07.xlsx")

data_ht %>% 
  dplyr::filter(!stringr::str_detect(SampleID, 'CONTROL_SAMPLE')) %>% 
  olink_umap_plot(df = .,
                 color_g = "QC_Warning", byPanel = TRUE)  

npx_df <- 
  data_curated |>
  left_join(manifest |> 
              select(DAid, Disease), by = "DAid") |> 
  rename(SampleID = DAid) |> 
  dplyr::filter(!grepl("control", SampleID, ignore.case = TRUE))

disease <- 
  npx_df |> 
  filter(Disease %in% c("Healthy", "Pancreatic ductal adenocarcinoma"))

ttest_results <- olink_ttest(
  df = disease,
  variable = "Disease",
  alternative = "two.sided")

  gsea_results <- olink_pathway_enrichment(data = disease, test_results = ttest_results)
  
  ora_results <- olink_pathway_enrichment(
    data = npx_data1,
    test_results = ttest_results, method = "ORA")
}, silent = TRUE)


From OlinkAnalyze tests:
- Not possible to import with read_NPX because it is not the original format, but works with read_tsv of course
- Need to rename the DAid column to SampleID
- Some functions like UMAP fail but also fails for the raw data?
- Managed to run t-test -> so I think we are good? 
