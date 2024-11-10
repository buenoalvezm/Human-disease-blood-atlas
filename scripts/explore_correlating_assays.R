
library(tidyverse)

source("scripts/functions/functions_utility.R")
source("scripts/functions/functions_visualization.R")
source("scripts/functions/themes_palettes.R")

data_phase2 <- read_csv("data/final_data/data_phase2_20241110.csv")  
excluded_proteins <- read_csv("../Human-disease-blood-atlas/data/processed/hooking.correlation.above.0.7.csv")
manifest <- import_df("data/samples_2024-04-19.xlsx")


# Profiles
plots_proteins <- 
  lapply(excluded_proteins$value, function(protein) {
    plot_boxplot(protein,
                 data = data_phase2,
                 metadata = manifest,
                 platform = "HT",
                 title = protein)
  })

pdf(savepath("excluded_proteins_phase2.pdf"), h = 6, w = 12)
plots_proteins
dev.off()

# Correlations all proteins - annotate correlating assays

