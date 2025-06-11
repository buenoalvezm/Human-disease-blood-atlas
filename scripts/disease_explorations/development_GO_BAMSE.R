 
library(tidyverse)

source("scripts/functions/functions_utility.R")

go <- read_tsv("../Human-blood-atlas/data/hpa/Ensembl103 GO terms.txt")

go |> 
  distinct(`GO term name`) |> 
  filter(grepl("development", `GO term name`, ignore.case = TRUE)) |> 
  write_tsv("go_terms_development.tsv")

675
