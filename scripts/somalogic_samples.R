
source("scripts/functions/functions_utility.R")
source("scripts/functions/themes_palettes.R")

library(tidyverse)

data_soma <- import_df("data/Disease_Atlas_Somascan_selection.xlsx")

dat <- 
  data_soma |> 
  mutate(Disease = case_when(Disease == "Lung cancer" & Contact == "Adil Mardinoglu" ~ "Lung cancer (Adil)",
                             Disease == "Lung cancer" & Contact == "Fredrik Pontén" ~ "Lung cancer (UCAN)",
                             Disease == "Breast cancer" & Contact == "Adil Mardinoglu" ~ "Breast cancer (Adil)",
                             Disease == "Breast cancer" & Contact == "Fredrik Pontén" ~ "Breast cancer (UCAN)",
                             Disease == "Healthy" & Contact == "Adil Mardinoglu" ~ "Healthy (Adil)",
                             Disease == "Healthy" & Contact == "Camilla Svensson" ~ "Healthy (Camilla)",
                             Disease == "Healthy" & Contact == "Caroline Ingre" ~ "Healthy (Caroline)",
                             Disease == "Healthy" & Contact == "Lars Klareskog" ~ "Healthy (Lars)",
                             Disease == "Healthy" & Contact == "Per Svenningsson" ~ "Healthy (Per)",
                             T ~ Disease),
         Class = ifelse(Class == "Neuro", "Neurologic", Class),
         Class = ifelse(Class == "Infectious", "Infection", Class),
         Class = factor(Class, names(pal_class_2))) |> 
  arrange(Class, `Number of samples`) 

dat |> 
  mutate(Disease = factor(Disease, levels = dat$Disease)) |> 
  ggplot(aes(Disease, `Number of samples`, fill = Class)) +
  geom_col() +
  geom_text(aes(label = `Number of samples`, y = `Number of samples` + 7)) +
  scale_fill_manual(values = pal_class_2) +
  theme_hpa(angled = T) +
  xlab("") +
  ggtitle("Somalogic - sample overview")
  
ggsave(savepath("sample_overview_somalogic.png"), width = 14, height = 8)
 #c("#EAB8D1", "#F2A900", "#A3C4F3", "#A3EBA3")
