resource_data <- read_csv("data/data_resource_20240604.csv")
resource_meta <- read_csv("data/metadata_resource_20240715.csv")

# Filter to match data & metadata samples
resource_meta <- 
  resource_meta |> 
  filter(DAid %in% resource_data$DAid)

resource_data <- 
  resource_data |> 
  filter(DAid %in% resource_meta$DAid)

resource_data_w <- 
  resource_data |> 
  select(DAid, Assay, NPX) |> 
  pivot_wider(names_from = Assay, 
              values_from = NPX)




resource_meta  |> 
  filter(is.na(Sex) | is.na(Age)) |> 
  filter(is.na(Sex)) |> 
  select(DAid, `Vial barcode`, `PI - sample owner`, Cohort, Disease, Age, Sex) #|> 
# write_csv(savepath_data("LIMS", "samples_missing_age_sex.csv"))

original_missing <- 
  read_delim(savepath_data("LIMS", "samples_missing_age_sex_original.csv"))

new_missing <- 
  resource_meta  |> 
  filter(is.na(Sex) | is.na(Age)) |> 
  filter(is.na(Sex)) |> 
  select(DAid, `Vial barcode`, `PI - sample owner`, Cohort, Disease, Age, Sex) 

rescued <- 
  original_missing |> 
  filter(!DAid %in% new_missing$DAid) #|> 
#count(Disease)


resource_meta |> 
  filter(DAid %in% rescued$DAid) |> 
  select(DAid, Disease, Age, Sex)


check_sex <- 
  resource_meta |> 
  filter(Cohort == "SPS5")

female_up <-  c("CGA", "LEP", "LPL", "XG")
female_down <-  c("MMP3", "OBP2B","PROK1", "PSPN")

resource_data |> 
  filter(DAid %in% check_sex$DAid, 
         Assay %in% c(female_up, female_down)) |> 
  mutate(Type = case_when(Assay %in% female_up ~ "Upregulated in females",
                          Assay %in% female_down  ~ "Downregulated in females")) |> 
  left_join(resource_meta |> 
              select(DAid, Sex), by = "DAid") |> 
  ggplot(aes(Sex, NPX, fill = Sex, color = Sex)) +
  geom_point() +
  geom_boxplot(alpha = 0.5, outlier.color = NA) +
  facet_wrap(Type ~ Assay, scales = "free", ncol = 4) +
  scale_color_manual(values = pal_sex) +
  scale_fill_manual(values = pal_sex) +
  theme_hpa() +
  ggtitle("SPS5 samples")

ggsave(savepath("SPS5_sex.png"))