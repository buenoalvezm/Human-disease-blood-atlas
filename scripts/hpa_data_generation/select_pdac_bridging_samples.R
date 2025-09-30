
# Bridging data

bridge_Samples<- npx_data1 %>% 
  # Excluding control samples. Naming convention may differ.
  dplyr::filter((stringr::str_detect(SampleID, "CONTROL", negate = TRUE))) %>%
  olink_bridgeselector(sampleMissingFreq = 0.1,
                     n = 16)

npx_data1 %>% 
  filter(!str_detect(SampleID, 'CONT')) %>%
  mutate(Bridge = ifelse(SampleID %in% bridge_Samples$SampleID, "Bridge", "Sample")) %>% 
  olink_pca_plot(color_g = "Bridge")






