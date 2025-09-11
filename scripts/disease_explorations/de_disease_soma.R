# SOMA x
manifest_soma <- 
  manifest |> 
        filter(!(Disease == "Lung cancer" & Cohort == "ANMU")) |> 

    filter(DAid %in% data_soma_all$DAid) |> 
    mutate(Disease = case_when(Disease == "Pancreatic ductal adenocarcinoma" ~ "Pancreatic cancer", 
  Diagnose == "OVC" ~ "Ovarian cancer", 
Diagnose == "BRC" ~ "Breast cancer",
Diagnose == "LUNG" ~ "Lung cancer",
Diagnose == "CRC" ~ "Colorectal cancer",
Diagnose == "PRC" ~ "Prostate cancer",
T ~ Disease),
Disease = ifelse(Cohort == "EPIL", "Epilepsy", Disease),
Disease = ifelse(Cohort %in% c("FIBR", "PARD"), Diagnose, Disease),
Disease = ifelse(Disease %in% c("Healthy control", "healthy", "Control"), "Healthy", Disease))  |> 
  filter(!is.na(Disease))

other_controls_soma <- 
  manifest_soma |> 
  filter(!Disease %in% all_cancers,
  !Disease %in% c("Control", "Healthy", "Healthy control", "No cognitive impairment", "No diagnosis", "healthy"),
!is.na(Disease)) |>
  distinct(DAid)

# Healthy comparison
meta_healthy_soma <- 
  manifest_soma |> 
  filter(Disease %in% c(disease, "Healthy")) |> 
  select(DAid, Age, Sex, BMI, Disease)

data_healthy_soma <- 
  data_soma_all |> 
  select(-EntrezGeneSymbol, -Target) |> 
  filter(DAid %in% meta_healthy_soma$DAid) |> 
  mutate(RFU = log2(RFU)) |> 
  pivot_wider(names_from = Assay, values_from = RFU)

de_healthy_soma <-     
    do_limma_disease(data_wide = data_healthy_soma,
                     metadata = meta_healthy_soma,
                     disease = disease,
                     controls = "Healthy",
                     correct = F,
                     cutoff = 0) |> 
  translate_soma()

volcano_healthy_soma <- plot_volcano(de_healthy_soma, multiple = T, nrow = 2, cutoff = 0) + ggtitle("Healthy")


# Pan-cancer comparison
meta_cancer_soma <- 
  manifest_soma |> 
  filter(Disease %in% all_cancers) |> 
  select(DAid, Age, Sex, BMI, Disease)

data_cancer_soma <- 
   data_soma_all |> 
  select(-EntrezGeneSymbol, -Target) |> 
  filter(DAid %in% meta_cancer_soma$DAid) |> 
  mutate(RFU = log2(RFU)) |> 
  pivot_wider(names_from = Assay, values_from = RFU)

de_cancer_soma<-     
    do_limma_disease(data_wide = data_cancer_soma,
                     metadata = meta_cancer_soma,
                     disease = disease,
                     controls = "Cancer",
                     correct = F,
                     cutoff = 0) |> 
  translate_soma()

volcano_cancer_soma <- plot_volcano(de_cancer_soma, multiple = T, nrow = 2, cutoff = 0) + ggtitle("Pan-cancer")

# Pan-disease comparison
meta_disease_soma <- 
  manifest_soma |> 
  filter(DAid %in% unique(other_controls_soma$DAid) | Disease == disease) |> 
  select(DAid, Age, Sex, BMI, Disease)

data_disease_soma <- 
  data_soma_all |> 
  select(-EntrezGeneSymbol, -Target) |> 
  filter(DAid %in% meta_disease_soma$DAid) |> 
  mutate(RFU = log2(RFU)) |> 
  pivot_wider(names_from = Assay, values_from = RFU)

de_disease_soma <-     
    do_limma_disease(data_wide = data_disease_soma,
                     metadata = meta_disease_soma,
                     disease = disease,
                     controls = "Disease",
                     correct = F,
                     cutoff = 0) |> 
  translate_soma()

volcano_disease_soma <- plot_volcano(de_disease_soma, multiple = T, nrow = 2, cutoff = 0) + ggtitle("Pan-disease")


volcano_healthy_soma + volcano_cancer_soma + volcano_disease_soma
ggsave(paste0("results/de_", disease, "_Soma.png"), h = 5, w = 12)


# Examples
de_combined_soma<- 
  de_healthy_soma |> 
  bind_rows(de_cancer_soma) |> 
  bind_rows(de_disease_soma)

top_assays <- 
  de_combined_soma |> 
  group_by(Control) |> 
  top_n(20, abs(logFC)) |> 
  ungroup() |> 
  distinct(Assay)

my_palette <- colorRampPalette(c("white", "#DDC49C", "red"))(100)

de_combined_soma |>
  filter(Assay %in% top_assays$Assay) |> 
  group_by(Assay, Control) |> 
  top_n(1, logFC) |> 
  select(Assay, logFC, Control) |> 
  pivot_wider(names_from = Control, values_from = logFC) |> 
  column_to_rownames("Assay") |> 
  pheatmap() |> 
  as.ggplot()

ggsave(paste0("results/de_pheatmap_", disease, "_Soma.png"), h = 8, w = 5)


top_5 <- 
  de_combined_soma |>
  filter(Assay %in% top_assays$Assay) |>
   group_by(Assay, Control) |> 
  top_n(1, logFC) |> 
  select(Assay, logFC, Control) |> 
  pivot_wider(names_from = Control, values_from = logFC) |> 
  mutate(rank_healthy = rank(-Healthy),
rank_disease = rank(-Disease),
rank_cancer = rank(-Cancer)) |> 
  group_by_all() |> 
mutate(sum_rank = sum(rank_healthy, rank_cancer, rank_disease)) |> 
  arrange(sum_rank) |> 
  head(5)


cancer_controls <- 
  manifest_soma |> 
  filter(Disease %in% all_cancers,
  Disease != disease) |> 
  distinct(Disease) |> 
  pull()

disease_controls <- 
  manifest_soma |> 
  filter(DAid %in% other_controls_olink$DAid) |> 
  distinct(Disease) |> 
  pull()


levels <- c(disease, "Healthy", cancer_controls, disease_controls)
pal_controls <- c(disease = "darkred", "Healthy" = "grey", "Pan-cancer controls" = "#DDC49C", "Pan-disease controls" = "#99B2A6")

manifest_soma_controls <- 
  manifest_soma |> 
  mutate(Type = case_when(Disease == disease ~ Disease,
  Disease == "Healthy" ~ Disease,
Disease %in% cancer_controls ~ "Pan-cancer controls",
Disease %in% disease_controls ~ "Pan-disease controls")) |> 
  mutate(Disease = factor(Disease, levels = c(disease, "Healthy", cancer_controls, disease_controls)))

data_soma_all |>
      filter(EntrezGeneSymbol %in% top_5$Assay,
      DAid %in% manifest_soma$DAid) |>
      select(DAid, EntrezGeneSymbol, (RFU)) |>
      left_join(manifest_soma_controls |>
                  select(DAid, Disease, Type),
                by = "DAid") |>
    mutate(Disease = factor(Disease, levels = c(disease, "Healthy", cancer_controls, disease_controls))) |> 
filter(!is.na(Disease)) |> 
        ggplot(aes(Disease,
                   log2(RFU),
                   fill = Type,
                   color = Type)) +
        geom_quasirandom(alpha = 0.7
) +
        geom_boxplot(
          alpha = 0.3,
          outlier.color = NA,
          color = "grey20",
          show.legend = F
        ) +
 scale_color_manual(
    values = pal_controls,
    na.value = "darkred")  +
 scale_fill_manual(
    values = pal_controls,
    na.value = "darkred") +  
  facet_wrap(~EntrezGeneSymbol, scales = "free_y", ncol = 1) +
        theme_hpa(angled = T) +
        xlab("") 

ggsave(paste0("results/de_examples_", disease, "_Soma.png"), h = 10, w = 10)
