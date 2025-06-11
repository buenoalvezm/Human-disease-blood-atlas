
library(tidyverse)
source("scripts/functions/themes_palettes.R")
source("scripts/functions/functions_visualization.R")
source("scripts/functions/functions_utility.R")


# Read all curated files
b1 <- read_tsv(file = "data/final_data/data_phase2_batch1_curated_20241217.tsv") 
b3 <- read_tsv(file = "data/final_data/data_phase2_batch2_curated_20250127.tsv") 
b4 <- read_tsv(file = "data/final_data/data_phase2_batch4_curated_20250425.tsv") 
b5 <- read_tsv(file = "data/final_data/data_phase2_batch5_curated_20250512.tsv") 
b2 <- read_tsv(file = "data/final_data/data_phase2_batch2_curated_20250512.tsv") 

manifest <- import_df("data/samples_2025-05-12.xlsx")

# Make overview plots


# Make UMAP plots
umap_b1 <- plot_umap(data = b1)
umap_b2 <- plot_umap(data = b2)
umap_b3 <- plot_umap(data = b3)
umap_b4 <- plot_umap(data = b4)
umap_b5 <- plot_umap(data = b5)

saveRDS(umap_b1, savepath_data("UMAP", "umap_b1.rds"))
saveRDS(umap_b2, savepath_data("UMAP", "umap_b2.rds"))
saveRDS(umap_b3, savepath_data("UMAP", "umap_b3.rds"))
saveRDS(umap_b4, savepath_data("UMAP", "umap_b4.rds"))
saveRDS(umap_b5, savepath_data("UMAP", "umap_b5.rds"))

# Make barplots
samples_b1 <- plot_barplot(data = b1, 
                           batch = "Batch 1")

samples_b2 <- plot_barplot(data = b2, 
                           batch = "Batch 2")

samples_b3 <- plot_barplot(data = b3, 
                           batch = "Batch 3")

samples_b4 <- plot_barplot(data = b4, 
                           batch = "Batch 4")

samples_b5 <- plot_barplot(data = b5, 
                           batch = "Batch 5")


samples_b1 + umap_b1$umap_plot 
ggsave(savepath("b1_summary.png"), h = 6, w = 12)

samples_b2 + umap_b2$umap_plot 
ggsave(savepath("b2_summary.png"), h = 6, w = 12)

samples_b3 + umap_b3$umap_plot 
ggsave(savepath("b3_summary.png"), h = 6, w = 12)

samples_b4 + umap_b4$umap_plot 
ggsave(savepath("b4_summary.png"), h = 6, w = 12)

samples_b5 + umap_b5$umap_plot 
ggsave(savepath("b5_summary.png"), h = 6, w = 12)


# Disease specific for B3 (UCAN)
b3 |> 
  distinct(DAid) |> 
  left_join(manifest, by = "DAid") |> 
  count(Diagnose) |> 
  ggplot(aes(Diagnose, n, fill = Diagnose)) +
  geom_col() +
  scale_fill_manual(values = pal_ucan) +  
  geom_text(aes(label = n), vjust = -0.5) +
  theme_hpa(angled = T) +
  ylab("Number of samples") +
  ggtitle("Batch 3") +

umap_b3$umap |> 
  rename(DAid = Sample) |> 
  left_join(manifest, by = "DAid") |>
  ggplot(aes(UMAP1, UMAP2, color = Diagnose)) +
  geom_point(alpha = 0.7, size = 0.8) +
  scale_color_manual(values = pal_ucan) +  
  theme_hpa()

ggsave(savepath("b3_summary_disease.png"), h = 6, w = 12)

# Disease specific for B2 (Turkey)
b2 |> 
  distinct(DAid) |> 
  left_join(manifest, by = "DAid") |> 
  count(Disease) |> 
  ggplot(aes(Disease, n, fill = Disease)) +
  geom_col(show.legend = F) +
 # scale_fill_manual(values = pal_ucan) +  
  geom_text(aes(label = n), vjust = -0.5) +
  theme_hpa(angled = T) +
  ylab("Number of samples") +
  ggtitle("Batch 2") +
  
  umap_b2$umap |> 
  rename(DAid = Sample) |> 
  left_join(manifest, by = "DAid") |>
  ggplot(aes(UMAP1, UMAP2, color = Disease)) +
  geom_point(alpha = 0.7, size = 0.8, show.legend = F) +
  #scale_color_manual(values = pal_ucan) +  
  theme_hpa()

ggsave(savepath("b2_summary_disease.png"), h = 8, w = 12)




plot_umap <- function(data) {
  
  input_data <- 
    data |>
    select(DAid, OlinkID, NPX) |>
    pivot_wider(names_from = "OlinkID", values_from = "NPX")
  
  umap_all_proteins <-
    do_umap(data = input_data,
            plots = F)
  
  umap_plot <- 
    umap_all_proteins |>
    rename(DAid = Sample) |>
    left_join(manifest, by = "DAid") |>
    #mutate(Cohort = paste(Cohort, `PI - sample owner`, sep = "_")) |>
    ggplot(aes(UMAP1, UMAP2, color = Cohort)) +
    geom_point(alpha = 0.7, size = 0.8) +
    scale_color_manual(values = brewer.pal(n = 11, name = "Set3")) +  
    theme_hpa()
  
  return(list(umap_plot = umap_plot, 
              umap = umap_all_proteins))
  
}

plot_barplot <- function(data, 
                         batch) { 
  data |> 
    distinct(DAid) |> 
    left_join(manifest, by = "DAid") |> 
    mutate(Cohort = paste(Cohort, `PI - sample owner`, sep = "_")) |>
    count(Cohort) |> 
    ggplot(aes(Cohort, n, fill = Cohort)) +
    geom_col() +
    scale_fill_manual(values = brewer.pal(n = 11, name = "Set3")) +  
    geom_text(aes(label = n), vjust = -0.5) +
    theme_hpa(angled = T) +
    ylab("Number of samples") +
    ggtitle(batch) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
}


# Batch 2 outliers
plate_layout <- read_tsv("data/Turkey_layout_cohort_plates_for_lims_2024-02-20.txt")
umap_b2 <- readRDS(savepath_data("UMAP", "umap_b2.rds"))

umap_b2 |> 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point()

# Read the Excel file
annotation_plates <- read_excel("data/annotation_plates.xlsx")

# Generate row-wise and column-wise well numbers
rows <- LETTERS[1:8]
cols <- 1:12

# Well names in row-wise and column-wise order
rowwise_wells <- paste0(rep(rows, each=12), rep(cols, times=8))   # A1, A2, ..., H12
colwise_wells <- paste0(rep(rows, times=12), rep(cols, each=8))   # A1, B1, ..., H12

# Create mapping: rowwise number (1:96) -> colwise number (1:96)
mapping <- setNames(match(rowwise_wells, colwise_wells), 1:96)

# Apply the mapping to convert to column-wise numbering
annotation_plates$position_columnwise <- mapping[as.character(annotation_plates$position)]

# Save the result
write_xlsx(annotation_plates, "annotation_plates_columnwise.xlsx")




umap_b2$umap |> 
  left_join(plate_layout, by = c("Sample" = "da_id")) |> 
  left_join(annotation_plates |> 
              select(-position) |> 
              rename(position = position_columnwise), by = c("plate_label", "position")) |> 
  #filter(plate_label == "Disease Atlas Phase 2 Batch 2 Olink P11/20") |> 
  #count(plate_label) |> 
  #arrange(n)
  select(Sample, UMAP1, UMAP2, plate_label, position, state) |> 
  ggplot(aes(UMAP1, UMAP2, color = state)) +
  geom_point() +
  theme_hpa()
