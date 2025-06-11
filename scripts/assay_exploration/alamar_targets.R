# Overlap HT alamar

targets_alamar <- import_df("data/alamar.xlsx")

alamr_count <- 
  targets_alamar |> 
  count(`Function/Classification`) 

alamr_count_sorted <- alamr_count |>
  arrange(desc(n)) |>  # sort by abundance
  mutate(
    percent = round(100 * n / sum(n), 1),
    label = paste0(percent, "%")
  )

# Pie chart with center label
ggplot(alamr_count_sorted, aes(x = "", y = n, fill = `Function/Classification`)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center of each slice
    size = 4,
    color = "black"
  ) +
  labs(title = "Function Distribution") +
  scale_fill_futurama() +
  theme_void()

ggsave(savepath("alamar_pie.png"), width = 6, height = 6)


# Convert to percentage for labels (optional)
alamr_count$percent <- round(100 * alamr_count$n / sum(alamr_count$n), 1)
alamr_count$label <- paste0(alamr_count$`Function/Classification`, " (", alamr_count$percent, "%)")
alamr_count$ypos <- cumsum(alamr_count$n) - 0.5 * alamr_count$n  # position for labels

# Create pie chart
alamr_count |> 
  ggplot(aes(x = "", y = n, fill = `Function/Classification`)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
 # geom_text(aes(y = ypos, label = label), color = "black", size = 4) +
  labs(title = "Function Distribution") +
  scale_fill_futurama() +
  theme_void() 

targets_olink <- 
  data_b4 |> 
  distinct(Assay, AssayType) |> 
  filter(AssayType == "assay")


y <- list("Olink HT" = targets_olink$Assay, 
          "Alamar" = targets_alamar$Gene)
plot(euler(y, shape = "ellipse"), quantities = TRUE, fills = c("#D49DC5", "#A6CEE3")) |> as.ggplot()

ggsave(savepath("overlap_olink_alamar.png"), h = 4, w = 4)


targets_alamar |> 
  filter(!Gene %in% targets_olink$Assay) |> 
  select(Gene) |> 
  write_csv(savepath("alamar_olink_diff.csv"))

            