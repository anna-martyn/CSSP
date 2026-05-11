# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "dplyr", "tidyr", "tibble", "ggh4x")
for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
design <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
  header = TRUE,
  sep = "\t"
)
taxonomy <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)
asv_table <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  skip = 1,
  comment.char = ""
)

# Cleaning up taxonomy file
taxonomy <- taxonomy %>% rename(ASVid = Feature.ID)
taxonomy <- taxonomy %>%
  separate(
    col = Taxon,
    into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    sep = "; ",
    fill = "right"
  ) %>%
  mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Converting ASV reads to relative abundances (RA)
asv_table_norm <- sweep(asv_table, 2, colSums(asv_table), "/")
df <- as.data.frame(asv_table_norm) %>%
  rownames_to_column(var = "ASVid") %>%
  left_join(taxonomy %>% select(ASVid, Order), by = "ASVid")

# Reshaping dataframe to long format and keeping only soil (bulk) samples
df_long <- df %>%
  pivot_longer(cols = -c(ASVid, Order), names_to = "sampleID", values_to = "RA") %>%
  left_join(design %>% select(SampleID, Soil, Genotype), by = c("sampleID" = "SampleID")) %>%
  filter(Genotype == "Soil")

# Summarising RAs by order per sample
df_long_order <- df_long %>%
  group_by(Order, sampleID, Soil) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

# Identifying top 20 bacterial orders per soil type
top_orders <- df_long_order %>%
  group_by(Soil, Order) %>%
  summarise(MeanRA = mean(RA), .groups = "drop") %>%
  group_by(Soil) %>%
  slice_max(MeanRA, n = 20) %>%
  ungroup() %>%
  pull(Order) %>%
  unique()

# Grouping remaining orders as "Other" and sorting alphabetically
df_long_order <- df_long_order %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other")) %>%
  mutate(Order = factor(Order, levels = c(sort(unique(Order[Order != "Other"])), "Other"))) %>%
  mutate(Soil = factor(Soil, levels = c("NPK","PK","UF")))

# Order colours
colors <- read.table(
  file = "../../../0_files/Bacterial_order_colors.csv",
  header = TRUE,
  sep = ",",
  comment.char = ""
)

# Barplot of order-level RAs by sample ----------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  legend.position = "right",
  legend.background = element_blank(),
  text = element_text(family = "sans", size = 6, color = "black")
)

bar_plot <- ggplot(df_long_order, aes(x = sampleID, y = RA, fill = Order)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_nested(~Soil, scales = "free_x", labeller = label_parsed) +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand = c(0, 0)) +
  main_theme +
  ylab("Relative abundance") +
  labs(fill = "Bacterial order") +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.title.x = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1))

# Saving plot
ggsave(
  filename = "2_figures/Hordeum_bulk_order_top20_RA_stackedbp.pdf",
  plot = bar_plot,
  width = 8,
  height = 17,
  units = "cm"
)

saveRDS(
  object = bar_plot,
  file = "1_rds_files/Hordeum_bulk_order_top20_RA_stackedbp.rds"
)

# Barplot of order-level mean RAs by soil type --------------------------------

# Collapsing non-top20 orders as "Other"
df_long2 <- df_long %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other"))

# Aggregating RA at order-level
df_sample_order <- df_long2 %>%
  group_by(sampleID, Soil, Order) %>%
  summarise(RA = sum(RA), .groups = "drop")

# Mean RA per soil per order
df_mean_order <- df_sample_order %>%
  group_by(Soil, Order) %>%
  summarise(RA = mean(RA), .groups = "drop") %>%
  mutate(Order = factor(Order, levels = c(sort(unique(Order[Order != "Other"])), "Other")))

# Plot
barplot_mean <- ggplot(df_mean_order, aes(x = Soil, y = RA, fill = Order)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand = c(0,0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill = "Bacterial order") + 
  xlab("") +
  guides(fill = guide_legend(ncol = 1))

# Save plot
ggsave(
  filename = "2_figures/Hordeum_bulk_order_top20_RA_mean_stackedbp.pdf",
  plot = barplot_mean,
  width = 6,
  height = 17,
  units = "cm"
)
saveRDS(
  object = barplot_mean,
  file = "1_rds_files/Hordeum_bulk_order_top20_RA_mean_stackedbp.rds"
)

# Save plot without legend
barplot_mean_no_legend <- barplot_mean + theme(legend.position = "none")

ggsave(
  filename = "2_figures/Hordeum_bulk_order_top20_RA_mean_stackedbp_no_legend.pdf",
  plot = barplot_mean_no_legend,
  width = 6,
  height = 6,
  units = "cm"
)
saveRDS(
  object = barplot_mean_no_legend,
  file = "1_rds_files/Hordeum_bulk_order_top20_RA_mean_stackedbp_no_legend.rds"
)
