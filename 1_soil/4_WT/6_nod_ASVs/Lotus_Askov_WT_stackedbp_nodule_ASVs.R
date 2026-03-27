# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directoryto source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading the required packages
pkg <- c("ggplot2", "dplyr", "tidyr", "tibble")
for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
design <- read.table(
  file = "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt",
  header = TRUE,
  sep = "\t"
)
asv_table <- read.table(
  file = "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1, 
  check.names = FALSE,
  comment.char = ""
)

taxonomy <- read.table(
  file = "../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

# Cleaning up taxonomy
colnames(taxonomy)[colnames(taxonomy) == "Feature.ID"] <- "ASVid"
taxonomy <- taxonomy %>%
  separate(
    col = Taxon,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = "; ", fill = "right"
  ) %>%
  mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Filtering ASV table to keep only samples in design and ASVs in taxonomy
samples_to_keep <- intersect(colnames(asv_table), design$SampleID)
asv_table <- asv_table[rownames(asv_table) %in% taxonomy$ASVid, samples_to_keep]

# ASV table with relative abundances (RA)
asv_table <- sweep(asv_table, 2, colSums(asv_table), "/")

# Keeping only nodule samples and adding Genus information to ASVs
nodule_samples <- design %>% filter(Compartment == "Nodules") %>% pull(SampleID)

Lotus_nod <- as.data.frame(asv_table) %>%
  rownames_to_column("ASVid") %>%
  select(ASVid, all_of(nodule_samples)) %>%
  left_join(taxonomy %>% select(ASVid, Genus), by = "ASVid")

# Mean RA per ASV by soil
sample_soil_map <- design %>%
  filter(SampleID %in% nodule_samples) %>%
  select(SampleID, Soil)

avg_per_asv <- Lotus_nod %>%
  pivot_longer(
    cols = -c(ASVid, Genus), names_to = "SampleID", values_to = "RA"
  ) %>%
  left_join(sample_soil_map, by = "SampleID") %>%
  group_by(ASVid, Genus, Soil) %>%
  summarise(Average = mean(RA, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Soil, values_from = Average, values_fill = 0)

# Identifying top 20 ASVs with >=1% RA for each soil
top20_npk <- avg_per_asv %>% 
  filter(NPK >= 0.01) %>% 
  arrange(desc(NPK)) %>% 
  slice(1:20) %>% 
  pull(ASVid)

top20_pk  <- avg_per_asv %>%
  filter(PK >= 0.01) %>%
  arrange(desc(PK)) %>%
  slice(1:20) %>%
  pull(ASVid)

top20_uf  <- avg_per_asv %>%
  filter(UF >= 0.01) %>%
  arrange(desc(UF)) %>%
  slice(1:20) %>%
  pull(ASVid)

top_asv <- unique(c(top20_npk, top20_pk, top20_uf))

# Keeping only top ASVs for plotting
avg_top_asv <- avg_per_asv %>% filter(ASVid %in% top_asv)

# Assigning the genus and numeric suffix alphabetically.
avg_top_asv <- avg_top_asv %>%
  arrange(Genus, ASVid) %>%
  group_by(Genus) %>%
  mutate(Isolate = paste0(Genus, "_", row_number())) %>%
  ungroup()

# Long form for plotting
df_long <- avg_top_asv %>%
  pivot_longer(
    cols = c("NPK", "PK", "UF"),
    names_to = "Soil_type",
    values_to = "RA"
  ) %>%
  mutate(
    Soil_type = factor(Soil_type, levels = c("NPK", "PK", "UF")),
    Isolate = factor(Isolate, levels = sort(unique(Isolate)))
  )

# Adding order and family information to isolates of unknown genus
unknown_asvs <- df_long %>% filter(Genus == "Unknown")
unknown_asv_ids <- unique(unknown_asvs$ASVid)

unknown_tax_info <- taxonomy %>%
  filter(ASVid %in% unknown_asv_ids) %>%
  select(ASVid, Order, Family, Genus)

levels(df_long$Isolate) <- c(levels(df_long$Isolate), "Rhizobiaceae_Unknown_1")
df_long$Isolate[df_long$Isolate == "Unknown_2"] <- "Rhizobiaceae_Unknown_1"

# Setting factor levels
custom_order <- c(
  "Mesorhizobium_1", "Mesorhizobium_2", "Mesorhizobium_3",
  "Mesorhizobium_4", "Mesorhizobium_5", "Mesorhizobium_6",
  "Rhizobiaceae_Unknown_1", 
  "Pseudomonas_1",
  "Unknown_1"
)

df_long$Isolate <- factor(df_long$Isolate, levels = custom_order)

# Barplot ---------------------------------------------------------------------
# Defining custom colour palette
custom_colors <- data.frame(
  group = c(
    "Pseudomonas_1", "Mesorhizobium_1", "Mesorhizobium_2", "Mesorhizobium_3",
    "Mesorhizobium_4", "Mesorhizobium_5", "Mesorhizobium_6", 
    "Rhizobiaceae_Unknown_1", "Unknown_1"
  ),
  colors = c(
    "#88CCAA", "#fff9e7", "#ffeeb6", "#ffe385", "#ffdd6c",
    "#FFCC00", "#FFB300", "#F58700", "darkgrey"
  )
)

# Mapping colours to isolates
color_map <- setNames(custom_colors$colors, custom_colors$group)
fill_colors <- color_map[as.character(df_long$Isolate)]

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 8, color = "black"),
  legend.text = element_text(size = 8, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 8),
  legend.position = "right",
  legend.background = element_blank(),
  text = element_text(family = "sans", size = 8, color = "black")
)

# Stacked barplot
bar_plot <- ggplot(df_long, aes(x = Soil_type, y = RA, fill = Isolate)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0))+
  labs(x = "", y = "Relative Abundance", title = "Nodule ASVs (top20 ≥1% RA)") +
  main_theme +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.25, 'cm'),
    legend.margin = margin(l = -8)
  )

# Saving plot
ggsave(
  filename = "2_figures/Lotus_Askov_WT_stackedbp_NoduleASVs.pdf",
  plot = bar_plot,
  width = 10,
  height = 6,
  units = "cm"
)
saveRDS(
  object = bar_plot,
  file = "1_rds_files/Lotus_Askov_WT_stackedbp_NoduleASVs.rds"
)
