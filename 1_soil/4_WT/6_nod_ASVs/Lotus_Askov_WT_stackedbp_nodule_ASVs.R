# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directoryto source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

# Load the input files.
design <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", header=T, sep="\t"
)
asv_table <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = ""
)

taxonomy <- read.table("../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv", sep="\t", header=TRUE, fill=TRUE)

# Clean up the taxonomy layout.
colnames(taxonomy)[colnames(taxonomy) == "Feature.ID"] <- "ASVid"
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Filter ASV table to samples in design and ASVs in taxonomy file.
samples_to_keep <- intersect(colnames(asv_table), design$SampleID)
asv_table <- asv_table[rownames(asv_table) %in% taxonomy$ASVid, samples_to_keep]

# Convert the ASV counts to relative abundances in the ASV table.
asv_table <- sweep(asv_table, 2, colSums(asv_table), "/")

# Filter for the nodule samples and add the genus info from the taxonomy file.
nodule_samples <- design %>% filter(Compartment == "Nodules") %>% pull(SampleID)

Lotus_nod <- as.data.frame(asv_table) %>%
  rownames_to_column("ASVid") %>%
  select(ASVid, all_of(nodule_samples)) %>%
  left_join(taxonomy %>% select(ASVid, Genus), by = "ASVid")

# Calculate the average relative abundance (RA) per ASV per soil.
sample_soil_map <- design %>%
  filter(SampleID %in% nodule_samples) %>%
  select(SampleID, Soil)

avg_per_asv <- Lotus_nod %>%
  pivot_longer(cols = -c(ASVid, Genus),
               names_to = "SampleID",
               values_to = "RA") %>%
  left_join(sample_soil_map, by = "SampleID") %>%
  group_by(ASVid, Genus, Soil) %>%
  summarise(Average = mean(RA, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Soil, values_from = Average, values_fill = 0)

# Filter for the top20 ASvs per soil with >=1% RA.
top20_npk <- avg_per_asv %>% filter(NPK >= 0.01) %>% arrange(desc(NPK)) %>% slice(1:20) %>% pull(ASVid)
top20_pk  <- avg_per_asv %>% filter(PK >= 0.01)  %>% arrange(desc(PK))  %>% slice(1:20) %>% pull(ASVid)
top20_uf  <- avg_per_asv %>% filter(UF >= 0.01)  %>% arrange(desc(UF))  %>% slice(1:20) %>% pull(ASVid)

topASVs <- unique(c(top20_npk, top20_pk, top20_uf))

# Keep only the top ASVs for plotting.
avg_topASVs <- avg_per_asv %>% filter(ASVid %in% topASVs)

# Assign the genus and numeric suffix alphabetically.
avg_topASVs <- avg_topASVs %>%
  arrange(Genus, ASVid) %>%
  group_by(Genus) %>%
  mutate(Isolate = paste0(Genus, "_", row_number())) %>%
  ungroup()

# Prepare a long format for plotting.
df_long <- avg_topASVs %>%
  pivot_longer(cols = c("NPK", "PK", "UF"),
               names_to = "Soil_type",
               values_to = "RA") %>%
  mutate(
    Soil_type = factor(Soil_type, levels = c("NPK", "PK", "UF")),
    Isolate = factor(Isolate, levels = sort(unique(Isolate)))
  )

# Two isolates have genus "Unknown", therefore we will check the order and family info for those and change the names for the plot.
unknown_asvs <- df_long %>% filter(Genus == "Unknown")
unknown_asv_ids <- unique(unknown_asvs$ASVid)

unknown_tax_info <- taxonomy %>%
  filter(ASVid %in% unknown_asv_ids) %>%
  select(ASVid, Order, Family, Genus)

unknown_tax_info

levels(df_long$Isolate) <- c(levels(df_long$Isolate), "Rhizobiaceae_Unknown_1")
df_long$Isolate[df_long$Isolate == "Unknown_2"] <- "Rhizobiaceae_Unknown_1"

# Define the new factor levels.
custom_order <- c(
  "Mesorhizobium_1", "Mesorhizobium_2", "Mesorhizobium_3",
  "Mesorhizobium_4", "Mesorhizobium_5", "Mesorhizobium_6",
  "Rhizobiaceae_Unknown_1", 
  "Pseudomonas_1",
  "Unknown_1"
)

df_long$Isolate <- factor(df_long$Isolate, levels = custom_order)

# Define a custom colour palette for the plot.
custom_colors <- data.frame(
  group = c("Pseudomonas_1", "Mesorhizobium_1", "Mesorhizobium_2", "Mesorhizobium_3",
            "Mesorhizobium_4", "Mesorhizobium_5", "Mesorhizobium_6", "Rhizobiaceae_Unknown_1", "Unknown_1"),
  colors = c("#88CCAA","#fff9e7","#ffeeb6","#ffe385","#ffdd6c","#FFCC00","#FFB300","#F58700","darkgrey")
)

# Map the colours to the isolates.
color_map <- setNames(custom_colors$colors, custom_colors$group)
fill_colors <- color_map[as.character(df_long$Isolate)]

# Set the main theme for the plot.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size=8, color="black"),
  legend.text=element_text(size=8, color="black"),
  legend.key=element_blank(),
  axis.title.y=element_text(size=8),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size=8, color="black")
)

# Make a stacked barplot.
p <- ggplot(df_long, aes(x = Soil_type, y = RA, fill = Isolate)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0))+
  labs(x = "", y = "Relative Abundance", title="Nodule ASVs (top20 ≥1% RA)") +
  main_theme +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.25, 'cm'),
    legend.margin = margin(l = -8)
  )

p

# Save the plot.
ggsave("Lotus_Askov_WT_stackedbp_NoduleASVs.pdf",p, width = 10, height = 6, units = "cm")
saveRDS(p, file = "Lotus_Askov_WT_stackedbp_NoduleASVs.rds")
saveRDS(p, file = "../7_final_figures/Lotus_Askov_WT_stackedbp_NoduleASVs.rds")
