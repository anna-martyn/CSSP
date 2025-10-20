# Script written by Anna Martyn, Oct25. amartyn@mpipz.mpg.de

# -----------------------------
# Clean up
# -----------------------------
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

# -----------------------------
# Load files
# -----------------------------
design <- read.table("Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", header=TRUE, sep="\t")
asv_table <- read.table("feature-table.tsv", sep="\t", header=TRUE, row.names=1,
                        check.names=FALSE, comment.char = "", skip = 1)
taxonomy <- read.table("taxonomy.tsv", sep="\t", header=TRUE, fill=TRUE)

# Remove spike-in ASV
asv_table <- asv_table[row.names(asv_table) != "85fa8bb918a926d97659d9b64ca6fedd", ]

# -----------------------------
# Clean-up taxonomy
# -----------------------------
colnames(taxonomy)[colnames(taxonomy) == "Feature.ID"] <- "ASVid"

taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# -----------------------------
# Filter ASV table to samples in design and ASVs in taxonomy
# -----------------------------
samples_to_keep <- intersect(colnames(asv_table), design$SampleID)
asv_table <- asv_table[rownames(asv_table) %in% taxonomy$ASVid, samples_to_keep]

# Convert to relative abundance
asv_table <- sweep(asv_table, 2, colSums(asv_table), "/")

# -----------------------------
# Keep only nodule samples and add Genus
# -----------------------------
nodule_samples <- design %>% filter(Compartment == "Nodules") %>% pull(SampleID)

Lotus_nod <- as.data.frame(asv_table) %>%
  rownames_to_column("ASVid") %>%
  select(ASVid, all_of(nodule_samples)) %>%
  left_join(taxonomy %>% select(ASVid, Genus), by = "ASVid")

# -----------------------------
# Compute average RA per ASV per soil
# -----------------------------
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

# -----------------------------
# Filter top 20 ASVs per soil with >=1% RA
# -----------------------------
top20_npk <- avg_per_asv %>% filter(NPK >= 0.01) %>% arrange(desc(NPK)) %>% slice(1:20) %>% pull(ASVid)
top20_pk  <- avg_per_asv %>% filter(PK >= 0.01)  %>% arrange(desc(PK))  %>% slice(1:20) %>% pull(ASVid)
top20_uf  <- avg_per_asv %>% filter(UF >= 0.01)  %>% arrange(desc(UF))  %>% slice(1:20) %>% pull(ASVid)

topASVs <- unique(c(top20_npk, top20_pk, top20_uf))

# Keep only topASVs for plotting
avg_topASVs <- avg_per_asv %>% filter(ASVid %in% topASVs)

# -----------------------------
# Assign Genus + numeric suffix alphabetically
# -----------------------------
avg_topASVs <- avg_topASVs %>%
  arrange(Genus, ASVid) %>%
  group_by(Genus) %>%
  mutate(Isolate = paste0(Genus, "_", row_number())) %>%
  ungroup()

# -----------------------------
# Prepare long-format for plotting
# -----------------------------
df_long <- avg_topASVs %>%
  pivot_longer(cols = c("NPK", "PK", "UF"),
               names_to = "Soil_type",
               values_to = "RA") %>%
  mutate(
    Soil_type = factor(Soil_type, levels = c("NPK", "PK", "UF")),
    Isolate = factor(Isolate, levels = sort(unique(Isolate)))
  )

# -----------------------------
# Define custom color palette
# -----------------------------
custom_colors <- data.frame(
  group = c("Pseudomonas_1", "Mesorhizobium_1", "Mesorhizobium_2", "Mesorhizobium_3",
            "Mesorhizobium_4", "Mesorhizobium_5", "Mesorhizobium_6", "Unknown_1", "Unknown_2"),
  colors = c("#88CCAA","#fff9e7","#ffeeb6","#ffe385","#ffdd6c","#FFCC00","#FFB300","lightgrey","darkgrey")
)

# Map colors to Isolate
color_map <- setNames(custom_colors$colors, custom_colors$group)
fill_colors <- color_map[as.character(df_long$Isolate)]  # only Isolates present

# -----------------------------
# Define main theme
# -----------------------------
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size=20, color="black"),
  legend.text=element_text(size=20, color="black"),
  legend.key=element_blank(),
  axis.title.y=element_text(size=20),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size=20, color="black")
)

# -----------------------------
# Plot stacked barplot
# -----------------------------
p <- ggplot(df_long, aes(x = Soil_type, y = RA, fill = Isolate)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0))+
  labs(x = "", y = "Relative Abundance", title="Nodule ASVs (top20 ≥1% RA)") +
  main_theme +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0),
    legend.title = element_text(size = 20)
  )


p

# Save plot.
ggsave("Soil_LotusWT_stackedbp_NoduleASVs.pdf", p, width = 10, height = 6)
saveRDS(p, file = "Soil_LotusWT_stackedbp_NoduleASVs.rds")
