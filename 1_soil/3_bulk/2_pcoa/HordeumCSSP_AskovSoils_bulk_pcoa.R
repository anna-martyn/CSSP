# Seup ------------------------------------------------------------------------
# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
pkg <- c("ggplot2", "dplyr", "tidyr", "vegan")
for(pk in pkg){
  library(pk, character.only = T)
}

# Load data
design <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
  sep = "\t",       
  header = TRUE,    
  check.names = FALSE
)
taxonomy <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)
asv_table <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  skip = 1,
  comment.char = ""
)

# Clean up the taxonomy file
taxonomy <- taxonomy %>% rename(ASVid = Feature.ID)
taxonomy <- taxonomy %>%
  separate(
    col = Taxon,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = "; ",
    fill = "right"
  ) %>%
  mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Keep soil samples
design_bulk <- design %>% filter(Genotype == "Soil")
asv_table_bulk <- asv_table[, design_bulk$SampleID]

# Define colours and shapes for the plot
colors <- data.frame(
  group = c("NPK", "PK", "UF"),
  colors = c("#6F944F", "#B2563C", "#3C7D82")
)

shapes <- data.frame(
  group = c("NPK", "PK", "UF"),
  shape = c(19, 19, 19)
)

colors <- colors[colors$group %in% design_bulk$Soil,]
shapes <- shapes[shapes$group %in% design_bulk$Soil,]

# Principal coordinate analysis -----------------------------------------------

# Convert ASV reads to relative abundances
asv_table_norm <- apply(asv_table_bulk, 2, function(x) x / sum(x))

# Bray-Curtis distances
bray_curtis <- vegdist(t(asv_table_norm), method = "bray")

# Principal coordinate analysis with Bray-Curtis distances
k <- 2
pcoa <- cmdscale(bray_curtis, k = k, eig = TRUE)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_bulk[match(rownames(points), design_bulk$SampleID), ])

# Figure ----------------------------------------------------------------------
colors <- colors[colors$group %in% points$Soil, ]
points$Soil <- factor(points$Soil, levels = colors$group)

shapes <- shapes[shapes$group %in% points$Soil, ]
points$Soil <- factor(points$Soil, levels = shapes$group)

# Calculate centroids by soil
centroids <- aggregate(cbind(points$x, points$y) ~ Soil, data = points, FUN = mean)

# Join the centroids back to the points so each sample has its group centroid
segments <- merge(
  x = points, 
  y = setNames(centroids, c('Soil','seg_x','seg_y')),
  by = "Soil",
  sort = FALSE
)

# Set main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color="black"),
  axis.line.y = element_line(color="black"),
  axis.ticks = element_line(color="black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.title = element_text(size = 6, color = "black", hjust = 0.5, margin = margin(b = -2)),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  text = element_text(size = 6, color="black"),
  legend.position = "right",
  legend.margin = margin(l = -10),
  legend.background = element_blank(),
  plot.title = element_text(size = 6, hjust=0.9),
  legend.key.spacing.y = unit(-0.3, 'cm')
)

# Plot
p <- ggplot(points, aes(x = x, y = y, colour = Soil)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_segment(
    data = segments,
    mapping = aes(x = x, y = y, xend = seg_x, yend = seg_y, colour = Soil,),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values=as.character(colors$colors)) +
  labs(
    x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
    y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = "")
  ) +
  main_theme+
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.2))+
  NULL

p

# Save
ggsave("Hordeum_bulk_PCoA.pdf", p, width = 5, height = 5, units = "cm")
saveRDS(p, file = "Hordeum_bulk_PCoA.rds")
saveRDS(p, file = "../5_final_figure/Hordeum_bulk_PCoA.rds")

