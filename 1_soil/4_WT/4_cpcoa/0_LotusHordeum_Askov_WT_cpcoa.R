# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
library("ggplot2")
library("vegan")

# Loading file with cpcoa functions
source("cpcoa.func.R")

# Loading Lotus data
Lotus_design <- read.table(
  file = "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt",
  header = TRUE,
  sep = "\t"
)
Lotus_asv_table <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1, 
  check.names = FALSE,
  comment.char = ""
)

# Loading Hordeum data
Hordeum_design <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
  header = TRUE,
  sep = "\t"
)
Hordeum_asv_table <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1, 
  check.names = FALSE,
  comment.char = "",
  skip = 1
)

# Keeping only samples present in both design and ASV table
idx <- intersect(Lotus_design$SampleID, colnames(Lotus_asv_table))
Lotus_design <- Lotus_design[match(idx, Lotus_design$SampleID), ]
Lotus_asv_table <- Lotus_asv_table[, idx]

idx <- intersect(Hordeum_design$SampleID, colnames(Hordeum_asv_table))
Hordeum_design <- Hordeum_design[match(idx, Hordeum_design$SampleID), ]
Hordeum_asv_table <- Hordeum_asv_table[, idx]

# Keeping only WT samples
Lotus_wt_samples <- Lotus_design$SampleID[Lotus_design$Genotype == "WT"]
Lotus_design <- Lotus_design[Lotus_design$SampleID %in% Lotus_wt_samples, ]
Lotus_asv_table <- Lotus_asv_table[, Lotus_wt_samples]

Hordeum_wt_samples <- Hordeum_design$SampleID[Hordeum_design$Genotype == "WT"]
Hordeum_design <- Hordeum_design[Hordeum_design$SampleID %in% Hordeum_wt_samples, ]
Hordeum_asv_table <- Hordeum_asv_table[, Hordeum_wt_samples]

# Bray-Curtis distances
Lotus_asv_table_norm <- apply(Lotus_asv_table, 2, function(x) x / sum(x))
bray_curtis_Lotus <- vegdist(t(Lotus_asv_table_norm), method = "bray")

Hordeum_asv_table_norm <- apply(Hordeum_asv_table, 2, function(x) x / sum(x))
bray_curtis_Hordeum <- vegdist(t(Hordeum_asv_table_norm), method = "bray")

# Colours and shapes
colors <- data.frame(
  group = c("NPK", "PK", "UF"),
  colors = c("#6F944F", "#B2563C", "#3C7D82")
)

shapes <- data.frame(group = c("Rhizosphere", "Root", "Nodules"), shapes = 15:17)

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  text = element_text(size = 6, color = "black"),
  legend.position = "none",
  legend.background = element_blank(),
  plot.title = element_text(size = 6, hjust = 0.9)
)

# Lotus CPCoA and plot --------------------------------------------------------
# CPCoA
CPCoA_Lj <- capscale(
  t(Lotus_asv_table_norm) ~ Compartment*Soil,
  data = Lotus_design,
  add = FALSE, 
  sqrt.dist = TRUE, 
  distance = "bray"
)

# PERMANOVA
set.seed(1773665952)
perm_anova <- anova.cca(CPCoA_Lj)
print(perm_anova)

# Variability tables and confidence intervals for the variance
var_tbl <- variability_table(CPCoA_Lj)
eig <- CPCoA_Lj$CCA$eig
variance <- var_tbl["constrained", "proportion"]
p_val <- perm_anova[1, 4]

## Extracting the weighted average (sample) scores
points <- CPCoA_Lj$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points <- cbind(points, Lotus_design[match(rownames(points), Lotus_design$SampleID), ])

points$Soil <- factor(points$Soil, levels = colors$group)
points$Compartment <- factor(points$Compartment, levels = shapes$group)

# Centroids by soil-compartment combination
centroids <- aggregate(cbind(x,y) ~ Soil + Compartment, data = points, FUN = mean)
segments <- merge(
  x = points, 
  y = setNames(centroids, c("Soil","Compartment","seg_x","seg_y")),
  by = c("Soil", "Compartment"),
  sort = FALSE
)

# Plot
CPCoA_plot_Lj <- ggplot(
  points,
  aes(x = x, y = y, color = Soil, shape = Compartment)
) +
  geom_segment(
    data = segments,
    mapping = aes(x = x, y = y, xend = seg_x, yend = seg_y, color = Soil),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_point(size = 3, alpha = 0.7) +
  scale_colour_manual(values = as.character(colors$color)) +
  scale_shape_manual(values = shapes$shape) +
  labs(
    x = paste(
      "CPCo 1 (",
      format(100 * eig[1] / sum(eig), digits = 4),
      "%)",
      sep = ""
    ),
    y = paste(
      "CPCo 2 (",
      format(100 * eig[2] / sum(eig), digits = 4),
      "%)",
      sep = ""
    )
  ) +
  ggtitle(
    "Lotus",
    subtitle = paste(
      format(100 * variance, digits = 3),
      "% of variance; p<",
      format(p_val, digits = 2),
      sep = ""
    )
  ) +
  main_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    plot.title = element_text(face = "bold", size = 6, hjust = 0),
    plot.subtitle = element_text(size = 6, hjust = 0)
  )

# Saving plot
ggsave(
  filename = "2_figures/Lotus_Askov_WT_cpcoa.pdf",
  CPCoA_plot_Lj,
  width = 5,
  height = 5,
  units = "cm"
)

saveRDS(object = CPCoA_plot_Lj, file = "1_rds_files/Lotus_Askov_WT_cpcoa.rds")

# Hordeum CPCoA and plot ------------------------------------------------------
CPCoA_Hv <- capscale(
  t(Hordeum_asv_table_norm) ~ Compartment*Soil, 
  data = Hordeum_design, 
  add = FALSE, 
  sqrt.dist = TRUE, 
  distance = "bray"
)

# Permanova
perm_anova <- anova.cca(CPCoA_Hv)
print(perm_anova)

# Variability tables and confidence intervals for the variance
var_tbl <- variability_table(CPCoA_Hv)
eig <- CPCoA_Hv$CCA$eig
variance <- var_tbl["constrained", "proportion"]
p_val <- perm_anova[1, 4]

# Extracting the weighted average (sample) scores
points <- CPCoA_Hv$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, Hordeum_design[match(rownames(points), Hordeum_design$SampleID), ])

points$Soil <- factor(points$Soil, levels = colors$group)
points$Compartment <- factor(points$Compartment, levels = shapes$group)

# Identifying centroids for each soil-compartment combination
centroids <- aggregate(cbind(x,y) ~ Soil + Compartment, data = points, FUN = mean)
segments <- merge(
  x = points,
  y = setNames(centroids, c("Soil","Compartment","seg_x","seg_y")),
  by = c("Soil","Compartment"),
  sort = FALSE
)

# Plot
CPCoA_plot_Hv <- ggplot(
  points,
  aes(x = x, y = y, color = Soil, shape = Compartment)
) +
  geom_segment(
    data = segments,
    mapping = aes(x = x, y = y, xend = seg_x, yend = seg_y, color = Soil),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_point(size = 3, alpha = 0.7) +
  scale_colour_manual(values = as.character(colors$color)) +
  scale_shape_manual(values = shapes$shape) +
  labs(
    x = paste(
      "CPCo 1 (",
      format(100 * eig[1] / sum(eig), digits = 4),
      "%)",
      sep = ""
    ),
    y = paste(
      "CPCo 2 (",
      format(100 * eig[2] / sum(eig), digits = 4),
      "%)",
      sep = ""
    )
  ) +
  ggtitle(
    "Hordeum",
    subtitle = paste(
      format(100 * variance, digits = 3),
      "% of variance; p<",
      format(p_val, digits = 2),
      sep = ""
    )
  ) +
  main_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    plot.title = element_text(face = "bold", size = 6, hjust = 0),
    plot.subtitle = element_text(size = 6, hjust = 0)
  ) 

# Saving plot
ggsave(
  filename = "2_figures/Hordeum_Askov_WT_cpcoa.pdf",
  plot = CPCoA_plot_Hv,
  width = 5,
  height = 5
)
saveRDS(object = CPCoA_plot_Hv, file = "1_rds_files/Hordeum_Askov_WT_cpcoa.rds")
