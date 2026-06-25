# Seup ------------------------------------------------------------------------
pkg <- c("data.table", "ggplot2", "cowplot", "ggtext")
for(pk in pkg){
  library(pk, character.only = T)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Genotype colours
cols <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD",
  "control" = "#cecece"
)

# Loading Lotus data ----------------------------------------------------------
# Lotus feature table filtered for background features
metabolite_data_Lj <- fread(
  "../2_background_removal/1_tables/feature_table_Lotus_filtered.csv"
)

design_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_RootEx_metadata.txt",
  drop = 4:7
)

# Setting sample names in metadata
design_Lj[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Lj <- design_Lj[ Sample_ID %in% colnames(metabolite_data_Lj)]

# Removing control samples from metadata and feature table
non_control_samples <- design_Lj[Genotype != "control", Sample_ID]
design_Lj <- design_Lj[Sample_ID %in% non_control_samples]
metabolite_data_Lj <- metabolite_data_Lj[,
  c("Feature", non_control_samples),
  with = FALSE
]

# Loading Hordeum data --------------------------------------------------------
metabolite_data_Hv <- fread(
  "../2_background_removal/1_tables/feature_table_Hordeum_filtered.csv"
)
design_Hv <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_RootEx_metadata.txt",
  drop = c(2, 4:7)
)

# Setting sample names in metadata
design_Hv[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Hv <- design_Hv[ Sample_ID %in% colnames(metabolite_data_Hv)]

# Principal Component Analysis ------------------------------------------------
# Lagend labels
legend_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*",
  "control" = "control"
)

# PCA
pca_Lj <- prcomp(t(metabolite_data_Lj[, -1]), center = TRUE, scale. = TRUE)
pca_Hv <- prcomp(t(metabolite_data_Hv[, -1]), center = TRUE, scale. = TRUE)

# Variance explained for each principal component
## Lotus
var_exp_Lj <- (pca_Lj$sdev^2/sum(pca_Lj$sdev^2))[1:2]*100
var_exp_Lj <- paste0("(", round(var_exp_Lj, 2), "%)")

## Hordeum
var_exp_Hv <- (pca_Hv$sdev^2/sum(pca_Hv$sdev^2))[1:2]*100
var_exp_Hv <- paste0("(", round(var_exp_Hv, 2), "%)")

# Extracting points from PCA
## Lotus
pca_dt_Lj <- data.table(
  pca_Lj$x[, 1:2],
  Genotype = design_Lj$Genotype,
  Host = "Lotus"
)

## Hordeum
pca_dt_Hv <- data.table(
  pca_Hv$x[, 1:2],
  Genotype = design_Hv$Genotype,
  Host = "Hordeum"
)

# Centroids and segments
centroids_Lj <- pca_dt_Lj[,
  .(PC1_cent = mean(PC1), PC2_cent = mean(PC2)),
  list(Genotype, Host)
]
segments_Lj <- merge(pca_dt_Lj, centroids_Lj, by = c("Genotype", "Host"))

# Setting factor levels
pca_dt_Lj[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

pca_plot_Lj <- ggplot(data = pca_dt_Lj, mapping = aes(x = PC1, y = PC2, colour = Genotype)) +
  geom_point(size = 1.5, stroke = 0.25) +
  geom_segment(
    data = segments_Lj,
    aes(x = PC1, y = PC2, xend = PC1_cent, yend = PC2_cent, color = Genotype),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = cols, labels = legend_labels) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(colour = 'black', size = 6, face = "bold"),
    strip.background = element_rect(colour = NA),
    legend.position = "bottom",
    axis.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    legend.text = element_markdown(
      size = 6,
      margin = margin(l = -0.1, unit = "pt")
    ),
    legend.title = element_text(size = 6),
    legend.margin = margin(t = 0, r = 5, l = 5),
    legend.key = element_rect(fill = NA),
    legend.key.spacing.y = unit(-0.1, "cm"),
    legend.key.spacing.x = unit(0.001, "cm"),
    plot.title = element_text(size = 6, hjust = 0.5)
  ) +
  labs(x = paste("PC1", var_exp_Lj[1]), y = paste("PC2", var_exp_Lj[2])) +
  facet_wrap(~Host) +
  guides(fill = guide_legend(nrow = 1)) +
  NULL

centroids_Hv <- pca_dt_Hv[,
  .(PC1_cent = mean(PC1), PC2_cent = mean(PC2)),
  list(Genotype, Host)
]

segments_Hv <- merge(pca_dt_Hv, centroids_Hv, by = c("Genotype", "Host"))

pca_plot_Hv <- ggplot(
  data = pca_dt_Hv,
  mapping = aes(x = PC1, y = PC2, colour = Genotype)
) +
  geom_point(size = 1.5, stroke = 0.25) +
  scale_colour_manual(name = "Genotype", breaks = names(cols), values = cols) +
  geom_segment(
    data = segments_Hv,
    aes(x = PC1, y = PC2, xend = PC1_cent, yend = PC2_cent, color = Genotype),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(colour = 'black', size = 6, face = "bold"),
    strip.background = element_rect(colour = NA),
    legend.position = "none",
    axis.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    legend.text = element_text(
      size = 6,
      margin = margin(l = -0.1, unit = "pt")
    ),
    legend.title = element_text(size = 6),
    legend.margin = margin(t = 0, r = 5, l = 5),
    legend.key = element_rect(fill = NA),
    legend.key.spacing.y = unit(-0.1, "cm"),
    legend.key.spacing.x = unit(0.001, "cm"),
    plot.title = element_text(size = 6, hjust = 0.5)
  ) +
  labs(x = paste("PC1", var_exp_Hv[1]), y = paste("PC2", var_exp_Hv[2])) +
  facet_wrap(~Host) +
  guides(fill = guide_legend(nrow = 1)) +
  NULL

pca_legend <- get_plot_component(pca_plot_Lj, 'guide-box', return_all = TRUE)
pca_plot_Lj <- pca_plot_Lj + guides(colour = "none")

pca_all <- plot_grid(pca_plot_Lj, pca_plot_Hv, nrow = 2)
pca_all + guides(colour = NULL)

# Saving plot
saveRDS(object = pca_all, file = "1_rds_files/PCA_plot.rds")
saveRDS(object = pca_legend, file = "1_rds_files/PCA_plot_legend.rds")

ggsave(
  filename = "2_figures/PCA_plot.pdf",
  plot = pca_all,
  height = 7,
  width = 7,
  units = "cm"
)
