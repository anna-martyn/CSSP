# Seup ------------------------------------------------------------------------
pkg <- c("data.table", "ggplot2", "cowplot")
for(pk in pkg){
  library(pk, character.only = TRUE)
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

legend_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*",
  "control" = "control"
)

## Loading Lotus data ---------------------------------------------------------
# Lotus feature table filtered for background features
metabolite_data_Lj <- fread(
  "../2_background_removal/1_tables/feature_table_Lotus_filtered.csv"
)

design_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_rootex_metadata.txt",
  drop = 4:6
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

## Loading Hordeum data -------------------------------------------------------
metabolite_data_Hv <- fread(
  "../2_background_removal/1_tables/feature_table_Hordeum_filtered.csv"
)
design_Hv <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_rootex_metadata.txt",
  drop = c(2, 4:7)
)

# Setting sample names in metadata
design_Hv[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Hv <- design_Hv[ Sample_ID %in% colnames(metabolite_data_Hv)]

# Adding plant variable to hordeum metadata
design_Hv[,Plant := "Hordeum"]

# Volcano plots ----------------------------------------------------------------
# Loading p-values and log-fold changes from genotype effects
p_adj_Lj <- fread("../4_genotype_effects/1_tables/p_adj_Lj.csv")
p_adj_Hv <- fread("../4_genotype_effects/1_tables/p_adj_Hv.csv")
lfc_Lj <- fread("../4_genotype_effects/1_tables/lfc_Lj.csv")
lfc_Hv <- fread("../4_genotype_effects/1_tables/lfc_Hv.csv")

# Adding plant-variable
p_adj_Lj[, Plant := "Lotus"]
lfc_Lj[, Plant := "Lotus"]
p_adj_Hv[, Plant := "Hordeum"]
lfc_Hv[, Plant := "Hordeum"]

# Combining Lotus and Hordeum results
p_vals_dt <- rbind(p_adj_Lj, p_adj_Hv)
lfc_dt <- rbind(lfc_Lj, lfc_Hv)

# Long form
p_vals_dt <- melt(
  p_vals_dt,
  id.vars = c("Feature", "Plant"),
  variable.name = "Genotype",
  value.name = "p_adj"
)
lfc_dt <- melt(
  lfc_dt,
  id.vars = c("Feature", "Plant"),
  variable.name = "Genotype",
  value.name = "logFC"
)

# Cleaning up Genotype names
p_vals_dt[, Genotype := gsub("_p_adj", "", Genotype)]
lfc_dt[, Genotype := gsub("Lfc_", "", Genotype)]

# Combining log-fold changes and p-values
res_table <- merge(lfc_dt, p_vals_dt, by = c("Feature", "Genotype", "Plant"))

# Converting log-fold changes to log2-fold changes
res_table[, logFC := logFC / log(2)]

# Discriminating between depleted, enriched, and NS (non-significant) features
res_table[, Sig := p_adj < 0.05]
res_table[,
  diff := fcase(
    Sig == TRUE & logFC > 0 , "Enriched" ,
    Sig == TRUE & logFC < 0 , "Depleted" ,
    default = "NS"
  )
]

# Setting factor levels
res_table[, ":="(
  Genotype = factor(Genotype, levels = c("symrk", "ccamk", "nsp1", "nsp2")),
  Plant = factor(Plant, levels = c("Lotus", "Hordeum"))
)]

# Tallying the number of significant features by plant-genotype combination
text_data <- res_table[, .(N_diff = sum(Sig)), list(Plant, Genotype)]

# Plot
volcano_plot <- ggplot(
  data = res_table,
  aes(x = logFC, y = -log10(p_adj), colour = diff)
) +
  geom_point(size = 0.5) +
  xlab(expression("log"[2] * "FC vs WT")) +
  ylab(expression("-log"[10] * "p-value (adjusted for FDR)")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_light() +
  scale_color_manual(
    values = c("Enriched" = "#902121", "Depleted" = "darkblue", "NS" = "gray"),
    name = "DEM"
  ) +
  facet_grid(factor(Plant, levels = c("Lotus", "Hordeum")) ~ Genotype) +
  geom_label(
    data = text_data,
    aes(x = -25, y = 8.5, label = N_diff),
    colour = "black",
    fill = "grey",
    alpha = 0.2,
    size = 6 / .pt
  ) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    plot.title = element_text(size = 6, hjust = 0.5),
    strip.text = element_text(colour = 'black', size = 6, face = "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightgrey"),
    legend.key = element_rect(fill = NA),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.size = unit(5, "pt"),
    legend.box.spacing = unit(5, "pt"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6, hjust = 0.5, margin = margin(r = 5))
  ) +
  xlim(-30, 30) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  NULL

volcano_legend <- get_plot_component(volcano_plot, 'guide-box', return_all = TRUE)
volcano_plot <- volcano_plot + guides(colour = "none")

saveRDS(object = volcano_legend, file = "1_rds_files/volcano_legend.rds")
saveRDS(object = volcano_plot, file = "1_rds_files/volcano_plot.rds")

ggsave(
  filename = "2_figures/volcano_plot.pdf",
  plot = volcano_plot,
  width = 14,
  height = 7,
  units = "cm"
)
