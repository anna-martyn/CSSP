# Seup ------------------------------------------------------------------------
pkg <- c("data.table", "ggplot2", "ggh4x", "cowplot")
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

# Boxplots for highlighted features -------------------------------------------
## Lotus ----------------------------------------------------------------------
highlighted_features_Lj <- paste0(
  "Feature",
  c(269, 455, 285, 973, 1053, 1047, 945, 976, 1177, 1182, 1046, 1320)
)

# Setting more instructive names for highlighted features
name_change <- c(
  "Feature269" = "Coumaric\nacid (F269)",
  "Feature455" = "Ferulic acid\n (F455)",
  "Feature285" = "Vanillic acid\n (F285)",
  "Feature973" = "Naringenin\n (F973)",
  "Feature1053" = "Vestitone\n (F1053)",
  "Feature1047" = "BiochaninA/\nOlmelin (F1047)",
  "Feature945" = "Formononetin\n (F945)",
  "Feature976" = "Vestitol\n (F976)",
  "Feature1177" = "Dehydroquer-\ncetin (F1177)",
  "Feature1182" = "Diosmetin\n (F1182)",
  "Feature1046" = "Wogonin\n (F1046)",
  "Feature1320" = "Velutin\n (F1320)"
)

# Setting factor levels
highlighted_features_Lj <- factor(
  highlighted_features_Lj,
  levels = highlighted_features_Lj
)

# Long form feature table
metabolite_long_Lj <- melt(
  data = metabolite_data_Lj,
  id.vars = 1,
  variable.name = "Sample_ID",
  value.name = "Intensity"
)

# Keeping only highligted features in long form feature table
metabolite_long_Lj <- metabolite_long_Lj[Feature %in% highlighted_features_Lj]

# Adding metadata to long feature table
metabolite_long_Lj <- merge(metabolite_long_Lj, design_Lj, "Sample_ID")

# Loading p-values from genotype effects
p_vals_dt <- fread("../4_genotype_effects/1_tables/p_adj_Lj.csv")

# Cleaning up column names
colnames(p_vals_dt) <- gsub("_p_adj", "", colnames(p_vals_dt))

# Long form
p_vals_dt <- melt(
  data = p_vals_dt,
  id.vars = 1,
  variable.name = "Genotype",
  value.name = "p_adj"
)

# Keep only highligted features
p_vals_dt <- p_vals_dt[Feature %in% highlighted_features_Lj]

# Adding significance indicator (*)
p_vals_dt[,text:=ifelse(p_adj < 0.05, "*", "")]

# Adding y-position for significance indicator
y_pos <- metabolite_long_Lj[, .(y_pos = max(Intensity)), list(Feature, Genotype)]
p_vals_dt <- merge(p_vals_dt, y_pos, by = c("Feature", "Genotype"))

# Changing feature names according to name_change vector and setting factor levels
metabolite_long_Lj[, Feature := name_change[as.character(Feature)]]
metabolite_long_Lj[, Feature := factor(Feature, levels = name_change)]
p_vals_dt[, Feature := name_change[as.character(Feature)]]
p_vals_dt[, Feature := factor(Feature, levels = name_change)]

# Adjusting y-position of significance indicator
p_vals_dt[, y_pos := y_pos * 1.05]

# Setting factor level
metabolite_long_Lj[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

# Plot
box_plot_highlighted_Lj <- ggplot(
  metabolite_long_Lj,
  aes(x = Genotype, y = Intensity, fill = Genotype)
) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5, outlier.color = "red") +
  facet_wrap(~Feature, scales = "free", ncol = 3) +
  geom_label(
    data = p_vals_dt,
    aes(y = y_pos, label = text),
    label.size = NA,
    alpha = 0,
    size = 20 / .pt
  ) +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  theme_bw() +
  ggtitle("Lotus") +
  theme(
    legend.position = "right",
    strip.background = element_rect(colour = NA),
    axis.title.y = element_text(hjust = 0.2),
    axis.title = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6, colour = "black"),
    strip.text = element_text(colour = 'black', size = 6, face = "bold"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
    plot.margin = margin(l = 0.1, r = 0.1, t = 0.5, b = 0.5, unit = "line"),
    legend.margin = margin(t = 20, unit = "pt")
  ) +
  labs(x = NULL, title = "Lotus") +
  guides(fill = "none") +
  facetted_pos_scales(
    y = list(
      Feature == "BiochaninA/\nOlmelin (F1047)" ~
        scale_y_continuous(limits = c(0, 25000))
    )
  ) +
  NULL

## Hordeum --------------------------------------------------------------------
highlighted_features_Hv <- paste0(
  "Feature",
  c(2546, 2889, 495, 3095, 3069, 3288)
)

# Setting more instructive names for highlighted features
name_change <- c(
  "Feature2546" = "Gibberellin\n (F2546)",
  "Feature2889" = "Abscisic \nacid (F2889)",
  "Feature495" = "Esculetin\n (F495)",
  "Feature3095" = "Paeonin C\n (F3095)",
  "Feature3069" = "Isoorientin\n (F3069)",
  "Feature3288" = "Sapopharin\n (F3288)"
)

# Setting factor levels
highlighted_features_Hv <- factor(
  highlighted_features_Hv,
  levels = highlighted_features_Hv
)

# Long form feature table
metabolite_long_Hv <- melt(
  data = metabolite_data_Hv,
  id.vars = 1,
  variable.name = "Sample_ID",
  value.name = "Intensity"
)

# Keeping only highligted features in long form feature table
metabolite_long_Hv <- metabolite_long_Hv[Feature %in% highlighted_features_Hv]

# Adding metadata to long feature table
metabolite_long_Hv <- merge(metabolite_long_Hv, design_Hv, "Sample_ID")

# Loading p-values from genotype effects
p_vals_dt <- fread("../4_genotype_effects/1_tables/p_adj_Hv.csv")

# Cleaning up column names
colnames(p_vals_dt) <- gsub("_p_adj", "", colnames(p_vals_dt))

# Long form
p_vals_dt <- melt(
  data = p_vals_dt,
  id.vars = 1,
  variable.name = "Genotype",
  value.name = "p_adj"
)

# Keep only highligted features
p_vals_dt <- p_vals_dt[Feature %in% highlighted_features_Hv]

# Adding significance indicator (*)
p_vals_dt[, text := ifelse(p_adj < 0.05, "*", "")]

# Adding y-position for significance indicator
y_pos <- metabolite_long_Hv[,
  .(y_pos = max(Intensity)),
  list(Feature, Genotype)
]
p_vals_dt <- merge(p_vals_dt, y_pos, by = c("Feature", "Genotype"))

# Changing feature names according to name_change vector and setting factor levels
metabolite_long_Hv[, Feature := name_change[as.character(Feature)]]
metabolite_long_Hv[, Feature := factor(Feature, levels = name_change)]
p_vals_dt[, Feature := name_change[as.character(Feature)]]
p_vals_dt[, Feature := factor(Feature, levels = name_change)]

# Adjusting y-position of significance indicator
p_vals_dt[, y_pos := y_pos * 1.05]

# Setting factor level
metabolite_long_Hv[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

# Plot
box_plot_highlighted_Hv <- ggplot(
  metabolite_long_Hv,
  aes(x = Genotype, y = Intensity, fill = Genotype)
) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5, outlier.color = "red") +
  facet_wrap(~Feature, scales = "free", ncol = 3) +
  geom_label(
    data = p_vals_dt,
    aes(y = y_pos, label = text),
    label.size = NA,
    alpha = 0,
    size = 20 / .pt
  ) +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  theme_bw() +
  ggtitle("Hordeum") +
  theme(
    legend.position = "right",
    strip.background = element_rect(colour = NA),
    axis.title = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6, colour = "black"),
    strip.text = element_text(colour = 'black', size = 6, face = "bold"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
    plot.margin = margin(l = 0.1, r = 0.1, t = 0.5, b = 0.5, unit = "line"),
    legend.margin = margin(t = 20, unit = "pt")
  ) +
  labs(x = NULL, y = " ") +
  guides(fill = "none") +
  facetted_pos_scales(
    y = list(
      Feature == "Abscisic \nacid (F2889)" ~
        scale_y_continuous(limits = c(0, 2500))
    )
  ) +
  NULL

box_plot_highlighted <- plot_grid(
  box_plot_highlighted_Lj, box_plot_highlighted_Hv,
  ncol = 1,
  rel_heights = c(2/3, 1/3)
)

saveRDS(
  object = box_plot_highlighted,
  file = "1_rds_files/box_plot_highlighted.rds"
)

ggsave(
  filename = "2_figures/box_plot_highlighted.pdf",
  plot = box_plot_highlighted,
  width = 9,
  height = 18,
  units = "cm"
)
