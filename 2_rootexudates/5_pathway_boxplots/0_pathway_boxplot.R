# Seup ------------------------------------------------------------------------
pkg <- c("data.table", "ggplot2", "ggtext", "ggh4x", "scales")
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

dir_name <- paste(
  "..",
  "1_data",
  "1_Lotus",
  "LotusCSSP_RootEx_Apr26_stdAUnew_canopus_structure_summary.tsv",
  sep = "/"
)
annotation_Lj <- fread(dir_name, drop = 21)

# Setting sample names in metadata
design_Lj[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Lj <- design_Lj[ Sample_ID %in% colnames(metabolite_data_Lj)]

# Setting feature names in annotation table
setnames(annotation_Lj, "mappingFeatureId", "Feature")
setcolorder(annotation_Lj, "Feature")
annotation_Lj[, Feature := paste0("Feature", Feature)]

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

dir_name <- paste(
  "..",
  "1_data",
  "2_Hordeum",
  "HordeumCSSP_RootEx_Apr26_stdAUnew_canopus_structure_summary.tsv",
  sep = "/"
)
annotation_Hv <- fread(dir_name, drop = 21)

# Setting sample names in metadata
design_Hv[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Hv <- design_Hv[ Sample_ID %in% colnames(metabolite_data_Hv)]

# Adding plant variable to hordeum metadata
design_Hv[,Plant := "Hordeum"]

# Setting feature names in annotation table
setnames(annotation_Hv, "mappingFeatureId", "Feature")
setcolorder(annotation_Hv, "Feature")
annotation_Hv[, Feature := paste0("Feature", Feature)]

# Pathway-level boxplot -------------------------------------------------------
## Lotus ----------------------------------------------------------------------
# Removing features from annotation not in feature table
annotation_Lj <- annotation_Lj[Feature %in% metabolite_data_Lj$Feature]

# Merging metabolite data and annotations
pathway_Lj <- merge(metabolite_data_Lj, annotation_Lj, by = "Feature")

# Removing features with pathway probability less than 60%
pathway_Lj <- pathway_Lj[`NPC#pathway Probability` > 0.6]

# Aggregating on pathway level
sample_names_Lj <- colnames(metabolite_data_Lj)[-1]
pathway_Lj <- pathway_Lj[,
  lapply(.SD, sum),
  by = `NPC#pathway`,
  .SDcols = sample_names_Lj
]

# Long format
pathway_Lj <- melt(
  pathway_Lj,
  id.vars = 1,
  variable.name = "Sample_ID",
  value.name = "Intensity"
)
setnames(pathway_Lj, "NPC#pathway", "Pathway")

# Merging with with metadata and setting factor levels
pathway_Lj <- merge(pathway_Lj, design_Lj)
pathway_Lj[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

# Linear regression on log-intensity for each pathway with Genotype as covariate
pathways <- unique(pathway_Lj$Pathway)
n_ptw <- length(pathways)
p_vals_lst <- list()
for(i in 1:n_ptw){
  pathway_Lj[Pathway == pathways[i]]
  l <- lm(log(Intensity)~Genotype, data = pathway_Lj[Pathway == pathways[i]])
  p_vals <- c(1, coef(summary(l))[-1, 4])
  p_vals_lst[[i]] <- data.table(
    Genotype = c("WT", "symrk", "ccamk", "nsp1", "nsp2"),
    Plant = "Lotus",
    Pathway = pathways[i],
    p_val = p_vals
  )
}

# Data table with p-values
p_vals_Lj <- rbindlist(p_vals_lst)

## Hordeum --------------------------------------------------------------------
# Removing features from annotation not in feature table
annotation_Hv <- annotation_Hv[Feature %in% metabolite_data_Hv$Feature]

# Merging metabolite data and annotations
pathway_Hv <- merge(metabolite_data_Hv, annotation_Hv, by = "Feature")

# Removing features with pathway probability less than 60%
# pathway_Hv <- pathway_Hv[`NPC#pathway Probability` < 0.6]

# Aggregating on pathway level
sample_names_Hv <- colnames(metabolite_data_Hv)[-1]
pathway_Hv <- pathway_Hv[,
  lapply(.SD, sum),
  by = `NPC#pathway`,
  .SDcols = sample_names_Hv
]

# Long format
pathway_Hv <- melt(
  pathway_Hv,
  id.vars = 1,
  variable.name = "Sample_ID",
  value.name = "Intensity"
)
setnames(pathway_Hv, "NPC#pathway", "Pathway")

# Merging with with metadata and setting factor levels
pathway_Hv <- merge(pathway_Hv, design_Hv)
pathway_Hv[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

# Linear regression on log-intensity for each pathway with Genotype as covariate
pathways <- unique(pathway_Lj$Pathway)
n_ptw <- length(pathways)
p_vals_lst <- list()
for(i in 1:n_ptw){
  pathway_Hv[Pathway == pathways[i]]
  l <- lm(log(Intensity)~Genotype, data = pathway_Hv[Pathway == pathways[i]])
  p_vals <- c(1, coef(summary(l))[-1, 4])
  p_vals_lst[[i]] <- data.table(
    Genotype = c("WT", "symrk", "ccamk", "nsp1", "nsp2"),
    Plant = "Hordeum",
    Pathway = pathways[i],
    p_val = p_vals
  )
}

# Data table with p-values
p_vals_Hv <- rbindlist(p_vals_lst)

# Combining Lotus and Hordeum p_val tables
p_vals_dt <- rbind(p_vals_Lj, p_vals_Hv)

# Correcting p-values for multiple testing
p_vals_dt[Genotype != "WT", p_adj := p.adjust(p_val)]
p_vals_dt[is.na(p_adj), p_adj := 1]
p_vals_dt[, Label := ifelse(p_adj < 0.05, "*", "")]
p_vals_dt[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

# Combining results for Lotus and Hordeum
setcolorder(pathway_Hv, colnames(pathway_Lj))
pathway_full <- rbind(pathway_Lj, pathway_Hv)
pathway_full <- merge(
  pathway_full,
  p_vals_dt,
  by = c("Genotype", "Plant", "Pathway")
)
pathway_full[, Plant := factor(Plant, levels = c("Lotus", "Hordeum"))]

# Setting y-position of significance indicator
letter_position <- pathway_full[,
  .(Label_pos = max(Intensity)),
  .(Pathway, Plant, Genotype)
]

# Adding y-positions to pathway 
pathway_full <- merge(
  pathway_full,
  letter_position,
  by = c("Genotype", "Plant", "Pathway")
)

# Adjusting y-positions separetely by plant
pathway_full[, Label_pos := Label_pos + 10000]

# Adding line breaks to some pathway names for plot
pathway_full[,
  Pathway := fcase(
    Pathway == "Shikimates and Phenylpropanoids" , "Shikimates\nand Phenyl-\npropanoids",
    Pathway == "Amino acids and Peptides", "Amino acids\n and Peptides",
    Pathway == "Carbohydrates", "Carbohy-\ndrates",
    default = Pathway
  )
]

# Plot
box_plot <- ggplot(
  pathway_full,
  aes(x = Genotype, y = Intensity, fill = Genotype)
) +
  geom_boxplot(width = 0.3, alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  geom_text(
    data = pathway_full,
    aes(x = Genotype, y = Label_pos, label = Label),
    inherit.aes = FALSE,
    size = 20 / .pt
  ) +
  facet_grid2(Plant ~ Pathway, scales = "free_y", independent = "y") +
  scale_fill_manual(values = cols, labels = legend_labels) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(colour = 'black', size = 6, face = "bold"),
    strip.background = element_rect(colour = NA),
    legend.position = "bottom",
    axis.title = element_text(size = 6),
    axis.text.x = element_text(
      size = 6,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.text.y = element_text(size = 6, colour = "black"),
    legend.text = element_markdown(
      size = 6,
      margin = margin(l = -0.1, unit = "pt")
    ),
    legend.title = element_text(size = 6),
    plot.title = element_text(size = 6, hjust = 0.5)
  ) +
  scale_x_discrete(
    labels = c(
      "symrk" = expression(italic("symrk")),
      "ccamk" = expression(italic("ccamk")),
      "nsp1" = expression(italic("nsp1")),
      "nsp2" = expression(italic("nsp2"))
    )
  ) +
  scale_y_continuous(labels = scientific)+
  expand_limits(y = 0)+
  NULL

# Sving figure
ggsave(
  "../9_final_figures/Suppl_Fig5.pdf",
  box_plot,
  width = 18,
  height = 12,
  units = "cm"
)
