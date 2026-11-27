# Setup ------------------------------------------------------------------------
# Rarefaction curves (observed ASVs vs. sequencing depth) for the four datasets:
#   - Lotus soil (spike-in removed ASV table)
#   - Hordeum soil
#   - Lotus SynCom (ASVs matched to the SynCom only)
#   - Hordeum SynCom (ASVs matched to the SynCom only)
# Curves are computed on the non-rarefied ASV tables. The dashed line in each
# panel marks the depth of the rarefied ASV table on disk for that dataset.

# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "vegan")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Functions --------------------------------------------------------------------
# Reading an ASV table (ASVs as rows, samples as columns)
read_asv <- function(file, biom_header = FALSE){
  if(biom_header){
    # Tables exported from biom files start with a comment line
    read.table(
      file = file,
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE,
      comment.char = "",
      skip = 1
    )
  } else {
    read.table(
      file = file,
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )
  }
}

# Reading metadata
read_design <- function(file){
  read.table(
    file = file,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
}

# Depth of a rarefied ASV table (all samples have the same depth)
rarefied_depth <- function(file){
  depths <- unique(colSums(read_asv(file)))
  stopifnot(length(depths) == 1)
  return(depths)
}

# Rarefaction curves for one dataset, with metadata
curves <- function(asv_table, design, dataset){
  rc <- rarecurve(t(asv_table), step = 250, tidy = TRUE)
  colnames(rc) <- c("SampleID", "Reads", "ASVs")
  rc$SampleID <- as.character(rc$SampleID)
  rc$Compartment <- design[rc$SampleID, "Compartment"]
  rc$Dataset <- dataset
  return(rc)
}

# Loading data -----------------------------------------------------------------
# Lotus soil (spike-in removed)
lotus_soil <- read_asv(
  "../1_soil/1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
lotus_soil_design <- read_design(
  "../1_soil/1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt"
)
lotus_soil_depth <- rarefied_depth(
  "../1_soil/2_rarefication/1_Lotus/LotusCSSP_AskovSoils_ASVtable10_4_rfd.txt"
)

# Hordeum soil
hordeum_soil <- read_asv(
  "../1_soil/1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  biom_header = TRUE
)
hordeum_soil_design <- read_design(
  "../1_soil/1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt"
)
hordeum_soil_depth <- rarefied_depth(
  "../1_soil/2_rarefication/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4_rfd.txt"
)

# Lotus SynCom, keeping only ASVs matched to the SynCom (as in
# 3_syncoms/1_Lotus/4_rarefication/0_LotusSC_rarefication.R)
lotus_sc <- read_asv("../3_syncoms/1_Lotus/1_data/LotusSC_ASVtable_nospike.tsv")
lotus_sc <- lotus_sc[grepl("Lj", rownames(lotus_sc)), ]
lotus_sc_design <- read_design("../3_syncoms/1_Lotus/1_data/LotusSC_metadata.txt")
lotus_sc_depth <- rarefied_depth(
  "../3_syncoms/1_Lotus/4_rarefication/1_tables/LotusSC_rfd_nocontaminants.txt"
)

# Hordeum SynCom, keeping only ASVs matched to the SynCom (as in
# 3_syncoms/2_Hordeum/3_rarefication/HordeumSynCom_rarefication.R)
hordeum_sc <- read_asv(
  "../3_syncoms/2_Hordeum/1_data/HordeumSC_ASVtable.tsv",
  biom_header = TRUE
)
hordeum_sc <- hordeum_sc[grepl("_", rownames(hordeum_sc)), ]
hordeum_sc_design <- read_design("../3_syncoms/2_Hordeum/1_data/HordeumSC_metadata.txt")
hordeum_sc_depth <- rarefied_depth(
  "../3_syncoms/2_Hordeum/3_rarefication/1_tables/HordeumSC_rfd_nocontaminants.txt"
)

# Rarefaction curves -----------------------------------------------------------
datasets <- c("Lotus soil", "Hordeum soil", "Lotus SynCom", "Hordeum SynCom")

rare_curves <- rbind(
  curves(lotus_soil, lotus_soil_design, "Lotus soil"),
  curves(hordeum_soil, hordeum_soil_design, "Hordeum soil"),
  curves(lotus_sc, lotus_sc_design, "Lotus SynCom"),
  curves(hordeum_sc, hordeum_sc_design, "Hordeum SynCom")
)

# Rarefaction depths of the rarefied ASV tables on disk
thresholds <- data.frame(
  Dataset = datasets,
  Depth = c(lotus_soil_depth, hordeum_soil_depth, lotus_sc_depth, hordeum_sc_depth)
)
print(thresholds)

# Setting factor levels
rare_curves$Dataset <- factor(rare_curves$Dataset, levels = datasets)
thresholds$Dataset <- factor(thresholds$Dataset, levels = datasets)
rare_curves$Compartment <- factor(
  rare_curves$Compartment,
  levels = c("Bulk", "Rhizosphere", "Root", "Nodules")
)

# Compartment colours
colors <- c(
  "Bulk" = "gray40",
  "Rhizosphere" = "#B2563C",
  "Root" = "#3C7D82",
  "Nodules" = "#6F944F"
)

# Plot -------------------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  axis.title = element_text(size = 6),
  legend.text = element_text(size = 6, color = "black"),
  legend.title = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  legend.position = "bottom",
  strip.text = element_text(size = 6, color = "black", face = "bold"),
  legend.background = element_blank()
)

# Rarefaction curves, one line per sample and one panel per dataset.
# Dashed line marks the depth of the rarefied ASV table on disk.
rare_plot <- function(data){
  ggplot(data, aes(x = Reads, y = ASVs, group = SampleID, color = Compartment)) +
    geom_line(linewidth = 0.3, alpha = 0.7) +
    geom_vline(
      data = thresholds,
      aes(xintercept = Depth),
      linetype = "dashed",
      linewidth = 0.3
    ) +
    geom_text(
      data = thresholds,
      aes(x = Depth, y = Inf, label = paste0(" ", format(Depth, big.mark = ","), " reads")),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1.5,
      size = 6 / .pt
    ) +
    scale_color_manual(values = colors, drop = TRUE) +
    guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) +
    scale_x_continuous(labels = scales::label_comma()) +
    facet_wrap(~Dataset, scales = "free") +
    labs(x = "Number of reads", y = "Observed ASVs", color = "") +
    main_theme
}

# Full curves
rare_plot_full <- rare_plot(rare_curves)

# Zoomed version on the first 25,000 reads, where most curves plateau
rare_plot_zoom <- rare_plot(rare_curves[rare_curves$Reads <= 25000, ])

# Saving -----------------------------------------------------------------------
dir.create("1_rds_files", showWarnings = FALSE)
dir.create("2_figures", showWarnings = FALSE)

saveRDS(rare_plot_full, file = "1_rds_files/rarefaction_curves_combined.rds")
saveRDS(rare_plot_zoom, file = "1_rds_files/rarefaction_curves_combined_zoom.rds")

ggsave(
  filename = "2_figures/rarefaction_curves_combined.pdf",
  plot = rare_plot_full,
  width = 18,
  height = 14,
  units = "cm"
)
ggsave(
  filename = "2_figures/rarefaction_curves_combined_zoom.pdf",
  plot = rare_plot_zoom,
  width = 18,
  height = 14,
  units = "cm"
)
ggsave(
  filename = "Suppl_Figure_X_Rarefaction_curves.pdf",
  plot = rare_plot_zoom,
  width = 18,
  height = 14,
  units = "cm"
)
