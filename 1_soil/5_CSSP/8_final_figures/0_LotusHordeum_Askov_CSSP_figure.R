# Cleaning up
options(warn = -1)
rm(list = ls())

# Loading packages
pkg <- c("data.table", "ggplot2", "vegan", "ggtext", "ggpubr", "grid", "cowplot")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading individual files, legends etc. for plot
empty <- ggplot() + theme_void()
Lotus_boxplot <- readRDS(
  "../1_shoot_fw/1_rds_files/LotusCSSP_AskovSoils_shootfw.rds"
)
Hordeum_boxplot <- readRDS(
  "../1_shoot_fw/1_rds_files/HordeumCSSP_AskovSoils_shootfw.rds"
)

Lotus_CPCoA_PCoA <- readRDS(
  "../5_cpcoa_pcoa/1_rds_files/LotusCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.rds"
)
Hordeum_CPCoA_PCoA <- readRDS(
  "../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.rds"
)
Hordeum_CPCoA <- readRDS(
  "../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_cpcoa_with_legend.rds"
)

legend <- ggpubr::get_legend(Hordeum_CPCoA)

Heatmap <- readRDS(
  "../7_DA_analysis/3_rds_files/LotusHordeum_CSSP_Askov_DA_figure_combined.rds"
)

# Defining individual rows for final figure
row1 <- plot_grid(
  empty, Lotus_boxplot, Hordeum_boxplot,
  ncol = 3,
  rel_widths = c(0.33, 0.34, 0.32),
  labels = c("A","B",""),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row2 <- plot_grid(
  Lotus_CPCoA_PCoA, Hordeum_CPCoA_PCoA, legend,
  ncol = 3,
  rel_widths = c(0.45, 0.45, 0.1),
  labels = c("C","",""),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row3 <- Heatmap

# Combining rows in final figure
final_plot <- plot_grid(
  row1, row2, row3,
  ncol = 1,
  labels = c("","","D"),
  label_size = 12,
  label_fontface = "bold",
  rel_heights = c(0.2, 0.2, 0.6)
)

# Saving final figure
ggsave(
  filename = "Figure4_Askov_CSSP.pdf",
  plot = final_plot,
  width = 18,
  height = 22,
  units = "cm"
)
