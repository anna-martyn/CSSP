# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "cowplot")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading plots
p1 <- readRDS(
  "../3_chao1/1_Lotus/1_rds_files/LotusCSSP_AskovSoils_WT_chao1_rfd.rds"
)
p2 <- readRDS(
  "../3_chao1/2_Hordeum/1_rds_files/Hordeum_AskovSoils_WT_chao1_rfd.rds"
)
p3 <- readRDS(
  "../5_stacked_bp_heatmap/1_rds_files/LotusHordeum_Askov_WT_stackedbp_top20_meanRA.rds"
)
p4 <- readRDS(
  "../6_nod_ASVs/1_rds_files/Lotus_Askov_WT_stackedbp_NoduleASVs.rds"
)

# Defining rows of final figure
row1 <- plot_grid(
  p1, p2,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 8,
  rel_widths = c(3/5, 2/5)
)

row2 <- plot_grid(
  p3, p4,
  ncol = 2,
  labels = c("C", "D"),
  label_size = 8,
  rel_widths = c(2/3, 1/3)
)

# Combining rows in final figure
final_plot <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(0.5, 0.5)
)

# Saving combined plot
ggsave(
  filename = "Suppl_Figure1_Askov_WT.pdf",
  plot = final_plot,
  width = 18,
  height = 12,
  unit = "cm"
)
