# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
library(ggplot2)
library(patchwork)
library(cowplot)
library(magick)
library(grid)

# Read all plot files (rds files) into variables
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)  # remove .rds extension
  assign(plot_name, readRDS(f))
}

# Load the plots.
p1 <- Askov_Lotus_WT_chao1_rfd
p2 <- Askov_Hordeum_WT_chao1_rfd
p3 <- Soil_LotusWT_stackedbp_NoduleASVs

# 
final_plot <- plot_grid(
  p1, p2, p3,
  ncol = 3,
  rel_widths = c(1.2, 0.9, 1.2),
  labels = c("A","B","C"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

final_plot

# Save the combined plot as PDF
ggsave("Soil_WT_suppl_figure.pdf", final_plot, width=21, height=10)
