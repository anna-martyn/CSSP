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
  plot_name <- tools::file_path_sans_ext(f)
  assign(plot_name, readRDS(f))
}

# Load plots.
p1 <- HordeumSC_Shootfw_incl_uninoc
p1_aligned <- p1 + labs(title = " ", subtitle="") + 
  theme(plot.title = element_text(size = 20, color = NA))
p2 <- HordeumSC_chao1_allASVs
p3 <- HordeumSynCom_cpcoa_allASVs
p4 <- HordeumSC_chao1_filteredASVs
p5 <- HordeumSynCom_cpcoa_matchedASVsonly

# Combine in one plot.
## Define individual rows.
row1 <- plot_grid(
  p1_aligned, p3, p5,
  ncol = 3,
  rel_widths = c(1,1,1),
  labels = c("A","B","C"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row2 <- plot_grid(
  p2, p4,
  ncol = 2,
  rel_widths = c(1,1),
  labels = c("D","E"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row3 <- empty_plot <- ggplot() + theme_void()

## Combine all rows.
final_plot <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1,1,1.4)
)

final_plot

# Save the combined plot as PDF
ggsave("HordeumSC_suppl_figure.pdf", final_plot, width=21, height=29.7)
