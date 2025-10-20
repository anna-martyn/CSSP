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
p1 <- LotusSC_shootfw_incl_uninoc
p2 <- LotusSC_nodule_cts
p2_aligned <- p2 + labs(title = " ") +
  theme(plot.title = element_text(size = 20, color = NA))
p3 <- LotusSC_chao1_allASVs
p4 <- LotusSynCom_cpcoa_allASVs
p5 <- LotusSC_chao1_filteredASVs
p6 <- LotusSynCom_cpcoa_matchedASVsonly
p7 <- LotusSynCom_symbionts_RA

# Combine in one plot.
## Define individual rows.
row1 <- plot_grid(
  p1, p2_aligned,
  ncol = 2,
  rel_widths = c(2,1),
  labels = c("A","B"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row2 <- plot_grid(
  p3, p4,
  ncol = 2,
  rel_widths = c(2,1),
  labels = c("C","D"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row3 <- plot_grid(
  p5, p6,
  ncol = 2,
  rel_widths = c(2,1),
  labels = c("E","F"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row4 <- plot_grid(
  p7,
  ncol = 2,
  rel_widths = c(1,1),
  labels = c("G","", ""),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

## Combine all rows.
final_plot <- plot_grid(
  row1,
  row2,
  row3,
  row4,
  ncol = 1,
  rel_heights = c(1,1,1,1)
)

final_plot

# Save the combined plot as PDF
ggsave("LotusSC_suppl_figure.pdf", final_plot, width=21, height=29.7)
