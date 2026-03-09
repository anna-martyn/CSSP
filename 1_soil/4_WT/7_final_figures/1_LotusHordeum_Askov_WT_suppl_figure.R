# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
pkg <- c("ggplot2", "forcats", "cowplot", "patchwork", "magick", "grid")
for(pk in pkg){
  library(pk, character.only = T)
}

# Read all plot files (rds files) into variables.
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)  # remove .rds extension
  assign(plot_name, readRDS(f))
}

# Load the plots wanted for the final figure.
p1 <- LotusCSSP_AskovSoils_WT_chao1_rfd
p2 <- Hordeum_AskovSoils_WT_chao1_rfd
p3 <- LotusHordeum_Askov_WT_stackedbp_top20_meanRA
p4 <- Lotus_Askov_WT_stackedbp_NoduleASVs

# Define the individual rows for the final figure.
row1 <- plot_grid(
  p1, p2,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 15,
  rel_widths = c(3/5, 2/5)
)

row2 <- plot_grid(
  p3, p4,
  ncol = 2,
  labels = c("C", "D"),
  label_size = 15,
  rel_widths = c(2/3, 1/3)
)

# Combine rows in final figure.
final_plot <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(0.5, 0.5)
)

final_plot

# Save the combined plot as a PDF file.
ggsave("Suppl_Figure1_Askov_WT.pdf", final_plot, width=21, height=14, unit="cm")
