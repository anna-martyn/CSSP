# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library(ggplot2)
library(patchwork)
library(cowplot)
library(magick)
library(grid)
library(forcats)

# Read all plot files (rds files) into variables.
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)
  assign(plot_name, readRDS(f))
}

# Load the plots wanted for the supplementary figures.
p1 <- LotusCSSP_AskovSoils_chao1_rfd
genotype_labels <- c(
  "WT"    = "WT",
  "symrk" = "<i>symrk</i>",
  "ccamk" = "<i>ccamk</i>",
  "nsp1"  = "<i>nsp1</i>",
  "nsp2"  = "<i>nsp2</i>"
)
p1 <- p1 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

p2 <- HordeumCSSP_AskovSoils_chao1_rfd
p2 <- p2 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())


p3 <- LotusCSSP_AskovSoils_cpcoa_all
p4 <- HordeumCSSP_AskovSoils_cpcoa_all
p5 <- LotusCSSP_AskovSoils_pcoa_all
p6 <- HordeumCSSP_AskovSoils_pcoa_all
p7 <- LotusHordeum_Askov_stackedbp_top20_meanRA
p7 <- p7 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

p8 <- LotusHordeum_Askov_orders_heatmap
p8 <- p8 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

cpcoa_legend <- ggpubr::get_legend(HordeumCSSP_AskovSoils_cpcoa_with_legend)
pcoa_legend <- ggpubr::get_legend(
  HordeumCSSP_AskovSoils_cpcoa_with_legend +
    guides(color = guide_legend(title = "Genotype"),
           shape = "none",
           fill = "none")
)

# Define the individual rows for supplementary figure 1.
row1 <- plot_grid(
  p1, p2,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 15,
  rel_widths = c(1/2, 1/2)
)

row2 <- plot_grid(
  p3, p4, cpcoa_legend,
  ncol = 3,
  labels = c("C", "D", ""),
  label_size = 15,
  rel_widths = c(5/11, 5/11, 1/11)
)

row3 <- plot_grid(
  p5, p6, pcoa_legend,
  ncol = 3,
  labels = c("E", "F"),
  label_size = 15,
  rel_widths = c(5/11, 5/11, 1/11)
)

# Make combined supplementary figure 1.
final_plot <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1/4, 1/4, 2/4)
)

final_plot

# Save supplementary figure 1 as a PDF file.
ggsave("Suppl_Figure2_Askov_CSSP.pdf", final_plot, width=21, height=14, unit="cm")

# Define the individual rows for supplementary figure 2.
row4 <- plot_grid(
  p7,
  ncol = 1,
  labels = c("A"),
  label_size = 15,
  rel_widths = c(1/1)
)

row5 <- plot_grid(
  p8,
  ncol = 1,
  labels = c("B"),
  label_size = 15,
  rel_widths = c(1/1)
)

# Make combined supplementary figure 2.
final_plot2 <- plot_grid(
  row4,
  row5,
  ncol = 1,
  rel_heights = c(1/3, 2/3)
)

final_plot2

# Save supplementary figure 1 as a PDF file.
ggsave("Suppl_Figure3_Askov_CSSP.pdf", final_plot2, width=21, height=14, unit="cm")
