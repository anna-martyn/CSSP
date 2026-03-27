# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
pkg <- c("ggplot2", "grid", "cowplot", "patchwork", "magick", "forcats")
for(pk in pkg){
  library(pk, character.only = T)
}

# Read all plot files (rds files) into variables.
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)
  assign(plot_name, readRDS(f))
}

# Load the plots wanted for the supplementary figures.
p1 <- readRDS("../3_nod_cts/1_rds_files/LotusCSSP_AskovSoils_nod_cts.rds")

p2 <- readRDS(
  "../2_shoot_fw_axenic/1_rds_files/HordeumCSSP_axenic_shoot_fw.rds"
)

p3 <- readRDS(
  "../4_chao1/1_Lotus/1_rds_files/LotusCSSP_AskovSoils_chao1_rfd.rds"
)
genotype_labels <- c(
  "WT"    = "WT",
  "symrk" = "<i>symrk</i>",
  "ccamk" = "<i>ccamk</i>",
  "nsp1"  = "<i>nsp1</i>",
  "nsp2"  = "<i>nsp2</i>"
)
p3 <- p3 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

p4 <- readRDS(
  "../4_chao1/2_Hordeum/1_rds_files/HordeumCSSP_AskovSoils_chao1_rfd.rds"
)
p4 <- p4 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())


p5 <- readRDS("../5_cpcoa_pcoa/1_rds_files/LotusCSSP_AskovSoils_cpcoa_all.rds")
p6 <- readRDS(
  "../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_cpcoa_all.rds"
)
p7 <- readRDS("../5_cpcoa_pcoa/1_rds_files/LotusCSSP_AskovSoils_pcoa_all.rds")
p8 <- readRDS("../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_pcoa_all.rds")
p9 <- readRDS(
  "../6_stackedbp_heatmap/1_rds_files/LotusHordeum_Askov_stackedbp_top20_meanRA.rds"
)

p9 <- p9 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

p10 <- readRDS(
  "../6_stackedbp_heatmap/1_rds_files/LotusHordeum_Askov_orders_heatmap.rds"
)

p10 <- p10 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

legend_plot <- readRDS(
  "../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_cpcoa_with_legend.rds"
)
cpcoa_legend <- ggpubr::get_legend(legend_plot)
pcoa_legend <- ggpubr::get_legend(
  legend_plot +
    guides(color = guide_legend(title = "Genotype"),
           shape = "none",
           fill = "none")
)

# Define the individual rows for supplementary figure 2.
# Make combined supplementary figure 2.
row1 <- plot_grid(
  p1, p2,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 15,
  rel_widths = c(2/3, 1/3)
)

final_plot <- plot_grid(
  row1,
  # row2,
  # row3,
  ncol = 1,
  rel_heights = c(1/1)
)

final_plot

# Define the individual rows for supplementary figure 3.
row1 <- plot_grid(
  p3, p4,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 15,
  rel_widths = c(1/2, 1/2)
)

row2 <- plot_grid(
  p5, p6, cpcoa_legend,
  ncol = 3,
  labels = c("C", "D", ""),
  label_size = 15,
  rel_widths = c(5/11, 5/11, 1/11)
)

row3 <- plot_grid(
  p7, p8, pcoa_legend,
  ncol = 3,
  labels = c("E", "F"),
  label_size = 15,
  rel_widths = c(5/11, 5/11, 1/11)
)

# Make combined supplementary figure 3.
final_plot2 <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1.5/4, 0.75/4, 1.75/4)
)

final_plot2

# Define the individual rows for supplementary figure 4.
row4 <- plot_grid(
  p9,
  ncol = 1,
  labels = c("A"),
  label_size = 15,
  rel_widths = c(1/1)
)

row5 <- plot_grid(
  p10,
  ncol = 1,
  labels = c("B"),
  label_size = 15,
  rel_widths = c(1/1)
)

# Make combined supplementary figure 4.
final_plot3 <- plot_grid(
  row4,
  row5,
  ncol = 1,
  rel_heights = c(0.5, 0.5)
)

final_plot3

# Save supplementary figures as PDF files.
ggsave(
  "Suppl_Figure2_Askov_CSSP.pdf",
  final_plot,
  width = 21,
  height = 8,
  unit = "cm"
)
ggsave(
  "Suppl_Figure3_Askov_CSSP.pdf",
  final_plot2,
  width = 21,
  height = 29.7,
  unit = "cm"
)
ggsave(
  "Suppl_Figure4_Askov_CSSP.pdf",
  final_plot3,
  width = 21,
  height = 24,
  unit = "cm"
)
