# Cleaning up
options(warn = -1)
rm(list = ls())

# Loading packages
pkg <- c("ggplot2", "patchwork", "ggh4x", "ggtext", "ggpubr")

for(pk in pkg){
  library(pk, character.only = T)
}

# Setting directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Lotus plot
axis_title <- readRDS("1_Lotus/LotusCSSP_Askov_DA_axis_title.rds")
lotus_barplot <- readRDS("1_Lotus/LotusCSSP_Askov_DA_barplot.rds")
lotus_taxonomy <- readRDS("1_Lotus/LotusCSSP_Askov_DA_taxonomy.rds")
lotus_heatmap <- readRDS("1_Lotus/LotusCSSP_Askov_DA_heatmap.rds")

lotus_top_plot <- (lotus_taxonomy | axis_title) +
  plot_layout(widths = c(0.9, 0.1))
lotus_bottom_plot <- (lotus_heatmap | lotus_barplot) +
  plot_layout(widths = c(0.9, 0.1))
lotus_main_plot <- (lotus_top_plot / lotus_bottom_plot) +
  plot_layout(heights = c(0.05, 0.95))
lotus_main_plot

# Hordeum plot
empty <- ggplot() + theme_void()
hordeum_barplot <- readRDS("2_Hordeum/HordeumCSSP_Askov_DA_barplot.rds")
hordeum_taxonomy <- readRDS("2_Hordeum/HordeumCSSP_Askov_DA_taxonomy.rds")
hordeum_heatmap <- readRDS("2_Hordeum/HordeumCSSP_Askov_DA_heatmap.rds")

tax_leg <-  get_legend(hordeum_taxonomy)
heat_leg <- get_legend(hordeum_heatmap)

hordeum_taxonomy <- hordeum_taxonomy + theme(legend.position = "none")
hordeum_heatmap <- hordeum_heatmap + theme(legend.position = "none")

hordeum_top_plot <- (hordeum_taxonomy | empty) +
  plot_layout(widths = c(0.9, 0.1))
hordeum_bottom_plot <- (hordeum_heatmap | hordeum_barplot) +
  plot_layout(widths = c(0.9, 0.1))
hordeum_main_plot <- (hordeum_top_plot / hordeum_bottom_plot) +
  plot_layout(heights = c(0.05, 0.95))
hordeum_main_plot

# Combined plot
combined_plot <- ((lotus_taxonomy | axis_title) +
  plot_layout(widths = c(0.9, 0.1))) /
  ((lotus_heatmap | lotus_barplot) + plot_layout(widths = c(0.9, 0.1))) /
  ((hordeum_taxonomy | empty) + plot_layout(widths = c(0.9, 0.1))) /
  ((hordeum_heatmap | hordeum_barplot) + plot_layout(widths = c(0.9, 0.1))) +
  plot_layout(heights = c(0.025, 0.475, 0.025, 0.475))

plot_legend <- ggarrange(tax_leg, heat_leg, widths = c(0.8, 0.2))

combined_plot_with_legend <- ggarrange(
  combined_plot, plot_legend, 
  heights = c(0.84, 0.16),
  ncol = 1
)

# Saving plot
ggsave(
  filename = "4_figures/LotusHordeum_CSSP_Askov_DA_figure_combined2.pdf",
  plot = combined_plot_with_legend,
  width = 210,
  height = 180,
  units = "mm"
)

saveRDS(
  object = combined_plot_with_legend,
  file = "3_rds_files/LotusHordeum_CSSP_Askov_DA_figure_combined.rds"
)
