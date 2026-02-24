# Clean up.
options(warn=-1)
rm(list=ls())

## Load required packages.
pkg <- c("data.table", "magrittr", "ggplot2", "vegan", "Maaslin2", "patchwork",
         "RColorBrewer", "colorRamp2", "ggh4x", "ggtext", "ggpubr", "cowplot")

for(pk in pkg){
  library(pk, character.only = T)
}

## Set directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Make Lotus plot.
Axis_title <- readRDS("1_Lotus/LotusCSSP_Askov_DA_axis_title.rds")
Lotus_barplot <- readRDS("1_Lotus/LotusCSSP_Askov_DA_barplot.rds")
Lotus_taxonomy <- readRDS("1_Lotus/LotusCSSP_Askov_DA_taxonomy.rds")
Lotus_heatmap <- readRDS("1_Lotus/LotusCSSP_Askov_DA_heatmap.rds")

Lotus_top_plot <- (Lotus_taxonomy | Axis_title) + plot_layout(widths = c(0.9, 0.1))
Lotus_bottom_plot <- (Lotus_heatmap | Lotus_barplot) + plot_layout(widths = c(0.9, 0.1))
Lotus_main_plot <- (Lotus_top_plot/Lotus_bottom_plot) + 
  plot_layout(heights = c(0.05, 0.95))
Lotus_main_plot

# Make Hordeum plot.
empty <- ggplot() + theme_void()
Hordeum_barplot <- readRDS("2_Hordeum/HordeumCSSP_Askov_DA_barplot.rds")
Hordeum_taxonomy <- readRDS("2_Hordeum/HordeumCSSP_Askov_DA_taxonomy.rds")
Hordeum_heatmap <- readRDS("2_Hordeum/HordeumCSSP_Askov_DA_heatmap.rds")

tax_leg <-  ggpubr::get_legend(Hordeum_taxonomy)
heat_leg <- ggpubr::get_legend(Hordeum_heatmap)

Hordeum_taxonomy <- Hordeum_taxonomy + theme(legend.position = "none")
Hordeum_heatmap <- Hordeum_heatmap + theme(legend.position = "none")

Hordeum_top_plot <- (Hordeum_taxonomy | empty) + plot_layout(widths = c(0.9, 0.1))
Hordeum_bottom_plot <- (Hordeum_heatmap | Hordeum_barplot) + plot_layout(widths = c(0.9, 0.1))
Hordeum_main_plot <- (Hordeum_top_plot/Hordeum_bottom_plot) + 
  plot_layout(heights = c(0.05, 0.95))
Hordeum_main_plot

# Combine both plots.
combined_plot <- (
  (Lotus_taxonomy | Axis_title) + plot_layout(widths = c(0.9, 0.1))
)/(
  (Lotus_heatmap | Lotus_barplot) + plot_layout(widths = c(0.9, 0.1))
)/(
  (Hordeum_taxonomy | empty) + plot_layout(widths = c(0.9, 0.1))
)/(
  (Hordeum_heatmap | Hordeum_barplot) + plot_layout(widths = c(0.9, 0.1))
) + plot_layout(heights = c(0.025, 0.475, 0.025, 0.475))

plot_legend <- plot_grid(tax_leg, heat_leg, rel_widths = c(0.8, 0.2))

combined_plot_with_legend <- plot_grid(
  combined_plot, plot_legend, 
  rel_heights = c(0.86, 0.14),
  ncol = 1
)

# Combine final plot.
ggsave("LotusHordeum_CSSP_Askov_DA_figure_combined.pdf", combined_plot_with_legend,
       width = 210, height = 180, units = "mm")

saveRDS(combined_plot_with_legend, "LotusHordeum_CSSP_Askov_DA_figure_combined.rds")
saveRDS(combined_plot_with_legend, "../8_final_figures/LotusHordeum_CSSP_Askov_DA_figure_combined.rds")
