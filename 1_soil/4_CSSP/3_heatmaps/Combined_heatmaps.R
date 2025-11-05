pkg <- c("data.table", "magrittr", "ggplot2", "vegan", "Maaslin2", "patchwork",
         "RColorBrewer", "colorRamp2", "ggh4x", "ggtext", "ggpubr", "cowplot")

for(pk in pkg){
  library(pk, character.only = T)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Lotus plot
Axis_title <- readRDS("1_Lotus/Axis_title_Lj.rds")
Barplot_Lj <- readRDS("1_Lotus/RA_Barplot_Lj.rds")
Taxonomy_Lj <- readRDS("1_Lotus/Taxonomy_Lj.rds")
Heatmap_Lj <- readRDS("1_Lotus/Heatmap_Lj.rds")

top_plot_Lj <- (Taxonomy_Lj | Axis_title) + plot_layout(widths = c(0.9, 0.1))
bottom_plot_Lj <- (Heatmap_Lj | Barplot_Lj) + plot_layout(widths = c(0.9, 0.1))
main_plot_Lj <- (top_plot_Lj/bottom_plot_Lj) + 
  plot_layout(heights = c(0.05, 0.95))
main_plot_Lj

# Hordeum plot plot
empty <- ggplot() + theme_void()
Barplot_Hv <- readRDS("2_Hordeum/RA_Barplot_Hv.rds")
Taxonomy_Hv <- readRDS("2_Hordeum/Taxonomy_with_legend_Hv.rds")
Heatmap_Hv <- readRDS("2_Hordeum/Heatmap_with_legend_Hv.rds")

tax_leg <-  ggpubr::get_legend(Taxonomy_Hv)
heat_leg <- ggpubr::get_legend(Heatmap_Hv)

Taxonomy_Hv <- Taxonomy_Hv + theme(legend.position = "none")
Heatmap_Hv <- Heatmap_Hv + theme(legend.position = "none")

top_plot_Hv <- (Taxonomy_Hv | empty) + plot_layout(widths = c(0.9, 0.1))
bottom_plot_Hv <- (Heatmap_Hv | Barplot_Hv) + plot_layout(widths = c(0.9, 0.1))
main_plot_Hv <- (top_plot_Hv/bottom_plot_Hv) + 
  plot_layout(heights = c(0.05, 0.95))
main_plot_Hv

# Combined
combined_plot <- (
  (Taxonomy_Lj | Axis_title) + plot_layout(widths = c(0.9, 0.1))
)/(
  (Heatmap_Lj | Barplot_Lj) + plot_layout(widths = c(0.9, 0.1))
)/(
  (Taxonomy_Hv | empty) + plot_layout(widths = c(0.9, 0.1))
)/(
  (Heatmap_Hv | Barplot_Hv) + plot_layout(widths = c(0.9, 0.1))
) + plot_layout(heights = c(0.025, 0.475, 0.025, 0.475))

plot_legend <- plot_grid(tax_leg, heat_leg, rel_widths = c(0.8, 0.2))

combined_plot_with_legend <- plot_grid(
  combined_plot, plot_legend, 
  rel_heights = c(0.9, 0.1),
  ncol = 1
)

ggsave("Heatmaps_combined.pdf", combined_plot_with_legend,
       width = 210, height = 180, units = "mm")

saveRDS(combined_plot_with_legend, "Combined_heatmap.rds")
