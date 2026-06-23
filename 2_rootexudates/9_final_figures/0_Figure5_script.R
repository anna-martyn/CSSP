pkg <- c("data.table", "ggplot2", "cowplot")
for(pk in pkg){
  library(pk, character.only = T)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading legends from rds files
pca_legend <- readRDS("../3_PCA/1_rds_files/PCA_plot_legend.rds")
volcano_legend <- readRDS("../6_volcano_plots/1_rds_files/volcano_legend.rds")
combined_legend <- plot_grid(pca_legend, volcano_legend)

# Loading plots from rds files
pca_all <- readRDS("../3_PCA/1_rds_files/PCA_plot.rds")
volcano_plot <- readRDS("../6_volcano_plots/1_rds_files/volcano_plot.rds")
box_plot_highlighted <- readRDS(
  "../7_highlighted_features/1_rds_files/box_plot_highlighted.rds"
)
bubble_plot <- readRDS("../8_bubble_plots/1_rds_files/bubble_plot.rds")

blank <- ggplot() + geom_blank() + theme_void()
bubble_plot <- bubble_plot + theme(legend.position = "none")

bubble_plot_no_legend <- plot_grid(
  bubble_plot, blank,
  ncol = 1,
  rel_heights = c(0.9, 0.1)
)

# Combining plots
row1_no_legend <- plot_grid(
  pca_all, volcano_plot,
  rel_widths = c(0.2, 0.8),
  labels = c("A", "B"),
  label_size = 8,
  label_fontface = "bold"
)

row1 <- plot_grid(
  row1_no_legend, combined_legend,
  nrow = 2,
  rel_heights = c(0.95, 0.05)
)

row2 <- plot_grid(
  bubble_plot_no_legend,
  box_plot_highlighted,
  rel_widths = c(0.5, 0.5),
  labels = c("C", "D"),
  label_size = 8,
  label_fontface = "bold"
)

final_figure <- plot_grid(row1, row2, nrow = 2, rel_heights = c(0.35, 0.65))

# Saving figure
ggsave(
  filename = "Figure5_LotusHordeum_rootexudates.pdf",
  plot = final_figure,
  width = 16,
  height = 18,
  units = "cm"
)
