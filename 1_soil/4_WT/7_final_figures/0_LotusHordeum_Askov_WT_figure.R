# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
pkg <- c("ggplot2", "ggpubr", "cowplot")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Blank placeholder
p1 <- placeholder <- ggplot() +
  geom_blank() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "", size = 6, hjust = 0.5)

# Loading plots from rds files
p2 <- readRDS("../1_shoot_fw/1_rds_files/LotusHordeum_Askov_WT_shootfw_boxplots.rds")
p3 <- readRDS("../2_nodule_cts/1_rds_files/Lotus_Askov_WT_nodule_cts.rds")
p4 <- readRDS("../4_cpcoa/1_rds_files/Lotus_Askov_WT_cpcoa.rds")
p5 <- readRDS("../4_cpcoa/1_rds_files/Hordeum_Askov_WT_cpcoa.rds")
p6 <- readRDS("../5_stacked_bp_heatmap/1_rds_files/LotusHordeum_Askov_WT_orders_heatmap.rds")

# Modifying plots
p3_aligned <- p3 + labs(title = " ") +
  theme(plot.title = element_text(size = 6, color = NA))

# p4 <- p4 +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.margin = margin(),
#     legend.spacing.y = unit(0, "pt")
#   )

# p5 <- p5 + 
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical", 
#     legend.margin = margin(),
#     legend.spacing.y = unit(0, "pt"),
#     plot.margin = margin(t = 0.5, b = 1, l = 0.5, r = 0.5, unit = "lines")
#   )

cpcoa_legend <- ggpubr::get_legend(p4)
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")

heatmap_legend <- ggpubr::get_legend(p6)
p6 <- p6 + theme(legend.position = "none")

legends <- ggarrange(cpcoa_legend, heatmap_legend)

# Assembling rows
row1 <- plot_grid(
  p1, p2, p3_aligned,
  ncol = 3,
  rel_widths = c(1.2, 1.2, 0.6),
  labels = c("A","B","C"),
  label_size = 8,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

col1 <- ggarrange(
  p4, p5,
  ncol = 1,
  labels = c("D", ""),
  font.label = list(size = 8, face = "bold")
)

bottom_rows <- plot_grid(
  col1, p6,
  ncol = 2,
  labels = c("", "E"),
  label_size = 8,
  rel_widths = c(1/3, 2/3)
)

# Combining all plots
final_plot <- plot_grid(
  row1, bottom_rows, legends,
  ncol = 1,
  rel_heights = c(0.3, 0.6, 0.1)
)

# Saving figure
ggsave(
  filename = "Figure2_Askov_WT.pdf",
  plot = final_plot,
  width = 11,
  height = 11,
  units = "cm"
)

