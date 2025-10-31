# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
library(ggplot2)
library(ggpubr)
library(patchwork)
library(cowplot)
library(magick)
library(grid)

# Read all plot files (rds files) into variables
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)  # remove .rds extension
  assign(plot_name, readRDS(f))
}

# Assign plots to variables and modify if necessary.
p1 <- placeholder <- ggplot() +
  geom_blank() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "", size = 6, hjust = 0.5)

p2 <- Askov_WT_shootfw_boxplots  

p3 <- LotusWT_pink_nod
p3_aligned <- p3 + labs(title = " ") +
  theme(plot.title = element_text(size = 8, color = NA))

# p4 <- cpcoa_Lotus_WT + theme(legend.position = "none")
p4 <- cpcoa_Lotus_WT
# p4 <- cpcoa_Lotus_WT + 
#   theme(legend.position = "bottom",
#         legend.box="vertical", 
#         legend.margin=margin())+
#   guides(shape = guide_legend(title.position = "top", title.hjust=0.5))

p5 <- Soil_WT_stackedbp_meanRA

p6 <- cpcoa_Barley_WT + 
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.margin=margin())+
  guides(shape = guide_legend(title.position = "top", title.hjust=0.5))

p7 <- Soil_WT_heatmap_orders

# Assemble rows.
row1 <- plot_grid(
  p1, p2, p3_aligned,
  ncol = 3,
  rel_widths = c(1.2, 1.2, 0.6),
  labels = c("A","B","C"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row2 <- ggarrange(
  p4, p6,
  ncol = 2,
  labels = c("D", ""),
  font.label = list(size = 15, face = "bold"),
  common.legend = T,
  legend = "right"
) + theme(
  plot.margin = margin(r = 0.5, l = 0.5, t = 0.5, b = 0.5, unit = "lines")
)

# row2 <- plot_grid(
#   p4, p5,
#   ncol = 2,
#   rel_widths = c(1,2),
#   labels = c("D","E"),
#   label_size = 15,
#   label_fontface = "bold",
#   label_x = 0, label_y = 1
# )

# col1 <- ggarrange(
#   p4, p6,
#   ncol = 1,
#   labels = c("D", ""),
#   font.label = list(size = 15, face = "bold"),
#   common.legend = T,
#   legend = "bottom"
# )
# 
# col2 <- plot_grid(
#   p5, p7,
#   ncol = 1,
#   labels = c("E", "F"),
#   label_size = 15,
#   label_fontface = "bold",
#   label_x = 0, label_y = 1
# )

# bottom_rows <- plot_grid(
#   col1, col2,
#   ncol = 2,
#   rel_widths = c(1,2.1)
# )

# bottom_rows <- plot_grid(
#   col1, p5,
#   ncol = 2,
#   rel_widths = c(1,2.1)
# )

# row3 <- plot_grid(
#   p6, p7,
#   ncol = 2,
#   rel_widths = c(1,2),
#   labels = c("","F"),
#   label_size = 15,
#   label_fontface = "bold",
#   label_x = 0, label_y = 1
# )

# Combine all rows.
# final_plot <- plot_grid(
#   row1,
#   row2,
#   row3,
#   ncol = 1,
#   rel_heights = c(1.6,2,2)  # adjust relative heights if needed
# )

# final_plot <- plot_grid(
#   row1,
#   bottom_rows,
#   ncol = 1,
#   rel_heights = c(0.28, 0.72)  # adjust relative heights if needed
# )

final_plot <- plot_grid(
  row1,
  row2,
  p5,
  ncol = 1,
  labels = c("", "", "E"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  rel_heights = c(0.25, 0.35, 0.4)  # adjust relative heights if needed
)

final_plot

# Save the combined plot as PDF
ggsave("Soil_WT_figure_2.pdf", final_plot, 
       width = 21, height = 23, units = "cm")

