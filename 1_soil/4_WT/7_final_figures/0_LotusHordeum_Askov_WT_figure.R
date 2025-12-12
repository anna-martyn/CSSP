# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library(ggplot2)
library(ggpubr)
library(patchwork)
library(cowplot)
library(magick)
library(grid)

# Read all plot files (rds files) into variables.
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

p2 <- LotusHordeum_Askov_WT_shootfw_boxplots

p3 <- Lotus_Askov_WT_nodule_cts
p3_aligned <- p3 + labs(title = " ") +
  theme(plot.title = element_text(size = 8, color = NA))

p4 <- Lotus_Askov_WT_cpcoa +
  theme(legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin())+
  guides(shape = guide_legend(title.position = "top", title.hjust=0.5))
p4 <- p4 + theme(legend.spacing.y = unit(0.5, "lines"))

p5 <- Hordeum_Askov_WT_cpcoa + 
  theme(legend.position = "bottom",
        legend.box="vertical", 
        legend.margin=margin())+
  guides(shape = guide_legend(title.position = "top", title.hjust=0.5))
p5 <- p5 + theme(
  legend.spacing.y = unit(0.5, "lines"),
  plot.margin = margin(t = 0.5, b = 1, l = 0.5, r = 0.5, unit = "lines")
)

p6 <- LotusHordeum_Askov_WT_orders_heatmap

# Assemble the individual rows for the final figure.
row1 <- plot_grid(
  p1, p2, p3_aligned,
  ncol = 3,
  rel_widths = c(1.2, 1.2, 0.6),
  labels = c("A","B","C"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

col1 <- ggarrange(
  p4 + theme(legend.spacing.y = unit(0, "pt")),
  p5 + theme(legend.spacing.y = unit(0, "pt")),
  ncol = 1,
  labels = c("D", ""),
  font.label = list(size = 15, face = "bold"),
  common.legend = T,
  legend = "bottom"
) + 
  NULL

bottom_rows <- plot_grid(
  col1, p6,
  ncol = 2,
  labels = c("", "E"),
  label_size = 15,
  rel_widths = c(1/3, 2/3)
)

# Combine all rows in the final figure.
final_plot <- plot_grid(
  row1,
  bottom_rows,
  ncol = 1,
  rel_heights = c(0.3, 0.7)
)

final_plot

# Save the final figure as a PDF file.
ggsave("Figure2_Askov_WT.pdf", final_plot, 
       width = 21, height = 23, units = "cm")

