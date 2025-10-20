# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
library(ggplot2)
library(patchwork)
library(cowplot)
library(magick)
library(grid)

# Read workflow plot (PDF from BioRender) and convert to raster for cowplot
workflow_img <- magick::image_read_pdf("Bulk_workflow_2.pdf", density = 600)
workflow_grob <- grid::rasterGrob(as.raster(workflow_img), interpolate = TRUE)

# Read all plot files (rds files) into variables
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)  # remove .rds extension
  assign(plot_name, readRDS(f))
}

# Convert workflow rasterGrob to ggplot (panel A)
p1_gg <- ggplot() +
  annotation_custom(workflow_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void()

# Prepare the other plots (no internal tags!)
p2 <- Askov_barley_bulk_chao1_rfd_6603   # has a title
p3 <- Barley_bulk_PCoA                   # no title
p4 <- Barley_bulk_order_top20_RA_stackedbp
p5 <- Barley_barplot_bulk_top20_RA
p6 <- Venn_ASVs_bulk                     # no title
p7 <- Barley_bulk_ASV_overlap_piecharts  # has plot_annotation title

# Add invisible title to C so it matches B’s height
p3_aligned <- p3 + labs(title = " ") +
  theme(plot.title = element_text(size = 10, color = NA))

# Add invisible title to F so it matches G’s height
p6_aligned <- p6 + plot_annotation(
  title = " ",   # blank title to reserve space
  theme = theme(
    plot.title = element_text(size = 20, colour = NA)  # invisible text
  )
)

# Adjust D for legend
p4_adjusted <- p4 +
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.direction = "vertical",
    legend.box = "vertical"
  ) +
  guides(fill = guide_legend(ncol = 1))

p4_with_legend <- plot_grid(
  p4_adjusted + theme(legend.position = "none"),
  get_legend(p4_adjusted),
  ncol = 2,
  rel_widths = c(0.35, 0.3)
)

# === Assemble rows ===

# Row1: workflow (A) + row_BC (B,C)
row_BC <- plot_grid(
  p2, p3_aligned,
  ncol = 2,
  labels = c("B","C"),
  label_size = 20,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row1 <- plot_grid(
  p1_gg, row_BC,
  ncol = 2,
  rel_widths = c(1.7, 1.3),
  labels = c("A",""),
  label_size = 20,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

# Row2: D (with legend) + E
row2 <- plot_grid(
  p4_with_legend, p5,
  ncol = 2,
  rel_widths = c(1, 2),
  labels = c("D","E"),
  label_size = 20,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

# Row3: F + G (now aligned at top)
row3 <- plot_grid(
  p6_aligned, p7,
  ncol = 2,
  rel_widths = c(1, 3),
  labels = c("F","G"),
  label_size = 20,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

# Row4: spacer (optional)
row4 <- ggplot() + theme_void()

# === Combine all rows ===
combined <- plot_grid(
  row1,
  row2,
  row3,
  row4,
  ncol = 1,
  rel_heights = c(0.8, 1, 1, 0.7)
)

# Show result
combined

# Save the combined plot as PDF
ggsave("Bulk_figure_combined.pdf", combined, width=21, height=29.7)
