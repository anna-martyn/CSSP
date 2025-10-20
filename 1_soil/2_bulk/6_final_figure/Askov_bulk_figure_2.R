# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
library(ggplot2)
library(cowplot)
library(patchwork)
library(magick)
library(grid)

# Read RDS plots.
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)
  assign(plot_name, readRDS(f))
}

# Assign plots to variables and modify if necessary.
p1 <- placeholder <- ggplot() +
  geom_blank() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "", size = 6, hjust = 0.5)

p2 <- Askov_barley_bulk_chao1_rfd_6603
p3 <- Barley_bulk_PCoA
p3_aligned <- p3 + labs(title = " ") + theme(plot.title = element_text(size = 20, color = NA))
p4 <- Barley_bulk_order_top20_RA_mean_stackedbp
p4_aligned <- p4 + labs(title = " ") + theme(plot.title = element_text(size = 20, color = NA))
p4_nolegend <- p4_aligned + theme(legend.position="none")
p5 <- Barley_barplot_bulk_top20_RA_sign
p5_aligned <- p5 + labs(title = " ") + theme(plot.title = element_text(size = 20, color = NA))
p6 <- Venn_ASVs_bulk
p6_aligned <- ggdraw(p6) + 
  theme(plot.margin = margin(t = 20, r = 5, b = 5, l = 5))
p7 <- Barley_bulk_ASV_overlap_piecharts
p7_nolegend <- Barley_bulk_ASV_overlap_piecharts_nolegend

legend_p7 <- get_legend(
  p4_aligned + theme(
    legend.position = "right",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
)

# Assemble rows.
row1 <- plot_grid(
  p1, p2, p3_aligned,
  ncol = 3,
  rel_widths = c(1.5, 0.8, 1),
  labels = c("A","B","C"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  align = "v", axis = "b"
)

row2 <- plot_grid(
  p4_nolegend, p5_aligned,
  ncol = 2,
  rel_widths = c(0.6,2),
  labels = c("D","E"),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  align = "v", axis = "b"
)

row3 <- plot_grid(
  p6_aligned, p7_nolegend, legend_p7,
  ncol = 3,
  rel_widths = c(1,2.3,0.7),
  labels = c("F","G",""),
  label_size = 30,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)


row4 <- ggplot() + theme_void()

# Combine all rows.
combined <- plot_grid(
  row1,
  row2,
  row3,
  # row4,
  ncol = 1,
  # rel_heights = c(1.5, 2, 2, 1.5),
  rel_heights = c(1.5, 2, 2)
)

combined

# Save plot.
ggsave("Bulk_figure_combined_final_3.pdf", combined, width=21, height=29.7)

