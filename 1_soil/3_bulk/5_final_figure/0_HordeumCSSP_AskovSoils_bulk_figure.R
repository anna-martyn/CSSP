# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library(ggplot2)
library(cowplot)
library(patchwork)
library(magick)
library(grid)

# Read all RDS plots.
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)
  assign(plot_name, readRDS(f))
}

# Select the plots wanted for display and modify if necessary.
p1 <- placeholder <- ggplot() +
  geom_blank() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "", size = 6, hjust = 0.5) # workflow added manually later.

p2 <- HordeumCSSP_AskovSoils_bulk_chao1_rfd

p3 <- Hordeum_bulk_PCoA

p4 <- Hordeum_bulk_order_top20_RA_mean_stackedbp
p4_aligned <- p4 + labs(title = " ") + theme(plot.title = element_text(size = 6, color = NA))
p4_nolegend <- p4_aligned + theme(legend.position="none")

p5 <- Hordeum_barplot_bulk_top20_RA_sign
p5_aligned <- p5 + labs(title = " ") + theme(plot.title = element_text(size = 6, color = NA))

p6 <- Hordeum_bulk_Venn_ASVs
p6_aligned <- ggdraw(p6) + 
  theme(plot.margin = margin(t = 20, r = 5, b = 5, l = 5))

p7_nolegend <- HordeumCSSP_bulk_ASV_overlap_piecharts_nolegend
legend_p7 <- get_legend(
  p4_aligned + theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.25, 'cm')
  )
)

# Assemble the rows for the final figure.
row1 <- plot_grid(
  p1, p2, p3,
  ncol = 3,
  rel_widths = c(0.55, 0.15, 0.3),
  labels = c("A","B","C"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  axis = "b"
)

row2 <- plot_grid(
  p4_nolegend, p5_aligned,
  ncol = 2,
  rel_widths = c(0.2, 0.8),
  labels = c("D","E"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  axis = "b"
)

row3 <- plot_grid(
  p6_aligned, p7_nolegend, legend_p7,
  ncol = 3,
  rel_widths = c(0.325,0.5,0.175),
  labels = c("F","",""),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)


row4 <- ggplot() + theme_void()

# Combine all rows for the final figure.
combined <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(0.3, 0.35, 0.35)
)

combined

# Save the final figure and finally add a box around the piecharts for better visualisation.
pdf("Figure1_Askov_bulk.pdf", width = 21/2.54, height = 25/2.54)
combined
grid.draw(linesGrob(x = unit(c(0.28, 0.34), "npc"),
                    y = unit(c(0.25, 0.34), "npc")))
grid.draw(linesGrob(x = unit(c(0.28, 0.34), "npc"),
                    y = unit(c(0.12, 0.02), "npc")))
grid.draw(linesGrob(x = unit(c(0.34, 0.81), "npc"),
                    y = unit(c(0.02, 0.02), "npc")))
grid.draw(linesGrob(x = unit(c(0.34, 0.81), "npc"),
                    y = unit(c(0.34, 0.34), "npc")))
grid.draw(linesGrob(x = unit(c(0.81, 0.81), "npc"),
                    y = unit(c(0.02, 0.34), "npc")))
grid.draw(linesGrob(x = unit(c(0.34, 0.34), "npc"),
                    y = unit(c(0.02, 0.34), "npc")))
dev.off()
