# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
pkg <- c("ggplot2", "cowplot", "grid")
for(pk in pkg){
  library(pk, character.only = T)
}

# Placeholder for workflow to be added in image editor
p1 <- placeholder <- ggplot() +
  geom_blank() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "", size = 6, hjust = 0.5)

# Reading plots from RDS files
p2 <- readRDS("HordeumCSSP_AskovSoils_bulk_chao1_rfd.rds")
p3 <- readRDS("Hordeum_bulk_PCoA.rds")
p4 <- readRDS("Hordeum_bulk_order_top20_RA_mean_stackedbp.rds")
p5 <- readRDS("Hordeum_barplot_bulk_top20_RA_sign.rds")
p6 <- readRDS("Hordeum_bulk_Venn_ASVs.rds")
p7_nolegend <- readRDS("HordeumCSSP_bulk_ASV_overlap_piecharts_nolegend.rds")

# Modifying plots
p4_aligned <- p4 + 
  labs(title = " ") + 
  theme(plot.title = element_text(size = 6, color = NA))
p4_nolegend <- p4_aligned + 
  theme(legend.position = "none")

p5_aligned <- p5 + 
  labs(title = " ") + 
  theme(plot.title = element_text(size = 6, color = NA))

p6_aligned <- ggdraw(p6) + 
  theme(plot.margin = margin(t = 20, r = 5, b = 5, l = 5))

legend_p7 <- get_legend(
  p4_aligned + theme(
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.25, 'cm')
  )
)

# Assembling the rows for final figure
row1 <- plot_grid(
  p1, p2, p3,
  ncol = 3,
  rel_widths = c(0.55, 0.15, 0.3),
  labels = c("A","B","C"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  axis = "b"
)

row2 <- plot_grid(
  p4_nolegend, p5_aligned,
  ncol = 2,
  rel_widths = c(0.2, 0.8),
  labels = c("D","E"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  axis = "b"
)

row3 <- plot_grid(
  p6_aligned, p7_nolegend, legend_p7,
  ncol = 3,
  rel_widths = c(0.325,0.5,0.175),
  labels = c("F","",""),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)


row4 <- ggplot() + theme_void()

# Combine all rows
combined <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(0.3, 0.35, 0.35)
)

combined

# Saving final figure and add a box around the piecharts
pdf("Figure1_Askov_bulk.pdf", width = 18/2.54, height = 22/2.54)
  combined
  grid.draw(linesGrob(x = unit(c(0.28, 0.34), "npc"), y = unit(c(0.25, 0.34), "npc")))
  grid.draw(linesGrob(x = unit(c(0.28, 0.34), "npc"), y = unit(c(0.12, 0.02), "npc")))
  grid.draw(linesGrob(x = unit(c(0.34, 0.81), "npc"), y = unit(c(0.02, 0.02), "npc")))
  grid.draw(linesGrob(x = unit(c(0.34, 0.81), "npc"), y = unit(c(0.34, 0.34), "npc")))
  grid.draw(linesGrob(x = unit(c(0.81, 0.81), "npc"), y = unit(c(0.02, 0.34), "npc")))
  grid.draw(linesGrob(x = unit(c(0.34, 0.34), "npc"), y = unit(c(0.02, 0.34), "npc")))
dev.off()
