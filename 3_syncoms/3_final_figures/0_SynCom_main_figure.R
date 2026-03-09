# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
pkg <- c(
  "ggplot2", "patchwork", "cowplot", "magick", "grid", "ggtext", "ggpubr"
)
for(pk in pkg){
  library(pk, character.only = T)
}

# Read all plot files (rds files) into variables
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)
  assign(plot_name, readRDS(f))
}

# Load plots.
p0 <- ggplot() + 
  geom_blank() + 
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "", size = 6, fontface="bold", hjust=0.5, vjust=0.5)
# p1 <- LotusSC_shootfw_LjSC_only
p1 <- Plot_shootfw_boxplot
# p1_aligned <- p1 + labs(title = " ") +
#   theme(plot.title = element_text(size = 8, color = NA),
#         axis.text.x = element_text(angle = 90, hjust=1))
# p2 <- LotusSynCom_pcoa_matchedASVs_byCompartment
p2 <- Plot_PCoA_matched
# p3 <- LotusSC_order_RA_stackedbp
# p4 <- Lotus_order_RA_sign_orders
# p4 <- Lotus_order_RA_sign_orders_asterisks
p3 <- Plot_RA_sign_orders
p4 <- LotusSynCom_DA_noNodule
p5 <- HordeumSynCom_DA_noNodule
lgd <- readRDS("../2_Hordeum/7_DA_maaslin2/legend.rds")
# p7 <- Hordeum_order_RA_sign_orders_asterisks

# Adjust legends of p3 and p4.
# p3 <- p3 +
#   theme(legend.position = "bottom") +
#   guides(
#     fill = guide_legend(
#       ncol = 3, nrow = 3,
#       title.position = "top",
#       title.hjust = 0
#     )
#   )

# p4 <- p4 +
#   theme(legend.position = "bottom") +
#   guides(
#     fill = guide_legend(
#       ncol = 5, nrow = 1,
#       title.position = "top",
#       title.hjust = 0
#     )
#   )


# Combine in one plot.
## Define individual rows.
# row1 <- plot_grid(
#   p0, p1_aligned, p2,
#   ncol = 3,
#   rel_widths = c(0.33, 0.2, 0.47),
#   labels = c("A","B","C"),
#   label_size = 15,
#   label_fontface = "bold",
#   label_x = 0, label_y = 1
# )

# p2 <- p2 + guides(color = "none")

row1 <- ggarrange(
  p0, p1, p2,
  ncol = 3,
  widths = c(0.33, 0.2, 0.47),
  labels = c("A","B","C"),
  font.label = list(size = 12, color = "black", face = "bold"),
  align = "v",
  common.legend = T,
  legend = "bottom"
)

row2 <- plot_grid(
  p3,
  labels = "D",
  label_size = 12,
  label_fontface = "bold"
)

# row2 <- plot_grid(
#   # p3, p4,
#   p4, p7,
#   ncol = 2,
#   rel_widths = c(0.5, 0.5),
#   labels = c("D",""),
#   label_size = 15,
#   label_fontface = "bold",
#   label_x = 0, label_y = 1
# )

# p_RA_Hv <- readRDS("../../2_Hordeum/8_final_figures/p_RA_Hv.rds")
# p_tax_clean_Hv <- readRDS("../../2_Hordeum/8_final_figures/p_tax_clean_Hv.rds")
# p_bubble_Hv <- readRDS("../../2_Hordeum/8_final_figures/p_bubble_Hv.rds")
# 
# p_RA_Lj <- readRDS("p_RA_Lj.rds")
# p_tax_clean_Lj <- readRDS("p_tax_clean_Lj.rds")
# p_bubble_Lj <- readRDS("p_bubble_Lj.rds")
# 
# row3 <- p_RA_Lj + p_RA_Hv + p_tax_clean_Lj + p_tax_clean_Hv +
#   p_bubble_Lj + p_bubble_Hv + 
#   plot_layout(nrow = 3, ncol = 2, 
#               heights = c(0.54, 0.04, 0.42),
#               widths = c(0.57, 0.43))

# p5 <- p5 + guides(fill = "none")
# p6 <- p6 + guides(fill = "none")
# 
# row3 <- plot_grid(
#   p5, p6,
#   ncol = 2,
#   rel_widths = c(0.57, 0.43),
#   labels = c("F", ""),
#   label_size = 15,
#   label_fontface = "bold",
#   label_x = 0, label_y = 1
# )

# p5 <- p5 + theme(plot.margin = margin(c(0.5, 0, 0.5, 0), unit = "lines"))
# p6 <- p6 + theme(plot.margin = margin(c(0.5, 0, 0.5, 0), unit = "lines"))

row3 <- ggarrange(p4, p5,
                  nrow = 1,
                  labels = c("E", ""),
                  font.label = list(size = 12),
                  widths = c(0.57, 0.43),
                  # common.legend = T,
                  # legend = "bottom",
                  align = "h")

## Combine all rows.
final_plot <- plot_grid(
  row1,
  row2,
  row3,
  lgd,
  ncol = 1,
  # labels = c("", "", ""),
  # label_size = 15,
  rel_heights = c(0.3, 0.24, 0.4, 0.06)
)

final_plot

# Save the combined plot as PDF
ggsave("Figure_6_SynCom.pdf", final_plot, width = 18, height = 22, 
       units = "cm")

# pdf("LotusSC_figure.pdf", width = 21/2.54, height = 29.7/2.54)
# final_plot
# grid.draw(linesGrob(x = unit(c(0.15, 0.15), "npc"),
#                     y = unit(c(0, 0.4), "npc")))
# dev.off()


# pdf("Bulk_figure_combined_final_3.pdf", width = 21/2.54, height = 25/2.54)
# combined
# grid.draw(linesGrob(x = unit(c(0.28, 0.34), "npc"),
#                     y = unit(c(0.25, 0.34), "npc")))
# dev.off()