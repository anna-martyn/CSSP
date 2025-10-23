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
p1 <- LotusSC_shootfw_LjSC_only
p1_aligned <- p1 + labs(title = " ") +
  theme(plot.title = element_text(size = 8, color = NA),
        axis.text.x = element_text(angle=45, hjust=1))
p2 <- LotusSynCom_pcoa_matchedASVs_byCompartment
p3 <- LotusSC_order_RA_stackedbp
# p4 <- Lotus_order_RA_sign_orders
p4 <- Lotus_order_RA_sign_orders_asterisks
p5 <- LotusSynCom_DA

# Adjust legends of p3 and p4.
p3 <- p3 +
  theme(legend.position = "bottom") +
  guides(
    fill = guide_legend(
      ncol = 3, nrow = 3,
      title.position = "top",
      title.hjust = 0
    )
  )

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
row1 <- plot_grid(
  p0, p1_aligned, p2,
  ncol = 3,
  rel_widths = c(0.4,0.2,0.6),
  labels = c("A","B","C"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row2 <- plot_grid(
  p3, p4,
  ncol = 2,
  rel_widths = c(1,1.2),
  labels = c("D","E"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row3 <- plot_grid(
  p5,
  ncol = 1,
  rel_widths = c(1),
  labels = c("F"),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)


## Combine all rows.
final_plot <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(0.6,0.8,1.6)
)

final_plot

# Save the combined plot as PDF
ggsave("LotusSC_figure.pdf", final_plot, width=21, height=27, units = "cm")

