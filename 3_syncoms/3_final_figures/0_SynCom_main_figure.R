# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "cowplot", "ggpubr", "ggtext")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# loading plots
p1 <- readRDS("1_rds_files/Plot_shootfw_boxplot.rds")
p2 <- readRDS("1_rds_files/Plot_PCoA_matched.rds")
p3 <- readRDS("1_rds_files/Plot_RA_sign_orders.rds")
p4 <- readRDS(
  "../1_Lotus/8_DA_maaslin2/1_rds_files/LotusSynCom_DA_noNodule.rds"
)
p5 <- readRDS(
  "../2_Hordeum/7_DA_maaslin2/1_rds_files/HordeumSynCom_DA_noNodule.rds"
)
lgd <- readRDS("../2_Hordeum/7_DA_maaslin2/1_rds_files/legend.rds")

# Empty space
p0 <- ggplot() +
  geom_blank() +
  theme_void() +
  annotate(
    "text",
    x = 0.5,
    y = 0.5,
    label = "",
    size = 6,
    fontface = "bold",
    hjust = 0.5,
    vjust = 0.5
  )

row1 <- ggarrange(
  p0, p1, p2,
  ncol = 3,
  widths = c(0.33, 0.2, 0.47),
  labels = c("A","B","C"),
  font.label = list(size = 8, color = "black", face = "bold"),
  align = "v",
  common.legend = TRUE,
  legend = "bottom"
)

row2 <- plot_grid(
  p3,
  labels = "D",
  label_size = 8,
  label_fontface = "bold"
)

row3 <- ggarrange(
  p4,
  p5,
  nrow = 1,
  labels = c("E", ""),
  font.label = list(size = 8),
  widths = c(0.6, 0.4),
  align = "h"
)

# Combining rows
final_plot <- plot_grid(
  row1,
  row2,
  row3,
  lgd,
  ncol = 1,
  rel_heights = c(0.3, 0.24, 0.4, 0.06)
)

# Saving plot
ggsave("Figure_6_SynCom.pdf", final_plot, width = 16, height = 20, units = "cm")
