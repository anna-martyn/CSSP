# Clean up.
options(warn=-1)
rm(list=ls())

# Load packages and import plots. -----------------------------------------------
pkg <- c("data.table", "ggplot2", "vegan", "ggtext", "ggpubr", "grid", 
         "cowplot")

for(pk in pkg){
  library(pk, character.only = T)
}

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load individal files and legends etc. for plot.
empty <- ggplot() + theme_void()
Boxplot_Lj <- readRDS("LotusCSSP_AskovSoils_shootfw.rds")
Boxplot_Hv <- readRDS("HordeumCSSP_AskovSoils_shootfw.rds")

CPCoA_PCoA_Lj <- readRDS("LotusCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.rds")
CPCoA_PCoA_Hv <- readRDS("HordeumCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.rds")
CPCoA_Hv <- readRDS("HordeumCSSP_AskovSoils_cpcoa_with_legend.rds")

legend <- ggpubr::get_legend(CPCoA_Hv)

Heatmap <- readRDS("Combined_heatmap.rds")

row1 <- plot_grid(
  empty, Boxplot_Lj, Boxplot_Hv,
  ncol = 3,
  rel_widths = c(0.33, 0.34, 0.32),
  labels = c("A","B",""),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row2 <- plot_grid(
  CPCoA_PCoA_Lj, CPCoA_PCoA_Hv, legend,
  ncol = 3,
  rel_widths = c(0.45, 0.45, 0.1),
  labels = c("C","",""),
  label_size = 15,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

row3 <- Heatmap

final_plot <- plot_grid(
  row1, row2, row3,
  ncol = 1,
  labels = c("","","D"),
  label_size = 15,
  label_fontface = "bold",
  rel_heights = c(0.2, 0.2, 0.6)
)

ggsave("Figure3_Askov_CSSP.pdf", final_plot, 
       width = 21, height = 29.7, units = "cm")
