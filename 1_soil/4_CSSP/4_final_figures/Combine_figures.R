# Load packages and import plots -----------------------------------------------
pkg <- c("data.table", "ggplot2", "vegan", "ggtext", "ggpubr", "grid", 
         "cowplot")

for(pk in pkg){
  library(pk, character.only = T)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

empty <- ggplot() + theme_void()
Boxplot_Lj <- readRDS("../1_shoot_fw/Lotus_CSSP_fw_BOX_Lj.rds")
Boxplot_Hv <- readRDS("../1_shoot_fw/Lotus_CSSP_fw_BOX_Hv.rds")

CPCoA_PCoA_Lj <- readRDS("../2_cpcoa_pcoa/CPCoA_PCoA_plot_Lj.rds")
CPCoA_PCoA_Hv <- readRDS("../2_cpcoa_pcoa/CPCoA_PCoA_plot_Hv.rds")
CPCoA_Hv <- readRDS("../2_cpcoa_pcoa/CPCoA_plot_with_legend_Hv.rds")

legend <- ggpubr::get_legend(CPCoA_Hv)

Heatmap <- readRDS("../3_heatmaps/Combined_heatmap.rds")

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

ggsave("CSSP_fig.pdf", final_plot, 
       width = 21, height = 29.7, units = "cm")
