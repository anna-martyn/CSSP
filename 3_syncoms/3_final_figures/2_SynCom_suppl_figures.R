# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "ggtext", "ggpubr", "cowplot")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading plots
chao1_plot_Lj <- readRDS(
  "../1_Lotus/5_chao1/1_rds_files/LotusSC_chao1_matchedASVsonly_rfd_combined.rds"
)
chao1_plot_Hv <- readRDS(
  "../2_Hordeum/4_chao1/1_rds_files/HordeumSC_chao1_matchedASVsonly_rfd.rds"
)
cpcoa_plot_Lj <- readRDS(
  "../1_Lotus/6_cpcoa_pcoa/1_rds_files/LotusSC_cpcoa_matched_ASVs.rds"
)
cpcoa_plot_Hv <- readRDS(
  "../2_Hordeum/5_cpcoa_pcoa/1_rds_files/HordeumSC_cpcoa_matched_ASVs.rds"
)
bar_plot_Lj <- readRDS(
  "../1_Lotus/7_stackedbp_barplots/1_rds_files/LotusSC_order_RA_stackedbp.rds"
)
bar_plot_Hv <- readRDS(
  "../2_Hordeum/6_stackedbp_barplots/1_rds_files/HordeumSC_order_RA_stackedbp.rds"
)
box_plot_nodules <- readRDS(
  "../1_Lotus/3_nodule_cts/1_rds_files/LotusSC_nodule_cts.rds"
)
box_plot_symbionts <- readRDS(
  "../1_Lotus/9_symbionts/1_rds_files/LotusSynCom_symbionts_RA.rds"
)
bubble_plot <- readRDS(
  "../1_Lotus/8_DA_maaslin2/1_rds_files/LotusSynCom_DA_withNodule.rds"
)

# Extracting legends for for CPCoA plots
legend_cpcoa <- ggpubr::get_legend(
  cpcoa_plot_Lj +
    guides(
      color = guide_legend(title = "Genotype"),
      shape = guide_legend(title = "Compartment")
    ) +
    theme(legend.position = "right")
)

# Removing legends
cpcoa_plot_Lj <- cpcoa_plot_Lj + theme(legend.position = "none")
cpcoa_plot_Hv <- cpcoa_plot_Hv + theme(legend.position = "none")

# Defining rows for supplementary figure
row1 <- plot_grid(
  chao1_plot_Lj, chao1_plot_Hv,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 15,
  rel_widths = c(2.7/5, 2.3/5)
)

row2 <- plot_grid(
  box_plot_nodules, cpcoa_plot_Lj, cpcoa_plot_Hv, legend_cpcoa,
  ncol = 4,
  labels = c("C", "D", "E", ""),
  label_size = 15,
  rel_widths = c(1/3, 0.8/3, 0.8/3, 0.4/3)
)

row3 <- plot_grid(
  bar_plot_Lj, bar_plot_Hv,
  ncol = 2,
  labels = c("F", "G"),
  label_size = 15,
  rel_widths = c(1/2, 1/2)
)

final_plot1 <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1/3, 1/3, 1/3)
)

# Defining rows for another supplementary figure
row4 <- plot_grid(
  box_plot_symbionts,
  ncol = 1,
  labels = c("A"),
  label_size = 15,
  rel_widths = c(1/1)
)

row5 <- plot_grid(
  bubble_plot,
  ncol = 1,
  labels = c("B"),
  label_size = 15,
  rel_widths = c(1/1)
)

final_plot2 <- plot_grid(
  row4,
  row5,
  ncol = 1,
  rel_heights = c(1/5, 4/5)
)

# Saving supplementary figures
ggsave(
  filename = "Suppl_Figure6_SynCom.pdf",
  plot = final_plot1,
  width = 21,
  height = 29.7,
  unit = "cm"
)

ggsave(
  filename = "Suppl_Figure7_SynCom.pdf",
  plot = final_plot2,
  width = 21,
  height = 29.7,
  unit = "cm"
)
