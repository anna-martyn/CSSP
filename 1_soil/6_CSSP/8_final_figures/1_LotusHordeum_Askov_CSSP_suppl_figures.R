# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "cowplot", "ggpubr")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Mutant names in italics
genotype_labels <- c(
  "WT"    = "WT",
  "symrk" = "<i>symrk</i>",
  "ccamk" = "<i>ccamk</i>",
  "nsp1"  = "<i>nsp1</i>",
  "nsp2"  = "<i>nsp2</i>"
)

# Loading plots for supplementary figures
p1 <- readRDS("../3_nod_cts/1_rds_files/LotusCSSP_AskovSoils_nod_cts.rds")
p2 <- readRDS(
  "../2_shoot_fw_AMphenotyping_axenic/1_rds_files/HordeumCSSP_axenic_shoot_fw_highnutrients.rds"
)
p3 <- readRDS(
  "../2_shoot_fw_AMphenotyping_axenic/1_rds_files/HordeumCSSP_axenic_shoot_fw_lownutrients.rds"
)
p4 <- readRDS(
  "../2_shoot_fw_AMphenotyping_axenic/1_rds_files/HordeumCSSP_axenic_AMphenotyping_boxplot.rds"
)
p5 <- readRDS(
  "../4_chao1/1_Lotus/1_rds_files/LotusCSSP_AskovSoils_chao1_rfd.rds"
)
p6 <- readRDS(
  "../4_chao1/2_Hordeum/1_rds_files/HordeumCSSP_AskovSoils_chao1_rfd.rds"
)
p7 <- readRDS("../5_cpcoa_pcoa/1_rds_files/LotusCSSP_AskovSoils_cpcoa_all.rds")
p8 <- readRDS(
  "../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_cpcoa_all.rds"
)
p9 <- readRDS("../5_cpcoa_pcoa/1_rds_files/LotusCSSP_AskovSoils_pcoa_all.rds")
p10 <- readRDS("../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_pcoa_all.rds")
p11 <- readRDS(
  "../6_stackedbp_heatmap/1_rds_files/LotusHordeum_Askov_stackedbp_top20_meanRA.rds"
)
p12 <- readRDS(
  "../6_stackedbp_heatmap/1_rds_files/LotusHordeum_Askov_orders_heatmap.rds"
)

# Adjusting plot layouts
p5 <- p5 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

p6 <- p6 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

p11 <- p11 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

p12 <- p12 +
  scale_x_discrete(labels = genotype_labels) +
  theme(axis.text.x = ggtext::element_markdown())

# Extracting common legend for CPCoA plots
legend_plot <- readRDS(
  "../5_cpcoa_pcoa/1_rds_files/HordeumCSSP_AskovSoils_cpcoa_with_legend.rds"
)
cpcoa_legend <- ggpubr::get_legend(legend_plot)
pcoa_legend <- ggpubr::get_legend(
  legend_plot +
    guides(
      color = guide_legend(title = "Genotype"),
      shape = "none",
      fill = "none"
    )
)

# Individual rows for supplementary figure 2
row1 <- plot_grid(
  p1, p2, p3,
  ncol = 3,
  labels = c("A", "B", ""),
  label_size = 8,
  rel_widths = c(2/4, 1/4,1/4)
)

row2 <- plot_grid(
  p4,
  ncol = 1,
  labels =c("C"),
  label_size = 8,
  rel_widths = c("1/1")
)

# Combined supplementary figure 2
final_plot <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(1/2, 1/2)
)

# Individual rows for supplementary figure 3
row1 <- plot_grid(
  p5, p6,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 8,
  rel_widths = c(1/2, 1/2)
)

row2 <- plot_grid(
  p7, p8, cpcoa_legend,
  ncol = 3,
  labels = c("C", "D", ""),
  label_size = 8,
  rel_widths = c(5/11, 5/11, 1/11)
)

row3 <- plot_grid(
  p9, p10, pcoa_legend,
  ncol = 3,
  labels = c("E", "F"),
  label_size = 8,
  rel_widths = c(5/11, 5/11, 1/11)
)

# Combined supplementary figure 3
final_plot2 <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1.5/4, 0.75/4, 1.75/4)
)

# Individual rows for supplementary figure 4
row4 <- plot_grid(
  p11,
  ncol = 1,
  labels = c("A"),
  label_size = 8,
  rel_widths = c(1/1)
)

row5 <- plot_grid(
  p12,
  ncol = 1,
  labels = c("B"),
  label_size = 8,
  rel_widths = c(1/1)
)

# Combined supplementary figure 4
final_plot3 <- plot_grid(
  row4,
  row5,
  ncol = 1,
  rel_heights = c(0.5, 0.5)
)

# Saving figures
ggsave(
  "Suppl_Figure2_Askov_CSSP.pdf",
  final_plot,
  width = 18,
  height = 14,
  unit = "cm"
)
ggsave(
  "Suppl_Figure3_Askov_CSSP.pdf",
  final_plot2,
  width = 18,
  height = 22,
  unit = "cm"
)
ggsave(
  "Suppl_Figure4_Askov_CSSP.pdf",
  final_plot3,
  width = 18,
  height = 20,
  unit = "cm"
)
