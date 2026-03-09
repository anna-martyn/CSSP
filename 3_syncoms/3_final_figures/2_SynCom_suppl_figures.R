# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
pkg <- c(
  "ggplot2", "patchwork", "cowplot", "magick", "grid", "forcats", "ggtext", "ggpubr"
)
for(pk in pkg){
  library(pk, character.only = T)
}

# Read all plot files (rds files) into variables.
plot_files <- list.files(pattern = "\\.rds$")
for (f in plot_files) {
  plot_name <- tools::file_path_sans_ext(f)
  assign(plot_name, readRDS(f))
}

# Load the plots wanted for the supplementary figures.
p1 <- LotusSC_chao1_matchedASVsonly_rfd_combined

p2 <- HordeumSC_chao1_matchedASVsonly_rfd

p3 <- LotusSC_cpcoa_matched_ASVs
p3 <- p3 +
  theme(legend.position = "none")

p4 <- HordeumSC_cpcoa_matched_ASVs
p4 <- p4 +
  theme(legend.position = "none")

legend_cpcoa <- ggpubr::get_legend(
  p3 +
    guides(color = guide_legend(title = "Genotype"),
           shape = guide_legend(title = "Compartment")) +
    theme(legend.position = "right")
) 

p5 <- LotusSC_order_RA_stackedbp

p6 <- HordeumSC_order_RA_stackedbp

p7 <- LotusSC_nodule_cts

p8 <- LotusSynCom_symbionts_RA

p9 <- LotusSynCom_DA_with_Nodule

# Now define the rows for the first supplementary figure of this section and then make the final plot.
row1 <- plot_grid(
  p1, p2,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 15,
  rel_widths = c(2.7/5, 2.3/5)
)

row2 <- plot_grid(
  p7, p3, p4, legend_cpcoa,
  ncol = 4,
  labels = c("C", "D", "E", ""),
  label_size = 15,
  rel_widths = c(1/3, 0.8/3, 0.8/3, 0.4/3)
)

row3 <- plot_grid(
  p5, p6,
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

final_plot1

# Now do the same for the other supplementary figure.
row4 <- plot_grid(
  p8,
  ncol = 1,
  labels = c("A"),
  label_size = 15,
  rel_widths = c(1/1)
)

row5 <- plot_grid(
  p9,
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

final_plot2

# Save supplementary figures as PDF files.
ggsave("Suppl_Figure6_SynCom.pdf", final_plot1, width=21, height=29.7, unit="cm")
ggsave("Suppl_Figure7_SynCom.pdf", final_plot2, width=21, height=29.7, unit="cm")
