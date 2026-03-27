# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("data.table", "ggplot2")

for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
nod <- fread("../../1_data/1_Lotus/LotusCSSP_Askov_nodule_cts_all.txt")

# Setting genotype colours
colors <- data.frame(
  group = c("WT", "symrk","ccamk","nsp1", "nsp2"), 
  color = c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD")
)

# Set the factor levels for soils and genotypes.
nod[,":="(
  Soil_type = factor(Soil_type, levels = c("NPK", "PK", "UF")),
  Genotype = factor(Genotype, levels = colors$group)
)]

# Plot ------------------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text.x = element_text(
    size = 8,
    angle = 45,
    vjust = 1,
    hjust = 1,
    colour = "black"
  ),
  axis.text.y = element_text(size = 8, colour = "black"),
  legend.background = element_blank(),
  legend.key = element_blank(),
  text = element_text(family = "sans")
)

# Plot
box_plot <- ggplot(nod, aes(x = Genotype, y = pink, fill = Genotype)) +
  geom_boxplot(
    width = 0.5,
    position = position_dodge(width = 0.9),
    outlier.size = 0.5
  ) +
  scale_fill_manual(values = as.character(colors$color)) +
  facet_wrap(~Soil_type, nrow = 1) +
  main_theme +
  ylab("Pink nodule counts/plant") +
  scale_y_continuous(limits = c(0, 15)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, colour = "black"),
    strip.text = element_text(size = 8, colour = "black", face = "bold"),
    legend.text = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      vjust = 1,
      hjust = 1,
      colour = "black"
    ),
    legend.key.size = unit(0.5, "cm")
  ) +
  scale_x_discrete(
    labels = c(
      "WT" = "WT",
      "symrk" = expression(italic("symrk")),
      "ccamk" = expression(italic("ccamk")),
      "nsp1" = expression(italic("nsp1")),
      "nsp2" = expression(italic("nsp2"))
    )
  )

ggsave(
  filename = "2_figures/LotusCSSP_AskovSoils_nod_cts.pdf",
  plot = box_plot,
  width = 6,
  height = 5,
  units = "cm"
)

saveRDS(
  object = box_plot,
  file = "1_rds_files/LotusCSSP_AskovSoils_nod_cts.rds"
)
