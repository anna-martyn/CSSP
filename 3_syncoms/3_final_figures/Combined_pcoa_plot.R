# Set working directory and load packages --------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pkg <- c(
  "data.table", "ggplot2", "ggtext"
)
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Settings ---------------------------------------------------------------------
# Genotypes colors
colors <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

legend_labels <- c(
  "WT" = "WT",
  "symrk" = "*symrk*",
  "ccamk" = "*ccamk*",
  "nsp1" = "*nsp1*",
  "nsp2" = "*nsp2*"
)

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text.x = element_text(size = 6, colour = "black"),
  axis.text.y = element_text(size = 6, colour = "black"),
  legend.background = element_blank(),
  legend.key = element_blank(),
  text = element_text(family = "sans")
)

# Load data --------------------------------------------------------------------
pcoa_point_Lj <- fread(
  "../1_Lotus/6_cpcoa_pcoa/3_tables/LotusSC_PCoA_points_segments_matched_ASVs.csv",
  drop = 6:8
)

pcoa_point_Hv <- fread(
  "../2_Hordeum/5_cpcoa_pcoa/3_tables/HordeumSC_PCoA_points_segments_matched_ASVs.csv",
  drop = 6:8
)

text_dt_Lj <- fread(
  "../1_Lotus/6_cpcoa_pcoa/3_tables/LotusSC_PCoA_text_matched_ASVs.csv"
)
text_dt_Hv <- fread(
  "../2_Hordeum/5_cpcoa_pcoa/3_tables/HordeumSC_PCoA_text_matched_ASVs.csv"
)

pcoa_point <- rbind(pcoa_point_Lj, pcoa_point_Hv)
text_dt <- rbind(text_dt_Lj, text_dt_Hv)

pcoa_point[, ":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Genotype = factor(Genotype, levels = names(colors))
)]

text_dt[, ":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  text = gsub("-", "\n", text)
)]

text_dt[, text := gsub("52", "52.0", text)]

# Visualisation ----------------------------------------------------------------
cpcoa_plot <- ggplot(pcoa_point, aes(x = PCo1, y = PCo2, color = Genotype)) +
  geom_point(size = 1.5, alpha = 0.7) +
  facet_grid(Host ~ Compartment, switch = "y") +
  geom_segment(
    data = pcoa_point,
    aes(x = PCo1, y = PCo2, xend = seg_x, yend = seg_y, color = Genotype),
    alpha = 0.5
  ) +
  geom_label(
    data = text_dt,
    aes(x = -0.175, y = -0.25, label = text),
    colour = "black",
    fill = "grey",
    alpha = 0.2,
    size = 6 / .pt
  ) +
  scale_color_manual(values = colors, labels = legend_labels) +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  labs(
    x = "PCo 1",
    y = "PCo 2"
  ) +
  main_theme +
  theme(
    plot.title = element_text(face = "bold", size = 6, hjust = 0),
    legend.text = element_markdown(size = 6, color = "black"),
    strip.text = element_text(size = 6, colour = "black", face = "bold"),
    strip.placement = "outside",
    axis.title.x = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, colour = "black"),
    legend.key.size = unit(0.25, "cm")
  ) +
  NULL

ggsave(
  filename = "2_temp_figures/Plot_PCoA_matched.pdf",
  plot = cpcoa_plot,
  width = 5,
  height = 13,
  unit = "cm"
)
saveRDS(object = cpcoa_plot, file = "1_rds_files/Plot_PCoA_matched.rds")
