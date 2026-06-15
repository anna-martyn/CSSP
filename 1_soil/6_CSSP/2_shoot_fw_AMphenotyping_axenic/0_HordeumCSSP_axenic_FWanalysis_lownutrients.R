# Setup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages 
pkg <- c("data.table", "ggplot2", "multcompView")

for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data
weight <- fread("../../1_data/2_Hordeum/HordeumCSSP_axenic_fw_input_lownutrients.txt")

# One symrk biological replicate was missing, hence we will remove it from the dataframe.
weight <- weight[!is.na(Shoot_weight)]

# Setting genotype colours
colors <- data.frame(
  group = c("WT", "symrk","ccamk","nsp1", "nsp2"), 
  color = c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD")
)

# Setting factor levels
weight[,":="(
  Genotype = factor(Genotype, levels = colors$group),
  AM_inoculum = factor(AM_inoculum, levels = c("no", "yes"), labels = c("no AMF", "AMF"))
)]

# Plots -----------------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text.x = element_text(
    size = 6,
    angle = 45,
    vjust = 1,
    hjust = 1,
    colour = "black"
  ),
  axis.text.y = element_text(size = 6, colour = "black"),
  legend.background = element_blank(),
  legend.key = element_blank(),
  text = element_text(family = "sans")
)

# Hypothesis tests using ANOVA and TukeyHSD
shoot_labels_list <- as.list(rep(NA, length(unique(weight$AM_inoculum))))

for(i in seq_along(unique(weight$AM_inoculum))){
  
  current_inoculum <- unique(weight$AM_inoculum)[i]
  
  weight_sub <- weight[AM_inoculum == current_inoculum]
  
  ano <- aov(Shoot_weight ~ Genotype, data = weight_sub)
  pairwise <- TukeyHSD(ano)
  
  letters <- multcompLetters(pairwise$Genotype[, "p adj"])
  
  dt <- data.table(
    Genotype = names(letters$Letters),
    Letter = letters$Letters,
    AM_inoculum = current_inoculum
  )
  
  shoot_labels_list[[i]] <- dt
}

shoot_labels_dt <- rbindlist(shoot_labels_list)

shoot_ypos <- weight[,
  .(y_pos = max(Shoot_weight)),
  by = .(AM_inoculum, Genotype)
]

shoot_labels_dt <- merge(shoot_labels_dt, shoot_ypos)
shoot_labels_dt[, y_pos := y_pos + 0.06]

# Box plot
## with and without AMF
shoot_fw_plot <- ggplot(
  weight,
  aes(x = Genotype, y = Shoot_weight, fill = Genotype)
) +
  geom_text(
    data = shoot_labels_dt,
    aes(x = Genotype, y = y_pos, label = Letter),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
  geom_boxplot(
    width = 0.5,
    position = position_dodge(width = 0.9),
    linewidth = 0.2,
    outlier.size = 0.5
  ) +
  geom_jitter(position = position_jitter(width = 0), size = 1, alpha = 0.3) +
  scale_fill_manual(values = as.character(colors$color)) +
  facet_wrap(~AM_inoculum, nrow = 1) +
  main_theme +
  ylab("Shoot fresh weight/plant [g]") +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Low nutrient levels") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 6, face = "bold"),
    legend.title = element_text(size = 6, colour = "black"),
    strip.text = element_text(size = 6, colour = "black", face = "bold"),
    legend.text = element_text(size = 6, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(
      size = 6,
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

shoot_fw_plot 

ggsave(
  filename = "./2_figures/HordeumCSSP_axenic_shoot_fw_lownutrients_all.pdf",
  plot = shoot_fw_plot,
  width = 6,
  height = 6,
  units = "cm"
)
saveRDS(
  object = shoot_fw_plot,
  file = "./1_rds_files/HordeumCSSP_axenic_shoot_fw_lownutrients_all.rds"
)

## without AMF only
weight_noAMF <- weight[AM_inoculum == "no AMF"]

### Recompute ANOVA for subset
ano_noAMF <- aov(Shoot_weight ~ Genotype, data = weight_noAMF)
pairwise_noAMF <- TukeyHSD(ano_noAMF)

letters_noAMF <- multcompLetters(pairwise_noAMF$Genotype[, "p adj"])

labels_noAMF <- data.table(
  Genotype = names(letters_noAMF$Letters),
  Letter = letters_noAMF$Letters
)

### Define y positions
ypos_noAMF <- weight_noAMF[,
  .(y_pos = max(Shoot_weight, na.rm = TRUE)), by = Genotype
]

labels_noAMF <- merge(labels_noAMF, ypos_noAMF, by = "Genotype")
labels_noAMF[, y_pos := y_pos + 0.06]

### Make plot
shoot_fw_plot_noAMF <- ggplot(
  weight_noAMF,
  aes(x = Genotype, y = Shoot_weight, fill = Genotype)
) +
  geom_text(
    data = labels_noAMF,
    aes(x = Genotype, y = y_pos, label = Letter),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
  geom_boxplot(
    width = 0.5,
    position = position_dodge(width = 0.9),
    linewidth = 0.2,
    outlier.size = 0.5
  ) +
  geom_jitter(position = position_jitter(width = 0), size = 1, alpha = 0.3) +
  scale_fill_manual(values = as.character(colors$color)) +
  main_theme +
  ylab("Shoot fresh weight/plant [g]") +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Low nutrient levels") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 6, face = "bold"),
    legend.title = element_text(size = 6, colour = "black"),
    strip.text = element_text(size = 6, colour = "black", face = "bold"),
    legend.text = element_text(size = 6, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(
      size = 6,
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

shoot_fw_plot_noAMF

ggsave(
  filename = "./2_figures/HordeumCSSP_axenic_shoot_fw_lownutrients.pdf",
  plot = shoot_fw_plot_noAMF,
  width = 6,
  height = 6,
  units = "cm"
)
saveRDS(
  object = shoot_fw_plot_noAMF,
  file = "./1_rds_files/HordeumCSSP_axenic_shoot_fw_lownutrients.rds"
)
