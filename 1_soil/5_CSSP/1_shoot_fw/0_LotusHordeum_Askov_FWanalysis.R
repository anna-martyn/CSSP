# Seup ------------------------------------------------------------------------
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
weight <- fread("../../1_data/LotusHordeum_AskovSoils_shootfw.txt")

# Setting genotype colours
colors <- data.frame(
  group = c("WT", "symrk","ccamk","nsp1", "nsp2"), 
  color = c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD")
)

# Setting factor levels for soil and genotype
weight[,":="(
  Soil_type = factor(Soil_type, levels = c("NPK", "PK", "UF")),
  Genotype = factor(Genotype, levels = colors$group)
)]

# Splitting data by plant species
Lotus_weight <- weight[Plant_species == "Lotus"]
Hordeum_weight <- weight[Plant_species == "Hordeum"]

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
opt <- expand.grid(
  Plant = unique(weight$Plant_species),
  Soil = unique(weight$Soil_type)
)

labels_list <- as.list(rep(NA, nrow(opt)))
for (i in 1:nrow(opt)) {
  current_plant <- opt$Plant[i]
  current_soil <- opt$Soil[i]
  weight_sub <- weight[
    Plant_species == current_plant & Soil_type == current_soil
  ]
  ano <- aov(Fresh_weight ~ Genotype, data = weight_sub)
  pairwise <- TukeyHSD(ano)

  letters <- multcompLetters(pairwise$Genotype[, "p adj"])
  dt <- data.table(
    Genotype = names(letters$Letters),
    Letter = letters$Letters,
    Plant_species = current_plant,
    Soil_type = current_soil
  )
  labels_list[[i]] <- dt
}

labels_dt <- rbindlist(labels_list)

shootfw_by_group <- weight[,
  .(y_pos = max(Fresh_weight)),
  list(Plant_species, Soil_type, Genotype)
]

labels_dt <- merge(labels_dt, shootfw_by_group)
labels_dt[,y_pos:=y_pos + ifelse(Plant_species == "Lotus", 0.01, 0.2)]

# Box plots
## Lotus
box_plot_lotus <- ggplot(
  Lotus_weight,
  aes(x = Genotype, y = Fresh_weight, fill = Genotype)
) +
  geom_text(
    data = labels_dt[Plant_species == "Lotus"],
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
  scale_fill_manual(values = as.character(colors$color)) +
  facet_wrap(~Soil_type, nrow = 1) +
  main_theme +
  ylab("Shoot fresh weight/plant [g]") +
  scale_y_continuous(limits = c(0, 0.15)) +
  ggtitle("Lotus") +
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

ggsave(
  filename = "2_figures/LotusCSSP_AskovSoils_shootfw.pdf",
  plot = box_plot_lotus,
  width = 6,
  height = 5,
  units = "cm"
)
saveRDS(
  object = box_plot_lotus,
  file = "1_rds_files/LotusCSSP_AskovSoils_shootfw.rds"
)

## Hordeum
box_plot_hordeum <- ggplot(
  Hordeum_weight,
  aes(x = Genotype, y = Fresh_weight, fill = Genotype)
) +
  geom_text(
    data = labels_dt[Plant_species == "Hordeum"],
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
  scale_fill_manual(values = as.character(colors$color)) +
  facet_wrap(~Soil_type, scales = "free_x", nrow = 1) +
  main_theme +
  ylab("Shoot fresh weight/plant [g]") +
  scale_y_continuous() +
  ggtitle("Hordeum") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 6, face = "bold"),
    legend.title = element_text(size = 6, colour = "black"),
    strip.text = element_text(size = 6, colour = "black", face = "bold"),
    legend.text = element_text(size = 6, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
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

ggsave(
  filename = "2_figures/HordeumCSSP_AskovSoils_shootfw.pdf",
  plot = box_plot_hordeum,
  width = 6,
  height = 5,
  units = "cm"
)

saveRDS(
  object = box_plot_hordeum,
  file = "1_rds_files/HordeumCSSP_AskovSoils_shootfw.rds"
)
