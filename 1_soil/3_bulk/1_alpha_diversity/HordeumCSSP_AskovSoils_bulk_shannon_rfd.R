# Setup ------------------------------------------------------------------------
# Alpha diversity (Shannon entropy) of the unplanted (bulk) soil samples.

# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "dplyr", "multcompView", "vegan")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data
# Rarefied ASV table, from
# 2_rarefication/2_Hordeum/1_rarefication/Hordeum_Askov_ASV10_4_rarefication.R
asv_table_rarefied <- read.table(
  file = "../../2_rarefication/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4_rfd.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

design <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Alpha diversity --------------------------------------------------------------
# Shannon entropy (natural logarithm) per sample
shannon <- diversity(t(asv_table_rarefied), index = "shannon")

index <- data.frame(
  Shannon = shannon,
  design[match(names(shannon), row.names(design)), ],
  check.names = FALSE
)

# Keeping only unplanted soil samples
index_bulk <- subset(index, Plant == "Soil")

# Colours
colors <- data.frame(
  group = c("NPK", "PK", "UF"),
  color = c("#6F944F", "#B2563C", "#3C7D82")
)
soils <- c("NPK", "PK", "UF")
index_bulk$Soil <- factor(index_bulk$Soil, levels = soils)
colors <- colors[match(soils, colors$group), ]

# ANOVA and plot ---------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  legend.position = "none",
  legend.background = element_blank(),
  plot.title = element_text(size = 6, hjust = 1)
)

# Summary statistics
shannon_summary <- index_bulk %>%
  group_by(Soil) %>%
  summarise(
    n = n(),
    Mean = mean(Shannon),
    Max = max(Shannon),
    Min = min(Shannon),
    Median = median(Shannon),
    Std = sd(Shannon)
  )

# ANOVA and Tukey HSD
ano <- aov(Shannon ~ Soil, data = index_bulk)
pairwise <- TukeyHSD(ano)

# Compact letter display, matched to soils by name
tukey_p <- pairwise[["Soil"]][, "p adj"]
letters_soil <- multcompLetters(tukey_p)$Letters
labels <- data.frame(
  Soil = factor(names(letters_soil), levels = soils),
  Letters = unname(letters_soil)
)

# Label positions above the max value by group
label_df <- shannon_summary %>%
  left_join(labels, by = "Soil") %>%
  mutate(y_position = Max + 0.075 * (max(Max) - min(Min))) %>%
  select(Soil, Max, Letters, y_position)

# Boxplot
box_plot <- ggplot(index_bulk, aes(x = Soil, y = Shannon, fill = Soil)) +
  geom_boxplot(
    alpha = 0.7,
    position = position_dodge(width = 0.7),
    outlier.color = NA,
    width = 0.3,
    linewidth = 0.4
  ) +
  geom_jitter(
    position = position_jitter(width = 0.1, height = 0),
    size = 0.75,
    alpha = 1
  ) +
  scale_fill_manual(values = as.character(colors$color)) +
  labs(x = "", y = "Shannon index") +
  geom_text(
    data = label_df,
    aes(x = Soil, y = y_position, label = Letters),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
  main_theme +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    margins = margin(
      t = unit(5.5, "pt"),
      r = unit(1, "pt"),
      b = unit(-20, "pt"),
      l = unit(0, "pt")
    )
  ) +
  NULL

# Saving -----------------------------------------------------------------------
ggsave(
  filename = "2_figures/HordeumCSSP_AskovSoils_bulk_shannon_rfd.pdf",
  plot = box_plot +
    theme(
      margins = margin(
        t = unit(5.5, "pt"),
        r = unit(5.5, "pt"),
        b = unit(-10, "pt"),
        l = unit(5.5, "pt")
      )
    ),
  width = 3,
  height = 5,
  units = "cm"
)
saveRDS(
  object = box_plot,
  file = "1_rds_files/HordeumCSSP_AskovSoils_bulk_shannon_rfd.rds"
)

# Save Shannon values and Tukey HSD output
write.csv(
  x = data.frame(Comparison = rownames(pairwise[["Soil"]]), pairwise[["Soil"]], check.names = FALSE),
  file = "3_tables/HordeumCSSP_AskovSoils_bulk_shannon_ANOVA_TukeyHSD_rfd.csv",
  row.names = FALSE
)
