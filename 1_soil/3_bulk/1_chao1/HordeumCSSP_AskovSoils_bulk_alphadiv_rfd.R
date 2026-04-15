# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "dplyr", "multcompView")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data
alpha <- read.table(
  file = "../../2_rarefication_chao1/2_Hordeum/2_chao1/HordeumCSSP_AskovSoils_chao1.txt",
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

# Combining alpha diversity (chao1) and metadata
index <- cbind(alpha[, 1], design[match(row.names(alpha), row.names(design)), ] )
colnames(index)[1] <- "Chao1"

# Keeping only soil data
index_bulk <- subset(index, Genotype == "Soil")

# Colours
colors <- data.frame(
  group = c("NPK", "PK", "UF"),
  color = c("#6F944F", "#B2563C", "#3C7D82")
)
soils <- c("NPK", "PK", "UF")
index_bulk$Soil <- factor(index_bulk$Soil, levels = soils)
colors <- colors[match(soils, colors$group), ]

# ANOVA and plot --------------------------------------------------------------
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
chao_summary <- index_bulk %>%
  group_by(Soil) %>%
  summarise(
    Mean = mean(Chao1),
    Max = max(Chao1), 
    Min = min(Chao1), 
    Median = median(Chao1), 
    Std = sd(Chao1)
  )

# ANOVA and Tukey HSD
ano <- aov(Chao1 ~ Soil, data = index_bulk)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$Type=rownames(Tukey.labels)
  Tukey.labels <- Tukey.labels[order(Tukey.labels$Type) , ]
  return(Tukey.labels)
}
labels <- generate_label_df(pairwise , "Soil")

# Label positions above the max value by group
label_df <- chao_summary %>%
  mutate(Letters = labels$Letters, y_position = Max + 0.075 * (max(Max) - min(Min))) %>%
  select(Soil, Max, Letters, y_position)

# Overall ANOVA p-value
anova_p <- summary(ano)[[1]][["Pr(>F)"]][1]

# Boxplot
box_plot <- ggplot(index_bulk, aes(x = Soil, y = Chao1, fill = Soil)) +
  geom_boxplot(
    alpha = 0.7,
    position = position_dodge(width = 0.7),
    outlier.color = NA,
    width = 0.3,
    linewidth = 0.4
  ) +
  geom_jitter(
    position = position_jitter(width = 0, height = 0.17),
    size = 0.75,
    alpha = 1
  ) +
  scale_fill_manual(values = as.character(colors$color)) +
  labs(x = "", y = "Chao1 index") +
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

# Saving plot

ggsave(
  filename = "2_figures/HordeumCSSP_AskovSoils_bulk_chao1_rfd.pdf",
  plot = box_plot,
  width = 3,
  height = 5,
  units = "cm"
)
saveRDS(
  object = box_plot,
  file = "1_rds_files/HordeumCSSP_AskovSoils_bulk_chao1_rfd.rds"
)

# Save ANOVA and Tukey HSD output file
write.csv(
  x = labels,
  file = "3_tables/HordeumCSSP_AskovSoils_bulk_chao1_ANOVA_TukeyHSD_rfd.csv"
)

