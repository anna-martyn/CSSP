# Seup ------------------------------------------------------------------------
# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
pkg <- c("ggplot2", "dplyr", "multcompView")
for(pk in pkg){
  library(pk, character.only = T)
}

# Load chao1 and metadata
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

# Combine alpha diversity (chao1) and metadata
index <- cbind(alpha[, 1], design[match(row.names(alpha), row.names(design)), ] )
colnames(index)[1] <- "value"

# Keep only soil data
index_bulk <- subset(index, Genotype == "Soil")

# Colours
colors <- data.frame(
  group = c("NPK", "PK", "UF"),
  color = c("#6F944F", "#B2563C", "#3C7D82")
)
soils <- c("NPK", "PK", "UF")
index_bulk$Soil <- factor(index_bulk$Soil, levels = soils)
colors <- colors[match(soils, colors$group), ]

# ANOVA and figure ------------------------------------------------------------
# Set main theme
main_theme <- theme(
  panel.background=element_blank(),
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
    Mean = mean(value),
    Max = max(value), 
    Min = min(value), 
    Median = median(value), 
    Std = sd(value)
  )

# ANOVA and Tukey HSD
ano <- aov(value ~ Soil, data = index_bulk)
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

# label positions above the max value by group
label_df <- chao_summary %>%
  mutate(Letters = labels$Letters, y_position = Max + 0.075 * (max(Max) - min(Min))) %>%
  select(Soil, Max, Letters, y_position)

# Overall ANOVA p-value
anova_p <- summary(ano)[[1]][["Pr(>F)"]][1]

# Boxplot with significance letters and ANOVA p-value as title
box_plot <- ggplot(index_bulk, aes(x = Soil, y = value, fill = Soil)) +
  geom_boxplot(
    alpha = 0.7, position=position_dodge(width = 0.7), outlier.color = NA, width = 0.3
  ) +
  geom_jitter(position = position_jitter(width = 0, height = 0.17), size=1, alpha = 1) +
  scale_fill_manual(values = as.character(colors$color)) +
  labs(x = "", y = "Chao1 index") + 
  geom_text(
    data = label_df, 
    aes(x = Soil, y = y_position, label = Letters), 
    inherit.aes = FALSE, 
    size = 6/.pt
  ) +
  main_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  NULL

box_plot

# Save the plot
ggsave("HordeumCSSP_AskovSoils_bulk_chao1_rfd.pdf", box_plot, width = 3, height = 5, units = "cm")
saveRDS(box_plot, file = "HordeumCSSP_AskovSoils_bulk_chao1_rfd.rds")
saveRDS(box_plot, file = "../5_final_figure/HordeumCSSP_AskovSoils_bulk_chao1_rfd.rds")

# Save ANOVA and Tukey HSD output file
write.csv(labels, file = "HordeumCSSP_AskovSoils_bulk_chao1_ANOVA_TukeyHSD_rfd.csv")

