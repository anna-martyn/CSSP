# Seup ------------------------------------------------------------------------
# Clean up
options(warn = -1)
rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
pkg <- c("ggplot2", "dplyr", "multcompView")
for(pk in pkg){
  library(pk, character.only = T)
}

# Load data
alpha <- read.table(
  "../../../2_rarefication_chao1/1_Lotus/2_chao1/LotusCSSP_AskovSoils_chao1.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
design <- read.table(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Combining alpha diversity (chao1) and metadata
index <- cbind(alpha[,1], design[match(row.names(alpha), row.names(design)), ])
colnames(index)[1] <- "Chao1"

# Keep only WT samples
index <- index %>% filter(Genotype == "WT")

# Setting the soil and compartment factor levels
soils <- c("NPK", "PK", "UF")
index$Soil <- factor(index$Soil, levels = soils)
index$Compartment <- factor(index$Compartment, levels = c("Rhizosphere", "Root", "Nodules"))

# Colours
colors <- data.frame(group = soils, color = c("#6F944F", "#B2563C", "#3C7D82"))
colors <- colors[match(soils, colors$group), ]

# Hypothesis testing ----------------------------------------------------------

# Function generating significance letters
generate_label_df <- function(pairwise, variable){
  Tukey_levels <- pairwise[[variable]][,4]
  Tukey_labels <- data.frame(multcompLetters(Tukey_levels)['Letters'])
  Tukey_labels$Type <- rownames(Tukey_labels)
  Tukey_labels <- Tukey_labels[order(Tukey_labels$Type), ]
  return(Tukey_labels)
}

# Empty dataframe for significance
label_df <- data.frame()

# Apply above function by compartment
for(comp in levels(index$Compartment)){
  
  sub_df <- subset(index, Compartment == comp)
  
  # Summary dataframe for all outputs
  summary_df <- sub_df %>%
    group_by(Soil) %>%
    summarise(
      Mean = mean(Chao1),
      Max = max(Chao1),
      Min = min(Chao1),
      Median = median(Chao1),
      Std = sd(Chao1),
      .groups = "drop"
    )
  
  # ANOVA and Tukey
  ano <- aov(Chao1 ~ Soil, data = sub_df)
  pairwise <- TukeyHSD(ano)
  letters <- generate_label_df(pairwise, "Soil")
  
  # Letters and the y-positions
  # y_offset <- 0.075 * (max(summary_df$Max) - min(summary_df$Min))
  y_offset <- 50
  
  tmp <- summary_df %>%
    mutate(
      Letters = letters$Letters,
      y_position = Max + y_offset,
      Compartment = comp
    )
  
  label_df <- rbind(label_df, tmp)
}

# Plot ------------------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 8, color = "black"),
  legend.text = element_text(size = 8, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 8),
  legend.position = "none",
  strip.text = element_text(size = 8, color = "black"),
  legend.background = element_blank(),
  plot.title = element_text(size = 8, hjust = 1)
)

box_plot <- ggplot(index, aes(x=Soil, y = Chao1, fill = Soil)) +
  geom_boxplot(
    alpha = 0.7,
    position = position_dodge(width = 0.7),
    outlier.color = NA,
    width = 0.3
  ) +
  geom_jitter(position = position_jitter(width = 0, height = 0.17), size = 1, alpha = 1) +
  scale_fill_manual(values = as.character(colors$color)) +
  labs(x = "", y = "Chao1 index", title = "Lotus") +
  geom_text(
    data = label_df,
    mapping = aes(x = Soil, y = y_position, label = Letters),
    inherit.aes = FALSE,
    size = 8/.pt
  ) +
  facet_wrap(
    vars(factor(Compartment, levels = c("Rhizosphere", "Root", "Nodules"))),
    scales = "fixed"
  ) +
  ylim(0, 600) +
  main_theme +
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0))

box_plot

# Save the plot.
ggsave("LotusCSSP_AskovSoils_WT_chao1_rfd.pdf", box_plot, width=7, height=5, unit="cm")
saveRDS(box_plot, file = "LotusCSSP_AskovSoils_WT_chao1_rfd.rds")
saveRDS(box_plot, file = "../../7_final_figures/LotusCSSP_AskovSoils_WT_chao1_rfd.rds")

# Save ANOVA Tukey HSD output file.
write.csv(label_df, file =  "LotusCSSP_AskovSoils_WT_chao1_ANOVA_TukeyHSD_rfd.csv")
