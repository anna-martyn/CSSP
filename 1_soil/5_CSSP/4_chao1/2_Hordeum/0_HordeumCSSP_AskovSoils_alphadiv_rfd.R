# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
pkg <- c("dplyr", "ggplot2", "multcompView")
for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
alpha <- read.table(
  "../../../2_rarefication_chao1/2_Hordeum/2_chao1/HordeumCSSP_AskovSoils_chao1.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
design <- read.table(
  "../../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Combining alpha diversity (chao1) and metadata
index <- cbind(alpha[, 1], design[match(row.names(alpha), row.names(design)), ])
colnames(index)[1] <- "Chao1"

# Removing bulk soil
index_filtered <- index %>%
  filter(Compartment != "Bulk")

# Seting factor levels
index_filtered$Soil <- factor(
  index_filtered$Soil,
  levels = c("NPK", "PK", "UF")
)
index_filtered$Compartment <- factor(
  index_filtered$Compartment,
  levels = c("Rhizosphere", "Root")
)
index_filtered$Genotype <- factor(
  index_filtered$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

# Setting genotype colours
colors <- data.frame(
  group = c("WT", "symrk", "ccamk", "nsp1", "nsp2"),
  color = c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD")
)

# Vector for mutants in italic
genotype_labels <- c(
  "WT" = "WT",
  "symrk" = "italic(symrk)",
  "ccamk" = "italic(ccamk)",
  "nsp1" = "italic(nsp1)",
  "nsp2" = "italic(nsp2)"
)

# Hypothesis testing ----------------------------------------------------------

# Defining function for significance letters
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$Type = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$Type), ]
  return(Tukey.labels)
}

# Empty dataframe for results
label_df <- data.frame()

# Using function for each soil-compartment combination
for(s in levels(index_filtered$Soil)){
  for(c in levels(index_filtered$Compartment)){
    sub_df <- subset(index_filtered, Soil == s & Compartment == c)

    # Skip if fewer than 2 genotypes present (this is the case for the nodule compartment)
    if (length(unique(sub_df$Genotype)) < 2) {
      next
    }

    # Summary dataframe
    summary_df <- sub_df %>%
      group_by(Genotype) %>%
      summarise(
        Mean = mean(Chao1),
        Max = max(Chao1),
        Min = min(Chao1),
        Median = median(Chao1),
        Std = sd(Chao1),
        .groups = "drop"
      )

    # ANOVA and TukeyHSD
    ano <- aov(Chao1 ~ Genotype, data = sub_df)
    pairwise <- TukeyHSD(ano)
    letters <- generate_label_df(pairwise, "Genotype")
    letters <- letters[as.character(summary_df$Genotype), ]

    # Adding letters and defining y-positions
    y_offset <- 0.1 * (max(summary_df$Max) - min(summary_df$Min))

    tmp <- summary_df %>%
      mutate(
        Letters = letters$Letters,
        y_position = Max + y_offset,
        Soil = s,
        Compartment = c
      )

    label_df <- rbind(label_df, tmp)
  }
}

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

box_plots <- ggplot(
  index_filtered,
  aes(x = Genotype, y = Chao1, fill = Genotype)
) +
  geom_boxplot(
    alpha = 0.7,
    position = position_dodge(width = 0.7),
    outlier.color = NA,
    width = 0.3
  ) +
  geom_jitter(
    position = position_jitter(width = 0, height = 0.17),
    size = 1,
    alpha = 1
  ) +
  scale_fill_manual(values = as.character(colors$color)) +
  scale_x_discrete(labels = function(x) parse(text = genotype_labels[x])) +
  labs(x = "", y = "Chao1 index", title = "Hordeum") +
  geom_text(
    data = label_df,
    aes(x = Genotype, y = y_position, label = Letters),
    inherit.aes = FALSE,
    size = 8 / .pt
  ) +
  facet_grid(Compartment ~ Soil, scales = "free_y", drop = FALSE) +
  main_theme +
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Saving plot
ggsave(
  filename = "2_figures/HordeumCSSP_AskovSoils_chao1_rfd.pdf",
  plot = box_plots,
  width = 12,
  height = 10,
  unit = "cm"
)
saveRDS(
  object = box_plots,
  file = "1_rds_files/HordeumCSSP_AskovSoils_chao1_rfd.rds"
)

# Saving ANOVA Tukey HSD results
write.csv(
  x = label_df,
  file = "3_tables/HordeumCSSP_AskovSoils_chao1_ANOVA_TukeyHSD_rfd.csv"
)

