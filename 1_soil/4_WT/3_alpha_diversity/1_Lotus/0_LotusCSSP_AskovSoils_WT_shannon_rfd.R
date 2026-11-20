# Setup ------------------------------------------------------------------------
# Alpha diversity (Shannon entropy) of WT samples by compartment and soil.
# Rarefaction and alpha diversity are both computed here in R.

# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "dplyr", "multcompView", "phyloseq", "vegan")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data
asv_table <- read.table(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
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

# Rarefaction ------------------------------------------------------------------
min(colSums(asv_table)) > 1000
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = TRUE))
set.seed(1673967505)
phs_rarefied <- rarefy_even_depth(phs, verbose = FALSE)
asv_table_rarefied <- as.matrix(otu_table(phs_rarefied))

# Alpha diversity --------------------------------------------------------------
# Shannon entropy (natural logarithm) per sample
shannon <- diversity(t(asv_table_rarefied), index = "shannon")

# Combining alpha diversity (Shannon) and metadata
index <- data.frame(
  Shannon = shannon,
  design[match(names(shannon), row.names(design)), ],
  check.names = FALSE
)

# Keeping only WT samples
index <- index %>% filter(Genotype == "WT")

# Setting the soil and compartment factor levels
soils <- c("NPK", "PK", "UF")
compartments <- c("Rhizosphere", "Root", "Nodules")
index$Soil <- factor(index$Soil, levels = soils)
index$Compartment <- factor(index$Compartment, levels = compartments)

# Colours
colors <- data.frame(group = soils, color = c("#6F944F", "#B2563C", "#3C7D82"))
colors <- colors[match(soils, colors$group), ]

# Hypothesis testing -----------------------------------------------------------

# Offset of the significance letters above the highest point in each group
y_offset <- 0.075 * diff(range(index$Shannon))

# Empty dataframes for significance letters and full Tukey HSD output
label_df <- data.frame()
tukey_df <- data.frame()

# ANOVA and Tukey HSD by compartment
for(comp in levels(index$Compartment)){

  sub_df <- subset(index, Compartment == comp)

  # Summary dataframe for all outputs
  summary_df <- sub_df %>%
    group_by(Soil) %>%
    summarise(
      n = n(),
      Mean = mean(Shannon),
      Max = max(Shannon),
      Min = min(Shannon),
      Median = median(Shannon),
      Std = sd(Shannon),
      .groups = "drop"
    )

  # ANOVA and Tukey
  ano <- aov(Shannon ~ Soil, data = sub_df)
  cat("\n==", comp, "==\n")
  print(anova(ano))
  pairwise <- TukeyHSD(ano)

  # Compact letter display, matched to soils by name
  letters_soil <- multcompLetters(pairwise[["Soil"]][, "p adj"])$Letters
  letters <- data.frame(
    Soil = factor(names(letters_soil), levels = soils),
    Letters = unname(letters_soil)
  )

  tmp <- summary_df %>%
    left_join(letters, by = "Soil") %>%
    mutate(
      y_position = Max + y_offset,
      Compartment = comp,
      ANOVA_p = anova(ano)[["Pr(>F)"]][1]
    )
  label_df <- rbind(label_df, tmp)

  tukey_df <- rbind(
    tukey_df,
    data.frame(
      Compartment = comp,
      Comparison = rownames(pairwise[["Soil"]]),
      pairwise[["Soil"]],
      check.names = FALSE,
      row.names = NULL
    )
  )
}
print(tukey_df)

label_df$Compartment <- factor(label_df$Compartment, levels = compartments)

# Plot -------------------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  legend.position = "none",
  strip.text = element_text(size = 6, color = "black", face = "bold"),
  legend.background = element_blank(),
  plot.title = element_text(size = 6, hjust = 1)
)

box_plot <- ggplot(index, aes(x = Soil, y = Shannon, fill = Soil)) +
  geom_boxplot(
    alpha = 0.7,
    position = position_dodge(width = 0.7),
    outlier.color = NA,
    linewidth = 0.3,
    width = 0.3
  ) +
  geom_jitter(
    position = position_jitter(width = 0.1, height = 0),
    size = 1,
    alpha = 1
  ) +
  scale_fill_manual(values = as.character(colors$color)) +
  labs(x = "", y = "Shannon index", title = "Lotus") +
  geom_text(
    data = label_df,
    mapping = aes(x = Soil, y = y_position, label = Letters),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
  facet_wrap(vars(Compartment), scales = "fixed") +
  main_theme +
  theme(plot.title = element_text(face = "bold", size = 6, hjust = 0))

# Saving -----------------------------------------------------------------------
ggsave(
  filename = "2_figures/LotusCSSP_AskovSoils_WT_shannon_rfd.pdf",
  plot = box_plot,
  width = 7,
  height = 5,
  units = "cm"
)
saveRDS(
  object = box_plot,
  file = "1_rds_files/LotusCSSP_AskovSoils_WT_shannon_rfd.rds"
)

# Saving ANOVA and Tukey HSD results
write.csv(
  x = label_df,
  file = "3_tables/LotusCSSP_AskovSoils_WT_shannon_ANOVA_TukeyHSD_rfd.csv",
  row.names = FALSE
)
write.csv(
  x = tukey_df,
  file = "3_tables/LotusCSSP_AskovSoils_WT_shannon_TukeyHSD_pairwise_rfd.csv",
  row.names = FALSE
)
