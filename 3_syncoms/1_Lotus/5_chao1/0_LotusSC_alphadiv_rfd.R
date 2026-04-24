# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
pkg <- c("data.table", "ggplot2", "patchwork", "multcompView")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data
alpha_filt <- fread(
  "0_chao1_qiime/matchedASVsonly/LotusSC_matchedASVsonly_chao1.txt"
)
setnames(alpha_filt, "V1", "SampleID")
design <- fread("../1_data/LotusSC_metadata.txt")

# Adding Chao1 index to metadata
design <- merge(design, alpha_filt, by = "SampleID")

# Genotypes colours
colors <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

# Hypothesis testing -----------------------------------------------------------
# Function to perform ANOVA, tukey HSD, and extract letters
anova_letters <- function(comp){
  anova_comp <- aov(chao1~Genotype, data = design[Compartment == comp])
  tk_comp <- TukeyHSD(anova_comp)
  letters <- multcompLetters(tk_comp$Genotype[,"p adj"])
  letters <- letters$Letters
  return(data.table(
    Compartment = comp,
    Genotype = names(letters),
    Letter = letters
  ))
}

# Applyign function and merging results
anova_root <- anova_letters("Root")
anova_rhizo <- anova_letters("Rhizosphere")
anova_res <- rbind(anova_root, anova_rhizo)

design <- merge(
  design,
  anova_res,
  by = c("Compartment", "Genotype"),
  all.x = TRUE
)

# Adding y-position for letters
design_y_pos <- design[,.(y_pos = max(chao1)),.(Compartment, Genotype)]
design <- merge(
  design,
  design_y_pos,
  c("Compartment", "Genotype"),
  all.x = TRUE
)

# Plots ------------------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text.y = element_text(size = 6, color = "black"),
  axis.text.x = element_text(hjust = 1, angle = 30, size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  legend.position = "none",
  strip.text = element_text(size = 6, color = "black", face = "bold"),
  legend.background = element_blank(),
  plot.title = element_text(size = 6, hjust = 1)
)

# Setting factor levels
design[, ":="(
  Genotype = factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  ),
  Compartment = factor(
    Compartment,
    levels = c("Rhizosphere", "Root", "Nodules")
  )
)]

# Boxplot
box_plot <- ggplot(design, aes(x = Genotype, y = chao1, fill = Genotype)) +
  geom_boxplot(width = 0.3, outlier.color = NA, alpha = 0.7, linewidth = 0.3) +
  geom_jitter(position = position_jitter(width = 0), size = 1, alpha = 0.5) +
  geom_text(
    data = design,
    aes(x = Genotype, y = y_pos * 1.05, label = Letter),
    inherit.aes = FALSE,
    size = 6 /.pt
  ) +
  scale_fill_manual(values = colors) +
  main_theme +
  labs(y = "Chao1 index", x = "") +
  scale_x_discrete(
    labels = c(
      "WT" = "WT",
      "symrk" = expression(italic("symrk")),
      "ccamk" = expression(italic("ccamk")),
      "nsp1" = expression(italic("nsp1")),
      "nsp2" = expression(italic("nsp2"))
    )
  ) +
  facet_wrap(~Compartment, scales = "free_x", space = "free_x") +
  # expand_limits(y = 0)+
  NULL

# Saving object
saveRDS(
  object = box_plot,
  file = "1_rds_files/LotusSC_chao1_matchedASVsonly_rfd_combined.rds"
)

ggsave(
  filename = "2_figures/LotusSC_chao1_matchedASVsonly_rfd_combined.pdf",
  plot = box_plot,
  width = 15,
  height = 6,
  units = "cm"
)


