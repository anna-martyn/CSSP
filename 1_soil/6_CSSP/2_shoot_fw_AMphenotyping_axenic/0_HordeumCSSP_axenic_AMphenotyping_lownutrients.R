# Setup ------------------------------------------------------------------------

# Cleaning up.
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages -------------------------------------------------------------
pkg <- c("ggplot2", "multcompView", "dplyr", "tidyr")

for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data -----------------------------------------------------------------
am <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_axenic_AMphenotyping_input.txt",
  header = TRUE,
  sep = "\t"
)

# Setting genotype colours -----------------------------------------------------
colors <- c(
  "WT"    = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1"  = "#C0E4EF",
  "nsp2"  = "#6EA6CD"
)

# Setting factor levels for genotypes ------------------------------------------
am$Genotype <- factor(
  am$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

# Setting structure order ------------------------------------------------------
structures <- c(
  "spores",
  "extraradical_hyphae",
  "intraradical_hyphae",
  "intracellular_hyphae",
  "vesicles",
  "arbuscules"
)

structure_labels <- c(
  spores = "Spores",
  extraradical_hyphae = "Extraradical hyphae",
  intraradical_hyphae = "Intraradical hyphae",
  intracellular_hyphae = "Intracellular hyphae",
  vesicles = "Vesicles",
  arbuscules = "Arbuscules"
)

# Converting dataframe to long format ------------------------------------------
am_long <- am %>%
  pivot_longer(
    cols = all_of(structures),
    names_to = "Structure",
    values_to = "Presence"
  ) %>%
  mutate(
    Structure = factor(Structure, levels = structures)
  )

# Percentage of presence of each AM structure per root piece -------------------
am_rootpiece <- am_long %>%
  group_by(Genotype, Bio_rep, Root_piece, Structure) %>%
  summarise(
    Percentage = mean(Presence) * 100,
    .groups = "drop"
  )

# Mean percentage per biological replicate -------------------------------------
am_percent <- am_rootpiece %>%
  group_by(Genotype, Bio_rep, Structure) %>%
  summarise(
    Percentage = mean(Percentage),
    .groups = "drop"
  ) %>%
  mutate(
    Structure = factor(Structure, levels = structures)
  )

# ANOVA + Tukey ---------------------------------------------------------------
letters_list <- list()

for(struct in structures){
  
  temp <- filter(am_percent, Structure == struct)
  
  aov_res <- aov(Percentage ~ Genotype, data = temp)
  tukey <- TukeyHSD(aov_res)
  
  tukey_pvals <- tukey$Genotype[, "p adj"]
  names(tukey_pvals) <- rownames(tukey$Genotype)
  
  letters <- multcompLetters(tukey_pvals)$Letters
  
  letters_df <- data.frame(
    Genotype = names(letters),
    label = letters,
    Structure = struct,
    stringsAsFactors = FALSE
  )
  
  p_value <- summary(aov_res)[[1]][["Pr(>F)"]][1]
  
  letters_df$anova_label <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ "ns"
  )
  
  letters_list[[struct]] <- letters_df
}

letters_df <- bind_rows(letters_list)

# Labels ----------------------------------------------------------------------
label_positions <- am_percent %>%
  group_by(Structure, Genotype) %>%
  summarise(
    y_pos = max(Percentage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = c("Structure", "Genotype"))

anova_labels <- am_percent %>%
  group_by(Structure) %>%
  summarise(
    y_pos = 48,
    .groups = "drop"
  ) %>%
  left_join(
    letters_df %>% distinct(Structure, anova_label),
    by = "Structure"
  )

label_positions$Structure <- factor(label_positions$Structure, levels = structures)
anova_labels$Structure <- factor(anova_labels$Structure, levels = structures)

# Theme -----------------------------------------------------------------------
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
    hjust = 1,
    vjust = 1,
    colour = "black"
  ),
  axis.text.y = element_text(size = 6),
  legend.background = element_blank(),
  legend.key = element_blank(),
  text = element_text(family = "sans")
)

# Plot -------------------------------------------------------------------------
am_plot <- ggplot(
  am_percent %>%
    mutate(Structure = factor(Structure, levels = structures)),
  aes(x = Genotype, y = Percentage, fill = Genotype)) +
  geom_boxplot(width = 0.3, outlier.color = NA, alpha = 0.7, linewidth = 0.3) +
  # geom_jitter(width = 0.1,height = 0, size = 1.5, alpha = 0.8) +
  geom_jitter(position = position_jitter(width = 0), size = 1, alpha = 0.3) +
  geom_text(
    data = label_positions,
    aes(x = Genotype, y = pmin(y_pos + 5, 90), label = label),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
  geom_text(
    data = anova_labels,
    aes(x = 3, y = y_pos, label = anova_label),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
  
  facet_wrap(
    ~Structure,
    nrow = 1,
    scales = "fixed",
    drop = FALSE,
    labeller = labeller(Structure = structure_labels)
  ) +
  
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.02))
  ) +
  main_theme +
  ylab("Occurence of AM structures (%)") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  scale_x_discrete(
    labels = c(
      WT = "WT",
      symrk = expression(italic("symrk")),
      ccamk = expression(italic("ccamk")),
      nsp1 = expression(italic("nsp1")),
      nsp2 = expression(italic("nsp2"))
    )
  )

am_plot

# Save plot --------------------------------------------------------------------
write.csv(
  label_positions,
  file = "3_tables/HordeumCSSP_axenic_AMphenotyping_ANOVA.csv",
  row.names = FALSE
)

ggsave(
  filename = "2_figures/HordeumCSSP_axenic_AMphenotyping_boxplot.pdf",
  plot = am_plot,
  width = 20,
  height = 8,
  units = "cm"
)

saveRDS(
  am_plot,
  file = "1_rds_files/HordeumCSSP_axenic_AMphenotyping_boxplot.rds"
)


