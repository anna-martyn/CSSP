# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "multcompView", "dplyr")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data
weight <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_axenic_fw_input.txt",
  header = TRUE,
  sep = "\t"
)

# Setting colours and factor levels for genotypes
colors <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

weight$Genotype <- factor(
  weight$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

# Hypothesis tests ------------------------------------------------------------
# ANOVA and TukeyHSD
aov_res <- aov(Shoot_fw ~ Genotype, data = weight)
tukey <- TukeyHSD(aov_res)
tukey_pvals <- tukey$Genotype[, "p adj"]
names(tukey_pvals) <- rownames(tukey$Genotype)
letters <- multcompLetters(tukey_pvals)$Letters
letters_df <- data.frame(
  Genotype = names(letters),
  label = letters,
  stringsAsFactors = FALSE
)

# Global test for genotype effect
p_value <- summary(aov_res)[[1]][["Pr(>F)"]][1]
asterisk <- case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  TRUE ~ NA_character_
)

# Summary
weight_summary <- weight %>%
  group_by(Genotype) %>%
  summarise(y_pos = max(Shoot_fw, na.rm = TRUE), .groups = "drop") %>%
  left_join(letters_df, by = "Genotype")

# Plot ------------------------------------------------------------------------
# Dummy points for plot limits
upper_limit <- max(weight$Shoot_fw, na.rm=TRUE) * 1.2

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

# Plot
box_plot <- ggplot(weight, aes(x = Genotype, y = Shoot_fw, fill = Genotype)) +
  geom_boxplot(width = 0.3, outlier.color = NA, alpha = 0.7, linewidth = 0.3) +
  geom_jitter(position = position_jitter(width = 0), size = 1, alpha = 0.3) +
  geom_text(
    data = weight_summary,
    aes(x = Genotype, y = y_pos * 1.2, label = label),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
  annotate("text", x = 3, y = upper_limit * 1.05, label = asterisk, size = 6) +
  scale_fill_manual(values = colors) +
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6)
  ) +
  scale_x_discrete(
    labels = c(
      "WT" = "WT",
      "symrk" = expression(italic("symrk")),
      "ccamk" = expression(italic("ccamk")),
      "nsp1" = expression(italic("nsp1")),
      "nsp2" = expression(italic("nsp2"))
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 1))

write.csv(
  x = weight_summary,
  file = "3_tables/HordeumCSSP_axenic_shoot_fw_ANOVA.csv"
)

# Saving plot
ggsave(
  filename = "2_figures/HordeumCSSP_axenic_shoot_fw.pdf",
  plot = box_plot,
  width = 6,
  height = 6,
  units = "cm"
)

saveRDS(object = box_plot, file = "1_rds_files/HordeumCSSP_axenic_shoot_fw.rds")
