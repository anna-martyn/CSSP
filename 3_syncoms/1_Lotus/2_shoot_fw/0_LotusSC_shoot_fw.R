# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("dplyr", "ggplot2", "multcompView")
for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
weight <- read.table(
  "LotusSC_shootfw_input.txt",
  header = TRUE,
  sep = "\t",
  dec = ","
)

# Setting factor levels and colours for genotypes
weight$Genotype <- factor(
  weight$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

colors <- c(
  "WT"="#A9C289", "symrk"="#FEDA8B", "ccamk"="#FDB366",
  "nsp1"="#C0E4EF", "nsp2"="#6EA6CD"
)

# Hypothesis tests ------------------------------------------------------------
# Hypothesis testing with ANOVA and TukeyHSD
aov_res <- aov(Fresh_weight ~ Genotype, data = weight)
tukey <- TukeyHSD(aov_res)
tukey_pvals <- tukey$Genotype[, "p adj"]
names(tukey_pvals) <- rownames(tukey$Genotype)
letters <- multcompLetters(tukey_pvals)$Letters
letters_df <- data.frame(
  Genotype = names(letters),
  label = letters,
  stringsAsFactors = FALSE
)

# p-values and asterisks
p_value <- summary(aov_res)[[1]][["Pr(>F)"]][1]
asterisk <- case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  TRUE ~ NA_character_
)

# Summary dataframe
weight_summary <- weight %>%
  group_by(Genotype) %>%
  summarise(y_pos = max(Fresh_weight, na.rm = TRUE), .groups = "drop") %>%
  left_join(letters_df, by = "Genotype")

# Plot ------------------------------------------------------------------------
# Setting dummy points for plot limits
upper_limit <- max(weight$Fresh_weight, na.rm = TRUE) * 1.2

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text.x = element_text(size = 6, colour = "black"),
  axis.text.y = element_text(size = 6, colour = "black"),
  legend.background = element_blank(),
  legend.key = element_blank(),
  text = element_text(family = "sans")
)

# Boxplots
weight_summary$y_pos[2] <- weight_summary$y_pos[2] - 0.015
box_plot <- ggplot(weight, aes(x = Genotype, y = Fresh_weight, fill = Genotype)) +
  geom_boxplot(width = 0.3, outlier.color = NA, alpha = 0.7) +
  geom_text(
    data = weight_summary,
    aes(x = Genotype, y = y_pos * 1.2, label = label),
    inherit.aes = FALSE,
    size = 6 / .pt
  ) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 0.08))

# Saving plots and results of hypothesis tests
weight$Host <- "Lotus"
weight_summary$Host <- "Lotus"
write.csv(weight, "LotusSC_shoot_fw_ANOVA.csv")
write.csv(weight_summary, "LotusSC_shoot_fw_significance_letters.csv")
ggsave(
  "LotusSC_shootfw_boxplots.pdf",
  box_plot,
  width = 5,
  height = 6,
  units = "cm"
)
saveRDS(box_plot, "LotusSC_shootfw_boxplots.rds")
saveRDS(box_plot, "../../3_final_figures/LotusSC_shootfw_boxplots.rds")
