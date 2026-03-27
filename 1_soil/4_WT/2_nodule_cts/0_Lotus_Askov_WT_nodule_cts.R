# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "dplyr", "multcompView")
for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
nod <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_nodule_cts.txt",
  header = TRUE,
  sep = "\t"
)

# Factor levels
nod$Soil_type <- factor(nod$Soil_type, levels = c("NPK", "PK", "UF"))

# Soil colours
colors <- data.frame(
  group = c("NPK", "PK", "UF"),
  colors = c("#6F944F","#B2563C","#3C7D82")
)

# Hypothesis tests ------------------------------------------------------------
# Hypothesis testing using ANOVA and Tukey HSD with letters and asterisks
get_tukey_letters <- function(df){
  aov_res <- aov(pink ~ Soil_type, data = df)
  tukey <- TukeyHSD(aov_res)
  
  pvals <- tukey$Soil_type[, "p adj"]
  names(pvals) <- rownames(tukey$Soil_type)
  
  letters <- multcompLetters(pvals)$Letters
  letters[levels(df$Soil_type)]
}

letters <- get_tukey_letters(nod)

letters_df <- data.frame(
  Soil_type = levels(nod$Soil_type),
  label = letters,
  stringsAsFactors = FALSE
)

aov_res <- aov(pink ~ Soil_type, data = nod)
anova_p <- summary(aov_res)[[1]][["Pr(>F)"]][1]
asterisk <- case_when(
  anova_p < 0.001 ~ "***",
  anova_p < 0.01  ~ "**",
  anova_p < 0.05  ~ "*",
  TRUE ~ NA_character_
)

# Plot ------------------------------------------------------------------------
y_positions <- nod %>%
  group_by(Soil_type) %>%
  summarise(y_pos = max(pink, na.rm = TRUE), .groups = "drop")

letters_df <- left_join(letters_df, y_positions, by = "Soil_type")

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text.x = element_text(size = 6, colour = "black"),
  axis.text.y = element_text(size = 6), colour = "black",
  legend.background = element_blank(),
  legend.key = element_blank(),
  text = element_text(family = "sans")
)

letters_df$y_pos[c(1,3)] <- c(9.5, 10.5)
box_plot <- ggplot(nod, aes(x = Soil_type, y = pink, fill = Soil_type)) +
  geom_boxplot(width = 0.3, alpha = 0.7, outlier.size = 1.5)+
  scale_fill_manual(values = as.character(colors$color)) +
  geom_text(
    data = letters_df, 
    mapping = aes(x = Soil_type, y = y_pos*1.1, label = label),
    inherit.aes = FALSE,
    size = 6/.pt
  ) +
  main_theme +
  ylab("Pink nodule counts/plant")+
  scale_y_continuous()+
  theme(legend.position = "none", 
        strip.text.x = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6, colour = "black"),
        legend.key.size = unit(1,"cm"))

# Save plot
ggsave(
  filename = "2_figures/Lotus_Askov_WT_nodule_cts.pdf",
  plot = box_plot,
  width = 3,
  height = 6,
  units = "cm"
)

saveRDS(object = box_plot, file = "1_rds_files/Lotus_Askov_WT_nodule_cts.rds")
