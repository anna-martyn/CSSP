# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
pkg <- c(
  "ggplot2", "dplyr", "scales", "multcompView", "grid", "car"
)
for(pk in pkg){
  library(pk, character.only = T)
}

# Load the shoot fresh weight input file and fix the decimals.
weight <- read.table("LotusHordeum_AskovSoils_shootfw.txt",
                     header = TRUE, sep = "\t")
weight$Fresh_weight <- as.numeric(gsub(",", ".", weight$Fresh_weight))

# Set the factor levels for the plot.
weight$Soil_type <- factor(weight$Soil_type, levels=c("NPK","PK","UF"))
weight$Plant_species <- factor(weight$Plant_species, levels=c("Lotus","Hordeum"))

# Set the colours for the plot.
colors <- data.frame(group=c("NPK","PK","UF"),
                     colors=c("#6F944F","#B2563C","#3C7D82"))

# Define the function for significance analysis using ANOVA and Tukey and generate letters.
get_tukey_letters <- function(df){
  aov_res <- aov(Fresh_weight ~ Soil_type, data=df)
  tukey <- TukeyHSD(aov_res)

  # Make a vector for the p-values.
  tukey_pvals <- tukey$Soil_type[, "p adj"]
  names(tukey_pvals) <- rownames(tukey$Soil_type)

  letters <- multcompLetters(tukey_pvals)$Letters

  # Reorder the letters to match the factor levels.
  letters_ordered <- letters[levels(df$Soil_type)]

  data.frame(
    Soil_type = levels(df$Soil_type),
    label = letters_ordered,
    stringsAsFactors = FALSE
  )
}

# Decipher the Tukey letters for each plant species.
letters_df <- weight %>%
  group_by(Plant_species) %>%
  group_modify(~ get_tukey_letters(.x))

# Merge the letters with the plotting data.
weight_summary <- weight %>%
  group_by(Plant_species, Soil_type) %>%
  summarise(y_pos = max(Fresh_weight, na.rm=TRUE), .groups="drop") %>%
  left_join(letters_df, by=c("Plant_species","Soil_type"))

# Compute the overall ANOVA p-values for each plant species.
anova_pvals <- weight %>%
  group_by(Plant_species) %>%
  summarise(
    aov_res = list(aov(Fresh_weight ~ Soil_type, data = cur_data())),
    .groups="drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = summary(aov_res)[[1]][["Pr(>F)"]][1],
    asterisk = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ NA_character_
    ),
    y_position = max(weight$Fresh_weight[weight$Plant_species==Plant_species], na.rm=TRUE) * 1.15
  )

# Add dummy points for upper limits for each plant species to allow room for the letters and asterisks in the plot.
upper_limits <- data.frame(
  Plant_species = c("Lotus", "Hordeum"),
  Fresh_weight = c(
    max(weight$Fresh_weight[weight$Plant_species=="Lotus"], na.rm=TRUE) * 1.2,
    max(weight$Fresh_weight[weight$Plant_species=="Hordeum"], na.rm=TRUE) * 1.2
  ),
  Soil_type = "NPK"
)

# Set the factor levels again.
weight$Plant_species <- factor(weight$Plant_species, levels=c("Lotus","Hordeum"))
weight_summary$Plant_species <- factor(weight_summary$Plant_species, levels=c("Lotus","Hordeum"))
upper_limits$Plant_species <- factor(upper_limits$Plant_species, levels=c("Lotus","Hordeum"))
anova_pvals$Plant_species <- factor(anova_pvals$Plant_species, levels=c("Lotus","Hordeum"))

# Set the main theme for plotting.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text.x = element_text(size = 6, colour = "black"),
                    axis.text.y = element_text(size = 6, colour = "black"),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

# Make the plot.
weight_summary$y_pos[5:6] <- c(0.9, 0.6)
p <- ggplot(weight, aes(x=Soil_type, y=Fresh_weight, fill=Soil_type)) +
  geom_boxplot(width=0.3, alpha=0.7, outlier.size = 1.5) +
  geom_text(data=weight_summary, aes(x=Soil_type, y=y_pos * 1.2, label=label),
            inherit.aes=FALSE, size = 6/.pt) +
  geom_blank(data=upper_limits, aes(y=Fresh_weight)) +
  scale_fill_manual(values=colors$color) +
  facet_wrap(~Plant_species, scales="free_y") +
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 6, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6, colour = "black")
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  )

p

# Save the plot.
ggsave("LotusHordeum_Askov_WT_shootfw_boxplots.pdf", p, width=3.5, height=6, units = "cm")
saveRDS(p, file = "LotusHordeum_Askov_WT_shootfw_boxplots.rds")
saveRDS(p, file = "../7_final_figures/LotusHordeum_Askov_WT_shootfw_boxplots.rds")

# Function to perform Tukey & save results
perform_anova_tukey <- function(df, species_name){
  df_species <- df %>% filter(Plant_species==species_name)
  aov_res <- aov(Fresh_weight ~ Soil_type, data=df_species)
  print(anova(aov_res))
  pairwise <- TukeyHSD(aov_res)
  pairwise[["Soil_type"]] <- na.omit(pairwise[["Soil_type"]])
  letters <- multcompLetters(pairwise$Soil_type[,4])$Letters
  letters_df <- data.frame(Soil_type=names(letters), Letters=letters, stringsAsFactors=FALSE)
  write.csv(letters_df, paste0("TukeyHSD_", species_name, "_WT_CSSP.csv"), row.names=FALSE)
}

perform_anova_tukey(weight, "Lotus")
perform_anova_tukey(weight, "Hordeum")

