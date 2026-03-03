# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(ggplot2)
library(multcompView)
library(dplyr)

# Load input data.
weight <- read.table("LotusSC_shootfw_input.txt", header=T, sep="\t")
weight$Fresh_weight <- as.numeric(gsub(",", ".", weight$Fresh_weight))

# Set the factor levels and colours for the genotypes.
weight$Genotype <- factor(
  weight$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

colors <- c("WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366",
            "nsp1"="#C0E4EF","nsp2"="#6EA6CD")

# Perform significance analysis using ANOVA and TukeyHSD.
aov_res <- aov(Fresh_weight ~ Genotype, data = weight)
tukey <- TukeyHSD(aov_res)
tukey_pvals <- tukey$Genotype[, "p adj"]
names(tukey_pvals) <- rownames(tukey$Genotype)
letters <- multcompLetters(tukey_pvals)$Letters
letters_df <- data.frame(Genotype = names(letters), label = letters, stringsAsFactors = FALSE)

## Calculate the overall ANOVA p-values and corresponding asterisks for plotting.
p_value <- summary(aov_res)[[1]][["Pr(>F)"]][1]
asterisk <- case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  TRUE ~ NA_character_
)

## Make a summary dataframe for plotting.
weight_summary <- weight %>%
  group_by(Genotype) %>%
  summarise(y_pos = max(Fresh_weight, na.rm=TRUE), .groups="drop") %>%
  left_join(letters_df, by="Genotype")

## Set dummy points for the plot limits.
upper_limit <- max(weight$Fresh_weight, na.rm=TRUE) * 1.2

# Set the main theme for plotting.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text.x = element_text(size = 6, colour = "black"),
  axis.text.y = element_text(size = 6, colour = "black"),
  legend.background=element_blank(),
  legend.key=element_blank(),
  text=element_text(family="sans")
)

# Make the boxplots.
weight_summary$y_pos[2] <- weight_summary$y_pos[2] - 0.015
p1 <- ggplot(weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  geom_text(data=weight_summary, aes(x=Genotype, y=y_pos*1.2, label=label), 
            inherit.aes=FALSE, size = 6/.pt) +
  scale_fill_manual(values=colors) +
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(
    legend.position="none",
    plot.title = element_text(size = 6, face="bold", hjust=0.5),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size = 6)
  ) +
  scale_x_discrete(labels=c(
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  )) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0,0.08))

p1

# Save the plot and significance analysis outputs.
weight$Host <- "Lotus"
weight_summary$Host <- "Lotus"
write.csv(weight, "LotusSC_shoot_fw_ANOVA.csv")
write.csv(weight_summary, "LotusSC_shoot_fw_significance_letters.csv")
ggsave("LotusSC_shootfw_boxplots.pdf", p1, width=5, height=6, units = "cm")
saveRDS(p1, "LotusSC_shootfw_boxplots.rds")
saveRDS(p1, "../../3_final_figures/LotusSC_shootfw_boxplots.rds")
