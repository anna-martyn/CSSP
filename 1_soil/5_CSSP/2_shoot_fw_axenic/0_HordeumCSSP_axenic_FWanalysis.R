# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(ggplot2)
library(multcompView)
library(dplyr)

# Load data.
weight <- read.table("HordeumCSSP_axenic_fw_input.txt", header=T, sep="\t")

# Set colours and factor levels for genotypes.
colors <- c("WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366",
            "nsp1"="#C0E4EF","nsp2"="#6EA6CD")

weight$Genotype <- factor(
  weight$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

# ANOVA and TukeyHSD
aov_res <- aov(Shoot_fw ~ Genotype, data = weight)
tukey <- TukeyHSD(aov_res)
tukey_pvals <- tukey$Genotype[, "p adj"]
names(tukey_pvals) <- rownames(tukey$Genotype)
letters <- multcompLetters(tukey_pvals)$Letters
letters_df <- data.frame(Genotype = names(letters), label = letters, stringsAsFactors = FALSE)

## Overall ANOVA p-values and asterisks.
p_value <- summary(aov_res)[[1]][["Pr(>F)"]][1]
asterisk <- case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  TRUE ~ NA_character_
)

## Summary for plotting letters.
weight_summary <- weight %>%
  group_by(Genotype) %>%
  summarise(y_pos = max(Shoot_fw, na.rm=TRUE), .groups="drop") %>%
  left_join(letters_df, by="Genotype")

## Dummy points for plot limits.
upper_limit <- max(weight$Shoot_fw, na.rm=TRUE) * 1.2

# Set main theme for plot.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text.x = element_text(size = 8, angle = 45, vjust = 1, 
                             hjust=1, colour = "black"),
  axis.text.y = element_text(size = 8, colour = "black"),
  legend.background=element_blank(),
  legend.key=element_blank(),
  text=element_text(family="sans")
)

# Plot fresh weights.
p1 <- ggplot(weight, aes(x=Genotype, y=Shoot_fw, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  geom_jitter(position=position_jitter(width=0), size=1, alpha=0.3) +
  geom_text(data=weight_summary, aes(x=Genotype, y=y_pos*1.2, label=label), inherit.aes=FALSE, size=8/.pt) +
  annotate("text", x=3, y=upper_limit*1.05, label=asterisk, size=8) +
  scale_fill_manual(values=colors) +
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(
    legend.position="none",
    plot.title = element_text(size=8, face="bold", hjust=0.5),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=8)
  ) +
  scale_x_discrete(labels=c(
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  )) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0,1))

p1

write.csv(weight_summary, "HordeumCSSP_axenic_shoot_fw_ANOVA.csv")

# Save the plot.
ggsave("HordeumCSSP_axenic_shoot_fw.pdf", p1, width=6, height=6, units = "cm")
saveRDS(p1, "HordeumCSSP_axenic_shoot_fw.rds")
saveRDS(p1, "../8_final_figures/HordeumCSSP_axenic_shoot_fw.rds")
