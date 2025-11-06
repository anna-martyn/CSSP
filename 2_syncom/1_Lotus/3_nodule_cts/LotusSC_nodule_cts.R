# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(ggplot2)
library(multcompView)
library(dplyr)

# Load data.
counts <- read.table("LotusSC_nodule_counts.txt", header=T, sep="\t")

# Select pink nodule type and set factor levels for genotypes.
counts <- counts %>%
  filter(
    Nodule_type == "pink"
  ) %>%
  mutate(
    Genotype = factor(Genotype, levels=c("WT","symrk","ccamk","nsp1","nsp2"))
  )

# ANOVA and Tukey HSD.
get_tukey_letters <- function(df){
  if(length(unique(df$Number)) == 1){  # all identical
    letters <- rep("a", length(unique(df$Genotype)))
    return(data.frame(Genotype=levels(df$Genotype), label=letters, stringsAsFactors=FALSE))
  } else {
    aov_res <- aov(Number ~ Genotype, data=df)
    tukey <- TukeyHSD(aov_res)
    tukey_pvals <- tukey$Genotype[, "p adj"]
    names(tukey_pvals) <- rownames(tukey$Genotype)
    letters <- multcompLetters(tukey_pvals)$Letters
    letters <- letters[levels(df$Genotype)]
    return(data.frame(Genotype=levels(df$Genotype), label=letters, stringsAsFactors=FALSE))
  }
}

letters_df <- get_tukey_letters(counts)

# Determine overall ANOVA p-values.
if(length(unique(counts$Number)) == 1){
  asterisk <- NA
} else {
  aov_res <- aov(Number ~ Genotype, data=counts)
  p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
  asterisk <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else NA
}

# Set y_position for asterisks.
y_max <- max(counts$Number, na.rm=TRUE)
y_pos_asterisk <- y_max * 1.05

# Make summary for plotting letters.
letters_df$y_pos <- counts %>%
  group_by(Genotype) %>%
  summarise(y_pos = max(Number, na.rm=TRUE)) %>%
  pull(y_pos)

# Set main theme for plotting.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text.x = element_text(size = 20, colour = "black"),
  axis.text.y = element_text(size = 20, colour = "black"),
  legend.background=element_blank(),
  legend.key=element_blank(),
  text=element_text(family="sans")
)

# Set colour for genotypes.
colors <- c("WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366",
            "nsp1"="#C0E4EF","nsp2"="#6EA6CD")

# Additional info for asterisk position.
asterisk_df <- data.frame(
  x = 3,
  y = y_pos_asterisk,
  label = asterisk
  
)
# Make plot.
p <- ggplot(counts, aes(x=Genotype, y=Number, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  geom_jitter(position=position_jitter(width=0), size=3, alpha=0.3) +
  # geom_text(data=letters_df, aes(x=Genotype, y=y_pos*1.05, label=label), inherit.aes=FALSE, size=6) +
  # geom_text(data = asterisk_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 10)+
  scale_fill_manual(values=colors) +
  main_theme +
  ylab("Number of pink nodules/plant") +
  theme(
    legend.position="none",
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=20)
  ) +
  scale_x_discrete(labels=c(
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  )) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0, y_max*1.2))

p

# Save plot.
ggsave("LotusSC_nodule_cts.pdf", p, width=5, height=6)
saveRDS(p, "LotusSC_nodule_cts.rds")
