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
weight <- read.table("BarleyCSSP_cerealSC_diaz_shootfw.txt", header=T, sep="\t")

# Select samples for genotypes and treatments of interest. Rename uninoculated WT to WT_uninoc in Genotype column.
weight <- weight %>%
  filter(
    Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2"),
    Inoculum %in% c("uninoc", "SC_only")
  ) %>%
  mutate(
    Inoculum = case_when(
      Inoculum == "uninoc" ~ "uninoculated",
      Inoculum == "SC_only" ~ "HvSPHERE",
      TRUE ~ Inoculum
    ),
    Inoculum = factor(Inoculum, levels = c("uninoculated", "HvSPHERE")),
    Genotype = ifelse(Genotype == "WT" & Inoculum == "uninoculated", "WT_uninoc", Genotype)
  )

# Set colours and factor levels for genotypes.
colors <- c("WT_uninoc"="grey","WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366",
            "nsp1"="#C0E4EF","nsp2"="#6EA6CD")

desired_order <- c("WT_uninoc", "WT", "symrk", "ccamk", "nsp1", "nsp2")
weight$Genotype <- factor(weight$Genotype, levels = desired_order)

# Write function for significance analysis (ANOVA and Tukey HSD).
get_tukey_letters <- function(df){
  df$Genotype <- factor(df$Genotype) # ensure factor levels
  if(length(unique(df$Genotype)) < 2){
    return(data.frame(Genotype = levels(df$Genotype), label = NA_character_, stringsAsFactors=FALSE))
  }
  aov_res <- aov(Fresh_weight ~ Genotype, data=df)
  tukey <- TukeyHSD(aov_res)
  tukey_pvals <- tukey$Genotype[, "p adj"]
  names(tukey_pvals) <- rownames(tukey$Genotype)
  letters <- multcompLetters(tukey_pvals)$Letters
  letters <- letters[levels(df$Genotype)]
  data.frame(Genotype = levels(df$Genotype), label = letters, stringsAsFactors=FALSE)
}

get_anova_sig <- function(df){
  if(length(unique(df$Genotype)) < 2) return(NA)
  aov_res <- aov(Fresh_weight ~ Genotype, data=df)
  p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
  asterisk <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE ~ NA_character_
  )
  y_pos <- max(df$Fresh_weight, na.rm=TRUE) * 1.15
  data.frame(p_value=p_val, asterisk=asterisk, y_position=y_pos)
}

# Set main theme for plots.
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
  text=element_text(family="sans", size=20, colour="black")
)

# Plot all genotypes including uninoculated WT.
letters_all <- get_tukey_letters(weight)
anova_all <- get_anova_sig(weight)
weight_summary_all <- weight %>%
  group_by(Genotype) %>%
  summarise(y_pos = max(Fresh_weight, na.rm=TRUE), .groups="drop") %>%
  left_join(letters_all, by="Genotype")

p_all <- ggplot(weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  geom_jitter(position=position_jitter(width=0), size=3, alpha=0.3) +
  geom_text(data=weight_summary_all, aes(x=Genotype, y=y_pos*1.2, label=label), inherit.aes=FALSE, size=6) +
  geom_text(data=anova_all, aes(x=3.5, y=y_position, label=asterisk), inherit.aes=FALSE, size=10) +
  scale_fill_manual(values=colors) +
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(legend.position="none", axis.title.x=element_blank()) +
  scale_x_discrete(labels=c(
    "WT_uninoc"="WT_uninoc",
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  )) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0,0.8))

p_all

ggsave("HordeumSC_Shootfw_incl_uninoc.pdf", p_all, width=10, height=6)
saveRDS(p_all, "HordeumSC_Shootfw_incl_uninoc.rds")

# Plot without uninoculated WT.
weight_subset <- weight %>% filter(Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2"))
letters_subset <- get_tukey_letters(weight_subset)
anova_subset <- get_anova_sig(weight_subset)
weight_summary_subset <- weight_subset %>%
  group_by(Genotype) %>%
  summarise(y_pos = max(Fresh_weight, na.rm=TRUE), .groups="drop") %>%
  left_join(letters_subset, by="Genotype")

p_subset <- ggplot(weight_subset, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  geom_jitter(position=position_jitter(width=0), size=3, alpha=0.3) +
  geom_text(data=weight_summary_subset, aes(x=Genotype, y=y_pos*1.2, label=label), inherit.aes=FALSE, size=6) +
  geom_text(data=anova_subset, aes(x=3, y=y_position, label=asterisk), inherit.aes=FALSE, size=10) +
  scale_fill_manual(values=colors) +
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(legend.position="none", axis.title.x=element_blank()) +
  scale_x_discrete(labels=c(
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  )) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0,0.8))

p_subset

ggsave("HordeumSC_Shootfw_no_uninoc.pdf", p_subset, width=8, height=6)
saveRDS(p_subset, "HordeumSC_Shootfw_no_uninoc.rds")
