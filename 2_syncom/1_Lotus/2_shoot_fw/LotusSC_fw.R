# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
library(ggplot2)
library(multcompView)
library(dplyr)

# -----------------------------
# Load data
# -----------------------------
weight <- read.table("Input_shoot_freshweights_harvest_magentas.txt", header=T, sep="\t")
weight$X <- NULL
weight$Fresh_weight <- as.numeric(gsub(",", ".", weight$Fresh_weight))

# Select genotypes of interest
weight <- weight %>%
  filter(Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")) %>%
  mutate(
    Genotype = factor(Genotype, levels=c("WT","symrk","ccamk","nsp1","nsp2")),
    Treatment = ifelse(Treatment == "no_SC", "uninoculated",
                       ifelse(Treatment == "Lj_SC", "Lj-SPHERE", Treatment)),
    Treatment = factor(Treatment, levels=c("uninoculated", "Lj-SPHERE"))
  )

# Set colors for genotypes
colors <- c("WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366",
            "nsp1"="#C0E4EF","nsp2"="#6EA6CD")

# -----------------------------
# ANOVA & Tukey letters
# -----------------------------
get_tukey_letters <- function(df){
  aov_res <- aov(Fresh_weight ~ Genotype, data=df)
  tukey <- TukeyHSD(aov_res)
  tukey_pvals <- tukey$Genotype[, "p adj"]
  names(tukey_pvals) <- rownames(tukey$Genotype)
  letters <- multcompLetters(tukey_pvals)$Letters
  letters <- letters[levels(df$Genotype)]
  data.frame(Genotype=levels(df$Genotype), label=letters, stringsAsFactors=FALSE)
}

# Tukey letters per treatment
letters_df <- weight %>%
  group_by(Treatment) %>%
  group_modify(~ get_tukey_letters(.x))

# Summary for plotting letters
weight_summary <- weight %>%
  group_by(Treatment, Genotype) %>%
  summarise(y_pos = max(Fresh_weight, na.rm=TRUE), .groups="drop") %>%
  left_join(letters_df, by=c("Treatment","Genotype")) %>%
  mutate(Treatment = factor(Treatment, levels=c("uninoculated", "Lj-SPHERE")))

# Overall ANOVA p-values
anova_pvals <- weight %>%
  group_by(Treatment) %>%
  summarise(aov_res = list(aov(Fresh_weight ~ Genotype, data = cur_data())), .groups="drop") %>%
  rowwise() %>%
  mutate(
    p_value = summary(aov_res)[[1]][["Pr(>F)"]][1],
    asterisk = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ NA_character_
    ),
    y_position = max(weight$Fresh_weight[weight$Treatment==Treatment], na.rm=TRUE) * 1.15,
    Treatment = factor(Treatment, levels=c("uninoculated", "Lj-SPHERE"))
  )

# Dummy points for plot limits
upper_limits <- data.frame(
  Treatment = levels(weight$Treatment),
  Fresh_weight = sapply(levels(weight$Treatment), function(t) max(weight$Fresh_weight[weight$Treatment==t])*1.2)
) %>%
  mutate(Treatment = factor(Treatment, levels=c("uninoculated", "Lj-SPHERE")))

# -----------------------------
# Plot theme
# -----------------------------
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text.x = element_text(size = 8, colour = "black"),
  axis.text.y = element_text(size = 8, colour = "black"),
  legend.background=element_blank(),
  legend.key=element_blank(),
  text=element_text(family="sans")
)

# -----------------------------
# Plot: Uninoculated + Lj-SPHERE
# -----------------------------
p1 <- ggplot(weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  geom_jitter(position=position_jitter(width=0), size=0.5, alpha=0.3) +
  geom_text(data=weight_summary, aes(x=Genotype, y=y_pos*1.2, label=label), inherit.aes=FALSE, size=6) +
  # geom_text(data=anova_pvals, aes(x=3, y=y_position, label=asterisk), inherit.aes=FALSE, size=10) +
  geom_blank(data=upper_limits, aes(y=Fresh_weight), inherit.aes=FALSE) +
  scale_fill_manual(values=colors) +
  facet_wrap(~Treatment, scales="fixed") +
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(
    legend.position="none",
    strip.text.x = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8)
  ) +
  scale_x_discrete(labels=c(
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  )) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0,0.1))

# Display plot
p1

# Save plot
ggsave("LotusSC_shootfw_incl_uninoc.pdf", p1, width=10, height=6, units = "cm")
saveRDS(p1, file = "LotusSC_shootfw_incl_uninoc.rds")
saveRDS(p1, file = "../10_final_figures/LotusSC_shootfw_incl_uninoc.rds")

# -----------------------------
# Lj-SPHERE only (single plot)
# -----------------------------
weight_Lj <- weight %>% filter(Treatment=="Lj-SPHERE")
weight_summary_Lj <- weight_summary %>% filter(Treatment=="Lj-SPHERE")
anova_pvals_Lj <- anova_pvals %>% filter(Treatment=="Lj-SPHERE")
upper_limits_Lj <- upper_limits %>% filter(Treatment=="Lj-SPHERE")

p2 <- ggplot(weight_Lj, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  geom_jitter(position=position_jitter(width=0), size=0.5, alpha=0.3) +
  geom_text(data=weight_summary_Lj, aes(x=Genotype, y=y_pos*1.2, label=label),
            inherit.aes=FALSE, size=3) +
  # geom_text(data=anova_pvals_Lj, aes(x=3, y=y_position, label=asterisk), inherit.aes=FALSE, size=10) +
  geom_blank(data=upper_limits_Lj, aes(y=Fresh_weight), inherit.aes=FALSE) +
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
  scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0,0.08))

# Display plot
p2

# Save Lj-SPHERE plot
ggsave("LotusSC_shootfw_LjSC_only.pdf", p2, width=5, height=6, units = "cm")
saveRDS(p2, "LotusSC_shootfw_LjSC_only.rds")
saveRDS(p2, "../10_final_figures/LotusSC_shootfw_LjSC_only.rds")
