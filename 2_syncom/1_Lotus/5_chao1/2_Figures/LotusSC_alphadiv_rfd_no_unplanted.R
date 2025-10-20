# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load libraries
library(ggplot2)
library(dplyr)
library(multcompView)
library(forcats)

# Load chao1 and metadata files (note: we have one alpha diversity file where all ASVs were taken into account, and one where only matched ASVs were used (filt)).
alpha <- read.table("../1_Scripts/allASVs/alpha-diversity.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
alpha_filt <- read.table("../1_Scripts/matchedASVsonly/alpha-diversity.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
design <- read.table("LotusCSSP_LjSC_metadata.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)

# For both dataframes:
# 1. add metadata info
# 2. set genotype factor levels
# 3. filter for genotypes of interest
# 4. rename compartments 

create_index <- function(alpha_df) {
  idx <- cbind(alpha_df[,1], design[match(row.names(alpha_df), row.names(design)), ])
  colnames(idx)[1] <- "value"
  idx$Genotype <- factor(idx$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))
  idx <- idx %>% filter(Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2"))
  idx$Compartment <- factor(idx$Compartment,
                            levels = c("rhizo","endo","nod"),
                            labels = c("Rhizosphere","Root","Nodules"))
  return(idx)
}

index <- create_index(alpha)
index_filt <- create_index(alpha_filt)

# Define colors for genotypes and main theme for boxplots.
colors <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

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

# Next, perform ANOVA and Tukey HSD for rhizosphere and root compartment.
## Function:
get_tukey_letters <- function(df){
  counts <- table(df$Genotype)
  valid_genotypes <- names(counts[counts > 0])
  if(length(valid_genotypes) < 2){
    return(data.frame(Genotype = levels(df$Genotype), label = NA_character_, stringsAsFactors = FALSE))
  }
  tryCatch({
    aov_res <- aov(value ~ Genotype, data = df)
    tukey <- TukeyHSD(aov_res)
    tukey_pvals <- tukey$Genotype[, "p adj"]
    names(tukey_pvals) <- rownames(tukey$Genotype)
    letters <- multcompLetters(tukey_pvals)$Letters
    letters <- letters[levels(df$Genotype)]
    data.frame(Genotype = levels(df$Genotype), label = letters, stringsAsFactors = FALSE)
  }, error = function(e){
    data.frame(Genotype = levels(df$Genotype), label = NA_character_, stringsAsFactors = FALSE)
  })
}

compute_anova_summary <- function(df, compartments=c("Rhizosphere","Root")){
  df %>%
    filter(Compartment %in% compartments) %>%
    group_by(Compartment) %>%
    summarise(aov_res = list(aov(value ~ Genotype, data = cur_data())), .groups="drop") %>%
    rowwise() %>%
    mutate(
      p_value = tryCatch(summary(aov_res)[[1]][["Pr(>F)"]][1], error = function(e) NA_real_),
      asterisk = case_when(
        is.na(p_value) ~ NA_character_,
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ NA_character_
      ),
      y_position = max(df$value[df$Compartment == Compartment], na.rm = TRUE) * 1.05
    ) %>%
    ungroup()
}

prepare_plot_summary <- function(df, letters_df){
  df %>%
    group_by(Compartment, Genotype) %>%
    summarise(y_pos = if(all(is.na(value))) NA_real_ else max(value, na.rm=TRUE),
              n = sum(!is.na(value)),
              .groups = "drop") %>%
    left_join(letters_df, by = c("Compartment", "Genotype"))
}

plot_alpha <- function(df, plot_summary, anova_pvals, upper_limits, colors, title="Chao1 index"){
  ggplot(df, aes(x=Genotype, y=value, fill=Genotype)) +
    geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
    geom_jitter(position=position_jitter(width=0), size=3, alpha=0.3) +
    geom_text(data = plot_summary %>% filter(!is.na(label)), aes(x = Genotype, y = y_pos * 1.2, label = label), inherit.aes = FALSE, size = 6) +
    geom_text(data = anova_pvals %>% filter(!is.na(asterisk)), aes(x = 3, y = y_position, label = asterisk), inherit.aes = FALSE, size = 10) +
    facet_wrap(~Compartment, scales="fixed") +
    scale_fill_manual(values = colors) +
    main_theme +
    labs(y = title, x = "") +
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    ) +
    scale_x_discrete(labels=c(
      "WT"="WT",
      "symrk"=expression(italic("symrk")),
      "ccamk"=expression(italic("ccamk")),
      "nsp1"=expression(italic("nsp1")),
      "nsp2"=expression(italic("nsp2"))
    )) +
    expand_limits(y = 0) +
    geom_blank(data = upper_limits, aes(y = ymax), inherit.aes = FALSE)
}

# Expanded function for full workflow:
run_analysis_plot <- function(df, colors, title, output_pdf, output_rds){
  letters_df <- df %>%
    filter(Compartment %in% c("Rhizosphere","Root")) %>%
    group_by(Compartment) %>%
    group_modify(~ get_tukey_letters(.x)) %>%
    ungroup()
  
  plot_summary <- prepare_plot_summary(df, letters_df)
  anova_pvals <- compute_anova_summary(df)
  upper_limits <- df %>%
    group_by(Compartment) %>%
    summarise(ymax = if(all(is.na(value))) 0 else max(value, na.rm=TRUE) * 1.2, .groups="drop")
  
  p <- plot_alpha(df, plot_summary, anova_pvals, upper_limits, colors, title)
  
  print(p)
  ggsave(output_pdf, p, width = 15, height = 6)
  saveRDS(p, output_rds)
}

# Run function for both datasets and save output files.
run_analysis_plot(index,
                  colors,
                  title = "Chao1 index",
                  output_pdf = "LotusSC_chao1_allASVs.pdf",
                  output_rds = "LotusSC_chao1_allASVs.rds")

run_analysis_plot(index_filt,
                  colors,
                  title = "Chao1 index",
                  output_pdf = "LotusSC_chao1_filteredASVs.pdf",
                  output_rds = "LotusSC_chao1_filteredASVs.rds")
