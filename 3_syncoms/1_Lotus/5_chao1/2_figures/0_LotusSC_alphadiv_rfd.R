# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(dplyr)
library(multcompView)
library(patchwork)
library(ggtext)

# Load chao1 and metadata files.
# Note: we have one file for matched ASVs only, and one where all ASVs are included.
 alpha <- read.table("../1_scripts/allASVs/LotusSC_allASVs_chao1.txt",
                    sep="\t", header=TRUE, row.names=1, check.names=FALSE)
alpha_filt <- read.table("../1_scripts/matchedASVsonly/LotusSC_matchedASVsonly_chao1.txt",
                         sep="\t", header=TRUE, row.names=1, check.names=FALSE)
design <- read.table("../../1_data/LotusSC_metadata.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)

# Add metadata info and to set genotype and compartment factor levels for both datasets.
create_index <- function(alpha_df) {
  idx <- cbind(alpha_df[,1], design[match(row.names(alpha_df), row.names(design)), ])
  colnames(idx)[1] <- "value"
  idx$Genotype <- factor(idx$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))
  idx$Compartment <- factor(idx$Compartment, levels = c("Rhizosphere", "Root", "Nodules"))
  return(idx)
}

index <- create_index(alpha)
index_filt <- create_index(alpha_filt)

# Define the colours for the genotypes.
colors <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

# Set the main theme for plotting.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text = element_text(size = 8, color = "black"),
  legend.text = element_text(size=8, color = "black"),
  legend.key=element_blank(),
  axis.title.y = element_text(size = 8),
  legend.position="none",
  strip.text = element_text(size = 8, color = "black"),
  legend.background=element_blank(),
  plot.title = element_text(size=8, hjust=1)
)

# Perform significance analysis using ANOVA and Tukey HSD for each compartment.
# Note: Nodules only present in WT and therefore no significance analysis for this compartment.

## Write a function to get the Tukey HSD letters, and skip compartments with data for only 1 genotype.
get_tukey_letters <- function(df){
  valid_genotypes <- unique(df$Genotype)
  if(length(valid_genotypes) < 2){
    return(data.frame(Genotype = valid_genotypes, label = NA_character_, Compartment = df$Compartment[1], stringsAsFactors = FALSE))
  }
  tryCatch({
    aov_res <- aov(value ~ Genotype, data = df)
    tukey <- TukeyHSD(aov_res)
    tukey_pvals <- tukey$Genotype[, "p adj"]
    letters <- multcompLetters(tukey_pvals)$Letters
    letters_df <- data.frame(
      Genotype = levels(df$Genotype),
      label = letters[levels(df$Genotype)],
      Compartment = df$Compartment[1],
      stringsAsFactors = FALSE
    )
    return(letters_df)
  }, error = function(e){
    data.frame(Genotype = levels(df$Genotype), label = NA_character_, Compartment = df$Compartment[1], stringsAsFactors = FALSE)
  })
}

# Summarise the ANOVA output and match asterisks to significance levels.
compute_anova_summary <- function(df){
  df %>%
    group_by(Compartment) %>%
    summarise(
      n_genotypes = length(unique(Genotype)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      aov_res = if(n_genotypes > 1) list(aov(value ~ Genotype, data = df[df$Compartment == Compartment, ])) else list(NULL),
      p_value = if(n_genotypes > 1 && !is.null(aov_res[[1]])) {
        s <- summary(aov_res[[1]])
        if(length(s) > 0 && "Pr(>F)" %in% colnames(s[[1]])) s[[1]][["Pr(>F)"]][1] else NA_real_
      } else NA_real_,
      asterisk = case_when(
        is.na(p_value) ~ NA_character_,
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ NA_character_
      ),
      y_position = if(n_genotypes > 1) max(df$value[df$Compartment == Compartment], na.rm = TRUE) * 1.05 else NA_real_
    ) %>%
    ungroup() %>%
    select(Compartment, p_value, asterisk, y_position)
}

# Prepare a plot summary file.
prepare_plot_summary <- function(df, letters_df){
  df %>%
    group_by(Compartment, Genotype) %>%
    summarise(y_pos = if(all(is.na(value))) NA_real_ else max(value, na.rm=TRUE), .groups="drop") %>%
    left_join(letters_df, by=c("Compartment", "Genotype"))
}

# Define a function to make boxplots.
plot_facet <- function(df, plot_summary, anova_pvals, colors, y_limits=NULL, y_breaks=NULL){
  ggplot(df, aes(x=Genotype, y=value, fill=Genotype)) +
    geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
    geom_jitter(position=position_jitter(width=0), size=1, alpha=0.5) +
    geom_text(data = plot_summary %>% filter(!is.na(label)), 
              aes(x=Genotype, y = y_pos * 1.05, label=label), inherit.aes=FALSE, size=3) +
    geom_text(data = anova_pvals %>% filter(!is.na(asterisk)),
              aes(x = 2, y = y_position, label = asterisk), inherit.aes=FALSE, size=4) +
    scale_fill_manual(values = colors) +
    main_theme +
    labs(y="", x="") +
    scale_x_discrete(labels=c(
      "WT"="WT",
      "symrk"=expression(italic("symrk")),
      "ccamk"=expression(italic("ccamk")),
      "nsp1"=expression(italic("nsp1")),
      "nsp2"=expression(italic("nsp2"))
    )) +
    facet_wrap(~Compartment, scales="fixed") +
    {if(!is.null(y_limits)) scale_y_continuous(limits=y_limits, breaks=y_breaks) else NULL}
}

# Define a function that runs the workflow for both datasets, and that makes the final plots.
run_analysis_plot <- function(df, colors, title, output_pdf, output_rds){
  # Define letters.
  letters_df <- do.call(rbind, lapply(levels(df$Compartment), function(comp){
    sub_df <- df[df$Compartment==comp, ]
    if(nrow(sub_df)==0) return(NULL)
    get_tukey_letters(sub_df)
  }))
  
  # Write summary.
  plot_summary <- prepare_plot_summary(df, letters_df)
  
  # Compute ANOVA and p-values.
  anova_pvals <- compute_anova_summary(df)
  
  # Separate the data based on compartments.
  df_main <- df %>% filter(Compartment %in% c("Rhizosphere","Root"))
  df_nod <- df %>% filter(Compartment=="Nodules")
  
  plot_summary_main <- plot_summary %>% filter(Compartment %in% c("Rhizosphere","Root"))
  plot_summary_nod <- plot_summary %>% filter(Compartment=="Nodules")
  
  anova_main <- anova_pvals %>% filter(Compartment %in% c("Rhizosphere","Root"))
  anova_nod <- anova_pvals %>% filter(Compartment=="Nodules")
  
  # Determine a shared y axis.
  y_limits <- c(0, max(df_main$value, na.rm=TRUE) * 1.2)
  y_breaks <- pretty(y_limits)
  
  # Make the main plot (rhizosphere and root).
  p1 <- plot_facet(df_main, plot_summary_main, anova_main, colors, y_limits=y_limits, y_breaks=y_breaks)
  
  # Make a separate plot for the nodules compartment (as only WT data present).
  p2 <- plot_facet(df_nod, plot_summary_nod, anova_nod, colors, y_limits=y_limits, y_breaks=y_breaks) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  # Combine both plots horizontally.
  combined <- p1 + p2 + plot_layout(ncol=2, widths=c(2,0.2))
  
  # Save the plots.
  print(combined)
  ggsave(output_pdf, combined, width=15, height=6)
  saveRDS(combined, output_rds)
  saveRDS(combined, file = file.path("../../../3_final_figures", basename(output_rds)))
}

# Run the function for both datasets (matched ASVs only, and all ASVs).
run_analysis_plot(index, colors, "Chao1 index",
                  "LotusSC_chao1_allASVs_rfd_combined.pdf",
                  "LotusSC_chao1_allASVs_rfd_combined.rds")

run_analysis_plot(index_filt, colors, "Chao1 index",
                  "LotusSC_chao1_matchedASVsonly_rfd_combined.pdf",
                  "LotusSC_chao1_matchedASVsonly_rfd_combined.rds")
