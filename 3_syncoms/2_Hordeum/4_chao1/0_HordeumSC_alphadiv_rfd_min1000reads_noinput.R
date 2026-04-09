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
## Alpha diversity estimated from all ASVs
alpha <- read.table(
  file = "0_chao1_qiime/allASVs/HordeumSC_allASVs_chao1.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
## Alpha diversity estimated from only ASVs matched to SynCom sequences
alpha_filt <- read.table(
  file = "0_chao1_qiime/matchedASVsonly/HordeumSC_matchedASVsonly_chao1.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
design <- read.table(
  file = "../1_data/HordeumSC_metadata.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Function to merge diversity data and metadata
create_index <- function(alpha_df){
  idx <- cbind(
    alpha_df[, 1],
    design[match(row.names(alpha_df), row.names(design)), ]
  )
  colnames(idx)[1] <- "value"
  idx$Genotype <- factor(
    idx$Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
  idx$Compartment <- factor(idx$Compartment, levels = c("Rhizosphere", "Root"))
  return(idx)
}

index <- create_index(alpha)
index_filt <- create_index(alpha_filt)

# Genotypes colours
colors <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

# Plot ------------------------------------------------------------------------
# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text.y = element_text(size = 8, color = "black"),
  axis.text.x = element_text(hjust = 1, angle = 30, size = 8, color = "black"),
  legend.text = element_text(size = 8, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 8),
  legend.position = "none",
  strip.text = element_text(size = 8, color = "black"),
  legend.background = element_blank(),
  plot.title = element_text(size = 8, hjust = 1)
)

# Test for WT vs. mutant effects with ANOVA and Tukey HSD for each compartment
# Note: Nodules only present in WT, so no hypothesis testing is carried out in
# that compartment

# Function that produces letters ith Tukey HSD, and skip compartments with data
# for only 1 genotype.
get_tukey_letters <- function(df){
  valid_genotypes <- unique(df$Genotype)
  if (length(valid_genotypes) < 2) {
    return(data.frame(
      Genotype = valid_genotypes,
      label = NA_character_,
      Compartment = df$Compartment[1],
      stringsAsFactors = FALSE
    ))
  }
  tryCatch(
    {
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
    },
    error = function(e) {
      data.frame(
        Genotype = levels(df$Genotype),
        label = NA_character_,
        Compartment = df$Compartment[1],
        stringsAsFactors = FALSE
      )
    }
  )
}

# Summarising ANOVA results and match asterisks to significance levels
compute_anova_summary <- function(df){
  df %>%
    group_by(Compartment) %>%
    summarise(
      n_genotypes = length(unique(Genotype)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      aov_res = if (n_genotypes > 1) {
        list(aov(value ~ Genotype, data = df[df$Compartment == Compartment, ]))
      } else {
        list(NULL)
      },
      p_value = if (n_genotypes > 1 && !is.null(aov_res[[1]])) {
        s <- summary(aov_res[[1]])
        if (length(s) > 0 && "Pr(>F)" %in% colnames(s[[1]])) {
          s[[1]][["Pr(>F)"]][1]
        } else {
          NA_real_
        }
      } else {
        NA_real_
      },
      asterisk = case_when(
        is.na(p_value) ~ NA_character_,
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ NA_character_
      ),
      y_position = if (n_genotypes > 1) {
        max(df$value[df$Compartment == Compartment], na.rm = TRUE) * 1.05
      } else {
        NA_real_
      }
    ) %>%
    ungroup() %>%
    select(Compartment, p_value, asterisk, y_position)
}

# Preparing plot summary
prepare_plot_summary <- function(df, letters_df){
  df %>%
    group_by(Compartment, Genotype) %>%
    summarise(
      y_pos = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(letters_df, by = c("Compartment", "Genotype"))
}

# Function to produce boxplots
plot_facet <- function(
  df,
  plot_summary,
  anova_pvals,
  colors,
  y_limits = NULL,
  y_breaks = NULL
) {
  ggplot(df, aes(x = Genotype, y = value, fill = Genotype)) +
    geom_boxplot(width = 0.3, outlier.color = NA, alpha = 0.7) +
    geom_jitter(position = position_jitter(width = 0), size = 1, alpha = 0.5) +
    geom_text(
      data = plot_summary %>% filter(!is.na(label)),
      aes(x = Genotype, y = y_pos * 1.05, label = label),
      inherit.aes = FALSE,
      size = 3
    ) +
    geom_text(
      data = anova_pvals %>% filter(!is.na(asterisk)),
      aes(x = 2, y = y_position, label = asterisk),
      inherit.aes = FALSE,
      size = 4
    ) +
    scale_fill_manual(values = colors) +
    main_theme +
    labs(y = "", x = "") +
    scale_x_discrete(
      labels = c(
        "WT" = "WT",
        "symrk" = expression(italic("symrk")),
        "ccamk" = expression(italic("ccamk")),
        "nsp1" = expression(italic("nsp1")),
        "nsp2" = expression(italic("nsp2"))
      )
    ) +
    facet_wrap(~Compartment, scales = "fixed") +
    {
      if (!is.null(y_limits)) {
        scale_y_continuous(limits = y_limits, breaks = y_breaks)
      } else {
        NULL
      }
    }
}

# Function for runing workflow for both datasets, and that make final plots
run_analysis_plot <- function(df, colors, title, output_pdf, output_rds){
  # Letters
  letters_df <- do.call(
    rbind,
    lapply(levels(df$Compartment), function(comp) {
      sub_df <- df[df$Compartment == comp, ]
      if (nrow(sub_df) == 0) {
        return(NULL)
      }
      get_tukey_letters(sub_df)
    })
  )

  # Summary
  plot_summary <- prepare_plot_summary(df, letters_df)

  # ANOVA
  anova_pvals <- compute_anova_summary(df)

  # # Separating data by compartments
  # df_main <- df %>% filter(Compartment %in% c("Rhizosphere","Root"))
  #
  # plot_summary_main <- plot_summary %>% filter(Compartment %in% c("Rhizosphere","Root"))
  #
  # anova_main <- anova_pvals %>% filter(Compartment %in% c("Rhizosphere","Root"))

  # Determining shared y axis
  y_limits <- c(0, max(df$value, na.rm = TRUE) * 1.2)
  y_breaks <- pretty(y_limits)

  # Main plot (rhizosphere and root)
  p1 <- plot_facet(
    df,
    plot_summary,
    anova_pvals,
    colors,
    y_limits = y_limits,
    y_breaks = y_breaks
  ) +
    ggtitle("Hordeum") +
    theme(plot.title = element_text(face = "bold", size = 8, hjust = 0))

  # Saving plots
  print(p1)
  ggsave(output_pdf, p1, width = 15, height = 6)
  saveRDS(p1, output_rds)
}

# Running the function for both datasets (matched ASVs and all ASVs)
run_analysis_plot(
  df = index,
  colors = colors,
  title = "Chao1 index",
  output_pdf = "2_figures/HordeumSC_chao1_allASVs_rfd.pdf",
  output_rds = "1_rds_files/HordeumSC_chao1_allASVs_rfd.rds"
)

run_analysis_plot(
  df = index_filt,
  colors = colors,
  title = "Chao1 index",
  output_pdf = "2_figures/HordeumSC_chao1_matchedASVsonly_rfd.pdf",
  output_rds = "1_rds_files/HordeumSC_chao1_matchedASVsonly_rfd.rds"
)
