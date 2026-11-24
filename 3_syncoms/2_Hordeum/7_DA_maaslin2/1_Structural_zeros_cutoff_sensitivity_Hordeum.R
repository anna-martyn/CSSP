# Sensitivity analysis of the structural zero detection used in the SynCom
# differential abundance analysis (0_SC_analysis_Hordeum.R).
# Only the structural zero part is run (no MaAsLin2), for a grid of read and
# The cutoff used in the main analysis is min_reads = 20, min_present_reps = 2.
# As in the main analysis, structural zero calls are not restricted to ASVs
# passing MaAsLin2's filtering.

# Setup ------------------------------------------------------------------------
# Clean up
options(warn = -1)
rm(list = ls())

# Loading packages
pkg <- c("data.table", "ggplot2")

for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading file with structural zero function. The copy used for the soil
# analysis is sourced, as it is identical to the local Structural_zeros.R except
# for 'drop = FALSE' when subsetting. Without it Structural_zeros3 fails when a
# single ASV is a structural zero candidate, which happens for some cutoffs.
source("../../../1_soil/6_CSSP/7_DA_analysis/0_files/Structural_zeros.R")

# Loading data
design <- read.table(
  "../1_data/HordeumSC_metadata.txt",
  header = TRUE,
  sep = "\t"
)
asv_table <- read.table(
  "../1_data/HordeumSC_ASVtable.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)
taxonomy <- fread("../1_data/CerealSC_taxonomy_May23.txt", fill = TRUE)
taxonomy[is.na(Order) | Order == "", Order := "Unknown"]

# Order colours (same as in the main analysis figure)
colors_orders <- fread("../../../0_files/Bacterial_order_colors.csv")

# Keep only ASVs matched to SynCom
asv_table_matched <- as.matrix(
  asv_table[grepl("_", rownames(asv_table)), , drop = FALSE]
)

# Library sizes (needed for Structural_zeros3)
design$library_size <- colSums(asv_table_matched[, design$SampleID])

# Vector of mutants
mutants <- c("symrk", "ccamk", "nsp1", "nsp2")
compartments <- c("Rhizosphere", "Root")

# Cutoff grid
cutoffs <- CJ(min_reads = seq(0, 40, 10), min_present_reps = 1:3)
baseline <- list(min_reads = 20, min_present_reps = 2)

# Splitting data by compartment ------------------------------------------------
asv_tables_subsets <- lapply(compartments, function(current_compartment){
  samples <- design$SampleID[design$Compartment == current_compartment]

  asv_table_sub <- asv_table_matched[, samples, drop = FALSE]
  design_sub <- design[design$Compartment == current_compartment, ]
  rownames(design_sub) <- design_sub$SampleID
  design_sub$SampleID <- NULL
  design_sub$Genotype <- factor(design_sub$Genotype, levels = c("WT", mutants))

  # Mean relative abundance of each ASV within each genotype
  ra_table <- t(t(asv_table_sub) / colSums(asv_table_sub))
  mean_ra <- sapply(
    levels(design_sub$Genotype),
    function(g) rowMeans(ra_table[, design_sub$Genotype == g, drop = FALSE])
  )

  list(
    asv_table = asv_table_sub,
    design = design_sub,
    mean_ra = mean_ra,
    compartment = current_compartment
  )
})

# Structural zero analysis for each compartment and cutoff ---------------------
struc_zero_results <- lapply(asv_tables_subsets, function(X){
  res_cutoffs <- lapply(seq_len(nrow(cutoffs)), function(j){
    current_min_reads <- cutoffs$min_reads[j]
    current_min_reps <- cutoffs$min_present_reps[j]

    da_struc_zero <- Structural_zeros3(
      X$asv_table,
      X$design,
      group = "Genotype",
      ref = "WT",
      min_reads = current_min_reads,
      min_present_reps = current_min_reps
    )

    # One row per ASV detected as DA by structural zero, with direction
    # (absent in mutant = depleted, absent in WT = enriched)
    res <- lapply(mutants, function(current_mutant){
      asvs <- da_struc_zero$struc_zero_DA[[current_mutant]]
      if(length(asvs) == 0){
        return(NULL)
      }
      zero_tab <- da_struc_zero$struc_zero_table[
        asvs,
        c(current_mutant, "WT"),
        drop = FALSE
      ]
      data.table(
        ASVid = asvs,
        Genotype = current_mutant,
        Direction = ifelse(
          apply(zero_tab, 1, diff) == 1,
          "Enriched",
          "Depleted"
        ),
        WT_mean_RA = X$mean_ra[asvs, "WT"],
        Mutant_mean_RA = X$mean_ra[asvs, current_mutant]
      )
    })
    res <- rbindlist(res)

    if(nrow(res) == 0){
      res <- data.table(
        ASVid = character(0),
        Genotype = character(0),
        Direction = character(0),
        WT_mean_RA = numeric(0),
        Mutant_mean_RA = numeric(0)
      )
    }

    res[, ":="(
      min_reads = current_min_reads,
      min_present_reps = current_min_reps
    )]
    return(res)
  })

  res_cutoffs <- rbindlist(res_cutoffs)
  res_cutoffs[, Compartment := X$compartment]
  return(res_cutoffs)
})

struc_zero_res <- rbindlist(struc_zero_results)

# Adding taxonomy
tax_cols <- intersect(
  c("ASVid", "Phylum", "Class", "Order", "Family", "Genus"),
  colnames(taxonomy)
)
struc_zero_res <- merge(
  struc_zero_res,
  taxonomy[, ..tax_cols],
  by = "ASVid",
  all.x = TRUE
)
struc_zero_res[is.na(Order), Order := "Unknown"]
setcolorder(
  struc_zero_res,
  c("min_reads", "min_present_reps", "Compartment", "Genotype",
    "ASVid", "Direction", "WT_mean_RA", "Mutant_mean_RA")
)

# Summarising number and WT RA of DA ASVs per cutoff ---------------------------
# Full grid so that combinations without any DA ASVs appear with zero counts
summary_grid <- CJ(
  min_reads = cutoffs$min_reads,
  min_present_reps = cutoffs$min_present_reps,
  Compartment = compartments,
  Genotype = mutants,
  unique = TRUE
)

# Cumulative mean RA in WT of the DA ASVs. ASVs enriched in the mutant are
# absent from WT and contribute zero, so their cumulative mean RA in the mutant
# is given as well.
counts <- struc_zero_res[,
  .(
    DA_ASVs = .N,
    Enriched = sum(Direction == "Enriched"),
    Depleted = sum(Direction == "Depleted"),
    WT_cum_RA = sum(WT_mean_RA),
    Mutant_cum_RA_enriched = sum(Mutant_mean_RA[Direction == "Enriched"])
  ),
  by = .(min_reads, min_present_reps, Compartment, Genotype)
]

# Overlap with ASVs detected using the cutoff of the main analysis
baseline_res <- struc_zero_res[
  min_reads == baseline$min_reads &
    min_present_reps == baseline$min_present_reps,
  .(Compartment, Genotype, ASVid, in_baseline = TRUE)
]
overlap <- merge(
  struc_zero_res,
  baseline_res,
  by = c("Compartment", "Genotype", "ASVid"),
  all.x = TRUE
)[,
  .(Shared_with_baseline = sum(!is.na(in_baseline))),
  by = .(min_reads, min_present_reps, Compartment, Genotype)
]
baseline_n <- baseline_res[, .(Baseline_DA_ASVs = .N), by = .(Compartment, Genotype)]

struc_zero_summary <- merge(
  summary_grid,
  counts,
  by = c("min_reads", "min_present_reps", "Compartment", "Genotype"),
  all.x = TRUE
)
struc_zero_summary <- merge(
  struc_zero_summary,
  overlap,
  by = c("min_reads", "min_present_reps", "Compartment", "Genotype"),
  all.x = TRUE
)
struc_zero_summary <- merge(
  struc_zero_summary,
  baseline_n,
  by = c("Compartment", "Genotype"),
  all.x = TRUE
)
count_cols <- c(
  "DA_ASVs", "Enriched", "Depleted", "Shared_with_baseline", "Baseline_DA_ASVs"
)
struc_zero_summary[,
  (count_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0L, as.integer(x))),
  .SDcols = count_cols
]
ra_cols <- c("WT_cum_RA", "Mutant_cum_RA_enriched")
struc_zero_summary[,
  (ra_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)),
  .SDcols = ra_cols
]

struc_zero_summary[, ":="(
  Genotype = factor(Genotype, levels = mutants),
  Compartment = factor(Compartment, levels = compartments)
)]
setcolorder(
  struc_zero_summary,
  c("min_reads", "min_present_reps", "Compartment", "Genotype",
    count_cols, ra_cols)
)
struc_zero_summary <- struc_zero_summary[
  order(min_reads, min_present_reps, Compartment, Genotype)
]

# Order-level composition of DA ASVs per cutoff --------------------------------
# Number of DA ASVs and their cumulative mean RA in WT per order, along with
# the proportion of DA ASVs (within cutoff, compartment and genotype)
# belonging to each order
order_comp <- struc_zero_res[,
  .(
    DA_ASVs = .N,
    Enriched = sum(Direction == "Enriched"),
    Depleted = sum(Direction == "Depleted"),
    WT_cum_RA = sum(WT_mean_RA)
  ),
  by = .(min_reads, min_present_reps, Compartment, Genotype, Order)
]
order_comp[,
  Prop_DA_ASVs := DA_ASVs / sum(DA_ASVs),
  by = .(min_reads, min_present_reps, Compartment, Genotype)
]
order_comp[, ":="(
  Genotype = factor(Genotype, levels = mutants),
  Compartment = factor(Compartment, levels = compartments),
  Order = factor(Order, levels = colors_orders$Order)
)]
order_comp <- order_comp[
  order(min_reads, min_present_reps, Compartment, Genotype, Order)
]

# Same composition pooled across compartments and genotypes. The cumulative
# mean RA in WT is averaged over the compartment-genotype combinations
n_combinations <- length(compartments) * length(mutants)
order_comp_pooled <- struc_zero_res[,
  .(DA_ASVs = .N, WT_cum_RA = sum(WT_mean_RA) / n_combinations),
  by = .(min_reads, min_present_reps, Order)
]
order_comp_pooled[,
  Prop_DA_ASVs := DA_ASVs / sum(DA_ASVs),
  by = .(min_reads, min_present_reps)
]
order_comp_pooled[, Order := factor(Order, levels = colors_orders$Order)]
order_comp_pooled <- order_comp_pooled[order(min_reads, min_present_reps, Order)]

# Labels for cutoff combinations used on the x-axis of composition plots
cutoff_levels <- cutoffs[, paste0(min_reads, "/", min_present_reps)]
order_comp[, Cutoff := factor(
  paste0(min_reads, "/", min_present_reps),
  levels = cutoff_levels
)]
order_comp_pooled[, Cutoff := factor(
  paste0(min_reads, "/", min_present_reps),
  levels = cutoff_levels
)]
baseline_label <- paste0(baseline$min_reads, "/", baseline$min_present_reps)

# Exporting tables -------------------------------------------------------------
fwrite(
  struc_zero_summary,
  "3_tables/HordeumSynCom_structural_zeros_cutoff_overview.csv"
)
fwrite(
  struc_zero_res[order(min_reads, min_present_reps, Compartment, Genotype, ASVid)],
  "3_tables/HordeumSynCom_structural_zeros_cutoff_results.csv"
)
fwrite(
  order_comp[, !"Cutoff"],
  "3_tables/HordeumSynCom_structural_zeros_cutoff_order_composition.csv"
)
fwrite(
  order_comp_pooled[, !"Cutoff"],
  "3_tables/HordeumSynCom_structural_zeros_cutoff_order_composition_pooled.csv"
)

# Plots ------------------------------------------------------------------------
theme_cutoff <- theme_bw() +
  theme(
    text = element_text(size = 6, family = "Helvetica", colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    strip.text = element_text(size = 6, face = "bold"),
    strip.background = element_rect(colour = NA),
    legend.key.size = unit(0.3, "cm"),
    panel.grid.minor = element_blank()
  )

# Number of DA ASVs across cutoffs
p_cutoffs <- ggplot(
  data = struc_zero_summary,
  mapping = aes(
    x = min_reads,
    y = DA_ASVs,
    colour = factor(min_present_reps),
    group = factor(min_present_reps)
  )
) +
  geom_vline(xintercept = baseline$min_reads, linetype = "dashed", colour = "grey60") +
  geom_line() +
  geom_point(size = 0.8) +
  facet_grid(Compartment ~ Genotype, scales = "free_y") +
  scale_x_continuous(breaks = unique(cutoffs$min_reads)) +
  labs(
    x = "Minimum number of reads",
    y = "ASVs DA by structural zero",
    colour = "Minimum number\nof samples",
    title = "Hordeum SynCom"
  ) +
  theme_cutoff

# Cumulative mean RA in WT of DA ASVs across cutoffs
p_ra <- ggplot(
  data = struc_zero_summary,
  mapping = aes(
    x = min_reads,
    y = WT_cum_RA,
    colour = factor(min_present_reps),
    group = factor(min_present_reps)
  )
) +
  geom_vline(xintercept = baseline$min_reads, linetype = "dashed", colour = "grey60") +
  geom_line() +
  geom_point(size = 0.8) +
  facet_grid(Compartment ~ Genotype, scales = "free_y") +
  scale_x_continuous(breaks = unique(cutoffs$min_reads)) +
  labs(
    x = "Minimum number of reads",
    y = "Cumulative mean RA in WT of ASVs DA by structural zero",
    colour = "Minimum number\nof samples",
    title = "Hordeum SynCom"
  ) +
  theme_cutoff

# Order composition (proportion of DA ASVs) per compartment and genotype
p_order <- ggplot(
  data = order_comp,
  mapping = aes(x = Cutoff, y = Prop_DA_ASVs, fill = Order)
) +
  geom_col(width = 0.9) +
  facet_grid(Compartment ~ Genotype) +
  scale_fill_manual(
    values = colors_orders$Color,
    breaks = colors_orders$Order,
    drop = TRUE
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = paste0(
      "Minimum number of reads / minimum number of samples (main analysis: ",
      baseline_label, ")"
    ),
    y = "Proportion of ASVs DA by structural zero",
    fill = "Bacterial order",
    title = "Hordeum SynCom"
  ) +
  theme_cutoff +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

# Order composition pooled across compartments and genotypes, shown both as
# proportion of DA ASVs and as cumulative mean RA in WT
order_comp_pooled_long <- melt(
  order_comp_pooled,
  id.vars = c("Cutoff", "Order"),
  measure.vars = c("Prop_DA_ASVs", "WT_cum_RA"),
  variable.name = "Measure"
)
order_comp_pooled_long[, Measure := factor(
  fifelse(
    Measure == "Prop_DA_ASVs",
    "Proportion of DA ASVs",
    "Cumulative mean RA in WT (average)"
  ),
  levels = c("Proportion of DA ASVs", "Cumulative mean RA in WT (average)")
)]
p_order_pooled <- ggplot(
  data = order_comp_pooled_long,
  mapping = aes(x = Cutoff, y = value, fill = Order)
) +
  geom_col(width = 0.9) +
  facet_wrap(~ Measure, scales = "free_y") +
  scale_fill_manual(
    values = colors_orders$Color,
    breaks = colors_orders$Order,
    drop = TRUE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = paste0(
      "Minimum number of reads / minimum number of samples (main analysis: ",
      baseline_label, ")"
    ),
    y = NULL,
    fill = "Bacterial order",
    title = "Hordeum SynCom (pooled across compartments and genotypes)"
  ) +
  theme_cutoff +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

# Saving plots -----------------------------------------------------------------
ggsave(
  "2_figures/HordeumSynCom_structural_zeros_cutoff_sensitivity.pdf",
  p_cutoffs,
  width = 18,
  height = 9,
  units = "cm"
)
ggsave(
  "2_figures/HordeumSynCom_structural_zeros_cutoff_WT_RA.pdf",
  p_ra,
  width = 18,
  height = 9,
  units = "cm"
)
ggsave(
  "2_figures/HordeumSynCom_structural_zeros_cutoff_order_composition.pdf",
  p_order,
  width = 24,
  height = 11,
  units = "cm"
)
ggsave(
  "2_figures/HordeumSynCom_structural_zeros_cutoff_order_composition_pooled.pdf",
  p_order_pooled,
  width = 22,
  height = 10,
  units = "cm"
)
