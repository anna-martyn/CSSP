# Sensitivity analysis of the structural zero detection used in the
# genotype-specific differential abundance analysis (0_Maaslin2_DAA_Lotus.R).
# Only the structural zero part is run, for a grid of read and sample cutoffs
# The cutoff used in the main analysis is min_reads = 20, min_present_reps = 3.
# As in the main analysis tables, only ASVs that pass MaAsLin2's feature
# filtering (prevalence > 10% of samples and non-zero variance) are kept.
# For each cutoff the script summarises:
#   1. the number of ASVs DA by structural zero,
#   2. their cumulative mean relative abundance (RA) in WT (as in the barplot
#      of the main analysis) and
#   3. their order-level taxonomic composition.

# Setup ------------------------------------------------------------------------
# Clean up
options(warn = -1)
rm(list = ls())

# Loading packages
pkg <- c("data.table", "ggplot2")

for(pk in pkg){
  library(pk, character.only = T)
}

# Setting directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load file with function for structural zero analysis
source("../0_files/Structural_zeros.R")

# Loading data
asv_table <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
colnames(asv_table)[1] <- "ASVid"

design <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt",
  drop = c(5, 7, 8)
)

taxonomy <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv"
)

# Cleaning up taxonomy
setnames(taxonomy, "Feature ID", "ASVid")
taxa_levels <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)
taxonomy[,c(taxa_levels):= tstrsplit(Taxon, "; ", fill = "u__Unknown")]
taxonomy[,c(taxa_levels):= lapply(.SD, substr, 4, 1000), .SDcols = taxa_levels]
taxonomy[Kingdom == "ssigned", Kingdom:= "Unassigned"]
taxonomy[,Taxon:=NULL]

# Order colours and orders to display (same as in the main analysis figures).
# Orders not among these are grouped as 'Other'
colors_orders <- fread("../../../../0_files/Bacterial_order_colors.csv")
combined_top_orders <- readRDS("../0_files/Orders_to_display.rds")
taxonomy[, Order_display := fifelse(
  Order %in% combined_top_orders,
  Order,
  "Other"
)]

# Merging ASV table and metadata
asv_table_with_metadata <- transpose(
  asv_table,
  make.names = "ASVid",
  keep.names = "SampleID"
)
asv_table_with_metadata <- merge(
  design,
  asv_table_with_metadata,
  by = "SampleID"
)

# Vector of mutants
mutants <- c("symrk", "ccamk", "nsp1", "nsp2")

# Setting reference levels
asv_table_with_metadata[,
  Genotype := factor(
    Genotype,
    levels = c("WT", mutants)
  )
]

opt <- as.matrix(
  expand.grid(
    Compartment = c("Root", "Rhizosphere"),
    Soil = c("NPK", "PK", "UF")
  )
)

# Cutoff grid
cutoffs <- CJ(min_reads = seq(0, 40, 10), min_present_reps = 2:4)
baseline <- list(min_reads = 20, min_present_reps = 3)

# Subsetting data for each compartment-soil combination ------------------------
asv_tables_subsets <- apply(
  opt, 1,
  function(x){
    # Subsetting data for compartment-soil combination
    data_subset <- asv_table_with_metadata[Compartment == x[1] & Soil == x[2]]

    # Extracting metadata and dropping unused factor levels
    col_design <- c("Plant", "Genotype", "Compartment", "Soil")
    design <- data.frame(
      data_subset[, ..col_design],
      row.names = data_subset$SampleID
    )
    design$Genotype <- droplevels(design$Genotype)

    # Extracting and transposing ASV table
    asv_table <- t(data_subset[,-(1:5)])
    colnames(asv_table) <- data_subset$SampleID

    # Library sizes needed by Structural_zeros3
    design$library_size <- colSums(asv_table)

    # ASVs kept by MaAsLin2's filtering (default min_abundance = 0,
    # min_prevalence = 0.1, min_variance = 0). Structural zero ASVs outside
    # this set are not included in the main analysis tables.
    n_present <- rowSums(asv_table > 0)
    maaslin_features <- rownames(asv_table)[
      n_present > ncol(asv_table) * 0.1 & apply(asv_table, 1, var) > 0
    ]

    # Mean relative abundance of each ASV within each genotype
    ra_table <- t(t(asv_table) / colSums(asv_table))
    mean_ra <- sapply(
      levels(design$Genotype),
      function(g) rowMeans(ra_table[, design$Genotype == g, drop = FALSE])
    )

    out <- list(
      asv_table = asv_table,
      design = design,
      maaslin_features = maaslin_features,
      mean_ra = mean_ra
    )
    return(out)
  }
)

# Structural zero analysis for each subset and cutoff --------------------------
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
      asvs <- intersect(
        da_struc_zero$struc_zero_DA[[current_mutant]],
        X$maaslin_features
      )
      if(length(asvs) == 0){
        return(NULL)
      }
      zero_tab <- da_struc_zero$struc_zero_table[
        asvs,
        c(current_mutant, "WT"),
        drop = F
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
  res_cutoffs[, ":="(
    Compartment = X$design$Compartment[1],
    Soil = X$design$Soil[1]
  )]
  return(res_cutoffs)
})

struc_zero_res <- rbindlist(struc_zero_results)

# Adding taxonomy
struc_zero_res <- merge(
  struc_zero_res,
  taxonomy[, .(ASVid, Phylum, Class, Order, Family, Genus, Order_display)],
  by = "ASVid",
  all.x = TRUE
)
struc_zero_res[is.na(Order_display), Order_display := "Other"]
setcolorder(
  struc_zero_res,
  c("min_reads", "min_present_reps", "Soil", "Compartment", "Genotype",
    "ASVid", "Direction", "WT_mean_RA", "Mutant_mean_RA")
)

# Summarising number of DA ASVs per cutoff -------------------------------------
# Full grid so that combinations without any DA ASVs appear with zero counts
summary_grid <- CJ(
  min_reads = cutoffs$min_reads,
  min_present_reps = cutoffs$min_present_reps,
  Soil = c("NPK", "PK", "UF"),
  Compartment = c("Rhizosphere", "Root"),
  Genotype = mutants,
  unique = TRUE
)

# Cumulative mean RA in WT of the DA ASVs (as in the main analysis barplot).
# ASVs enriched in the mutant are absent from WT and contribute zero, so their
# cumulative mean RA in the mutant is given as well.
counts <- struc_zero_res[,
  .(
    DA_ASVs = .N,
    Enriched = sum(Direction == "Enriched"),
    Depleted = sum(Direction == "Depleted"),
    WT_cum_RA = sum(WT_mean_RA),
    Mutant_cum_RA_enriched = sum(Mutant_mean_RA[Direction == "Enriched"])
  ),
  by = .(min_reads, min_present_reps, Soil, Compartment, Genotype)
]

# Overlap with ASVs detected using the cutoff of the main analysis
baseline_res <- struc_zero_res[
  min_reads == baseline$min_reads &
    min_present_reps == baseline$min_present_reps,
  .(Soil, Compartment, Genotype, ASVid, in_baseline = TRUE)
]
overlap <- merge(
  struc_zero_res,
  baseline_res,
  by = c("Soil", "Compartment", "Genotype", "ASVid"),
  all.x = TRUE
)[,
  .(Shared_with_baseline = sum(!is.na(in_baseline))),
  by = .(min_reads, min_present_reps, Soil, Compartment, Genotype)
]
baseline_n <- baseline_res[, .(Baseline_DA_ASVs = .N), by = .(Soil, Compartment, Genotype)]

struc_zero_summary <- merge(
  summary_grid,
  counts,
  by = c("min_reads", "min_present_reps", "Soil", "Compartment", "Genotype"),
  all.x = TRUE
)
struc_zero_summary <- merge(
  struc_zero_summary,
  overlap,
  by = c("min_reads", "min_present_reps", "Soil", "Compartment", "Genotype"),
  all.x = TRUE
)
struc_zero_summary <- merge(
  struc_zero_summary,
  baseline_n,
  by = c("Soil", "Compartment", "Genotype"),
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
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Soil = factor(Soil, levels = c("NPK", "PK", "UF"))
)]
setcolorder(
  struc_zero_summary,
  c("min_reads", "min_present_reps", "Soil", "Compartment", "Genotype",
    count_cols, ra_cols)
)
struc_zero_summary <- struc_zero_summary[
  order(min_reads, min_present_reps, Soil, Compartment, Genotype)
]

# Order-level composition of DA ASVs per cutoff --------------------------------
# Number of DA ASVs and their cumulative mean RA in WT per order, along with
# the proportion of DA ASVs (within cutoff, soil, compartment and genotype)
# belonging to each order
order_comp <- struc_zero_res[,
  .(
    DA_ASVs = .N,
    Enriched = sum(Direction == "Enriched"),
    Depleted = sum(Direction == "Depleted"),
    WT_cum_RA = sum(WT_mean_RA)
  ),
  by = .(min_reads, min_present_reps, Soil, Compartment, Genotype,
         Order = Order_display)
]
order_comp[,
  Prop_DA_ASVs := DA_ASVs / sum(DA_ASVs),
  by = .(min_reads, min_present_reps, Soil, Compartment, Genotype)
]
order_comp[, ":="(
  Genotype = factor(Genotype, levels = mutants),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Soil = factor(Soil, levels = c("NPK", "PK", "UF")),
  Order = factor(Order, levels = colors_orders$Order)
)]
order_comp <- order_comp[
  order(min_reads, min_present_reps, Soil, Compartment, Genotype, Order)
]

# Same composition pooled across soils, compartments and genotypes. The
# cumulative mean RA in WT is averaged over the soil-compartment-genotype
# combinations
n_combinations <- nrow(opt) * length(mutants)
order_comp_pooled <- struc_zero_res[,
  .(DA_ASVs = .N, WT_cum_RA = sum(WT_mean_RA) / n_combinations),
  by = .(min_reads, min_present_reps, Order = Order_display)
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
  "2_tables/Lotus_structural_zeros_cutoff_overview.csv"
)
fwrite(
  struc_zero_res[
    order(min_reads, min_present_reps, Soil, Compartment, Genotype, ASVid),
    !"Order_display"
  ],
  "2_tables/Lotus_structural_zeros_cutoff_results.csv"
)
fwrite(
  order_comp[, !"Cutoff"],
  "2_tables/Lotus_structural_zeros_cutoff_order_composition.csv"
)
fwrite(
  order_comp_pooled[, !"Cutoff"],
  "2_tables/Lotus_structural_zeros_cutoff_order_composition_pooled.csv"
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
  facet_grid(Compartment + Soil ~ Genotype, scales = "free_y") +
  scale_x_continuous(breaks = unique(cutoffs$min_reads)) +
  labs(
    x = "Minimum number of reads",
    y = "Cumulative mean RA in WT of ASVs DA by structural zero",
    colour = "Minimum number\nof samples",
    title = "Lotus"
  ) +
  theme_cutoff

# Order composition (proportion of DA ASVs) per soil, compartment and genotype
p_order <- ggplot(
  data = order_comp,
  mapping = aes(x = Cutoff, y = Prop_DA_ASVs, fill = Order)
) +
  geom_col(width = 0.9) +
  facet_grid(Compartment + Soil ~ Genotype) +
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
    title = "Lotus"
  ) +
  theme_cutoff +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

# Order composition pooled across soils, compartments and genotypes, shown both
# as proportion of DA ASVs and as cumulative mean RA in WT
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
    title = "Lotus (pooled across soils, compartments and genotypes)"
  ) +
  theme_cutoff +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
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
  facet_grid(Compartment + Soil ~ Genotype, scales = "free_y") +
  scale_x_continuous(breaks = unique(cutoffs$min_reads)) +
  labs(
    x = "Minimum number of reads",
    y = "ASVs DA by structural zero",
    colour = "Minimum number\nof samples",
    title = "Lotus"
  ) +
  theme_cutoff

# Saving plots -----------------------------------------------------------------
ggsave(
  "../4_figures/Lotus_structural_zeros_cutoff_sensitivity.pdf",
  p_cutoffs,
  width = 18,
  height = 18,
  units = "cm"
)
ggsave(
  "../4_figures/Lotus_structural_zeros_cutoff_WT_RA.pdf",
  p_ra,
  width = 18,
  height = 18,
  units = "cm"
)
ggsave(
  "../4_figures/Lotus_structural_zeros_cutoff_order_composition.pdf",
  p_order,
  width = 24,
  height = 22,
  units = "cm"
)
ggsave(
  "../4_figures/Lotus_structural_zeros_cutoff_order_composition_pooled.pdf",
  p_order_pooled,
  width = 22,
  height = 10,
  units = "cm"
)
