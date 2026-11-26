# Two-step (hurdle) analysis of genotype effects and comparison with Tobit ----
#
# For every feature, the genotype effect (each mutant vs. WT) is assessed with
# a two-step model instead of the Tobit model used in 0_genotype_effects.R:
#
#   Step 1 (presence/absence): Firth-penalised logistic regression
#            1{intensity > 0} ~ Genotype
#          Firth's penalty is used since complete separation (a feature absent
#          from all replicates of a genotype) is common.
#
#   Step 2 (intensity | present): ordinary linear regression
#            log(intensity) ~ Genotype,   using only non-zero observations.
#          Genotypes without any non-zero observations are not estimable in
#          this step and the step is skipped completely if WT (the reference)
#          has no non-zero observations or if no residual degrees of freedom
#          remain.
#
# The two steps form the two factors of the hurdle likelihood, so the two Wald/
# LR statistics are asymptotically independent. For each contrast they are
# combined into one test with 2 degrees of freedom
#   X2 = qchisq(1 - p_presence, 1) + qchisq(1 - p_intensity, 1) ~ chi2(2)
# (two-part test, cf. Lachenbruch 2002). If neither WT nor the mutant contains
# zeros, the presence/absence step is uninformative and only the intensity test
# (1 df) is used. If the intensity step is not estimable for a contrast, only
# the presence/absence test (1 df) is used, which mirrors the fallback to
# logistf in the Tobit analysis.
#
# The results of the two-step approach are then compared to those of the Tobit
# model (Tobit_model2) fitted exactly as in 0_genotype_effects.R.

# Setup -----------------------------------------------------------------------
library(data.table)
library(VGAM)
library(logistf)
library(ggplot2)

# Setting working directory to source file location
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

dir.create("1_tables", showWarnings = FALSE)
dir.create("2_figures", showWarnings = FALSE)

# Loading Tobit regression function
source("Tobit_function.R")
genotypes <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
mutants <- genotypes[-1]
alpha <- 0.05

## Loading Lotus data ----------------------------------------------------------
# Lotus feature table filtered for background features
metabolite_data_Lj <- fread(
  "../2_background_removal/1_tables/feature_table_Lotus_filtered.csv"
)
design_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_RootEx_metadata.txt",
  drop = 4:7
)
annotation_Lj <- fread(
  "../1_data/1_Lotus/canopus_structure_summary.tsv",
  drop = 21
)

# Setting sample names in metadata
design_Lj[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Lj <- design_Lj[Sample_ID %in% colnames(metabolite_data_Lj)]

# Setting feature names in annotation table
setnames(annotation_Lj, "mappingFeatureId", "Feature")
setcolorder(annotation_Lj, "Feature")
annotation_Lj[, Feature := paste0("Feature", Feature)]

# Removing control samples from metadata and feature table
non_control_samples <- design_Lj[Genotype != "control", Sample_ID]
design_Lj <- design_Lj[Sample_ID %in% non_control_samples]
metabolite_data_Lj <- metabolite_data_Lj[,
  c("Feature", non_control_samples),
  with = FALSE
]

## Loading Hordeum data --------------------------------------------------------
metabolite_data_Hv <- fread(
  "../2_background_removal/1_tables/feature_table_Hordeum_filtered.csv"
)
design_Hv <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_RootEx_metadata.txt",
  drop = c(2, 4:7)
)
annotation_Hv <- fread(
  "../1_data/2_Hordeum/canopus_structure_summary.tsv",
  drop = 21
)

# Setting sample names in metadata
design_Hv[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Hv <- design_Hv[Sample_ID %in% colnames(metabolite_data_Hv)]

# Setting feature names in annotation table
setnames(annotation_Hv, "mappingFeatureId", "Feature")
setcolorder(annotation_Hv, "Feature")
annotation_Hv[, Feature := paste0("Feature", Feature)]

# Functions -------------------------------------------------------------------
# Converts the data tables to the data frames used by the model functions and
# keeps features present in at least 3 samples (as in 0_genotype_effects.R)
prepare_data <- function(metabolite_data, design) {
  design[, Genotype := factor(Genotype, levels = genotypes)]
  metabolites <- data.frame(
    metabolite_data[, -1],
    row.names = metabolite_data$Feature
  )
  metabolites <- metabolites[, design$Sample_ID]
  design_df <- data.frame(design[, -1], row.names = design$Sample_ID)

  feat_keep_idx <- rowSums(metabolites != 0) >= 3
  metabolites <- metabolites[feat_keep_idx, ]

  list(metabolites = metabolites, design = design_df)
}

# Two-step model: Firth logistic regression for presence/absence and linear
# regression on log-intensities of the non-zero observations
two_step_model <- function(metabolites, meta_data, ref = "WT") {
  p <- nrow(metabolites)
  contrasts <- paste0("Genotype", setdiff(levels(meta_data$Genotype), ref))

  empty_matrix <- function() {
    matrix(
      NA_real_,
      nrow = p,
      ncol = length(contrasts),
      dimnames = list(rownames(metabolites), contrasts)
    )
  }

  beta_pa <- empty_matrix()
  p_pa <- empty_matrix()
  beta_int <- empty_matrix()
  se_int <- empty_matrix()
  p_int <- empty_matrix()
  sigma_int <- rep(NA_real_, p)
  # TRUE if neither the reference nor the mutant group contains zeros, in
  # which case the presence/absence step carries no information
  no_zeros <- matrix(
    FALSE,
    nrow = p,
    ncol = length(contrasts),
    dimnames = list(rownames(metabolites), contrasts)
  )

  for (i in seq_len(p)) {
    cat(i, "\r")
    y <- as.numeric(metabolites[i, rownames(meta_data)])
    df <- data.frame(
      Genotype = meta_data$Genotype,
      pres = as.integer(y > 0),
      log_int = ifelse(y > 0, log(y), NA_real_)
    )

    all_present <- tapply(df$pres, df$Genotype, function(x) all(x == 1))
    no_zeros[i, ] <- all_present[ref] & all_present[gsub("Genotype", "", contrasts)]

    # Step 1: presence/absence (Firth's penalised logistic regression)
    fit_pa <- logistf(pres ~ Genotype, data = df)
    beta_pa[i, ] <- fit_pa$coefficients[contrasts]
    p_pa[i, ] <- fit_pa$prob[contrasts]

    # Step 2: log-intensity of non-zero observations
    df_pos <- df[df$pres == 1, ]
    df_pos$Genotype <- droplevels(df_pos$Genotype)
    if (ref %in% df_pos$Genotype && nlevels(df_pos$Genotype) > 1) {
      fit_int <- lm(log_int ~ Genotype, data = df_pos)
      if (fit_int$df.residual > 0) {
        cf <- coef(summary(fit_int))
        idx <- intersect(rownames(cf), contrasts)
        beta_int[i, idx] <- cf[idx, "Estimate"]
        se_int[i, idx] <- cf[idx, "Std. Error"]
        p_int[i, idx] <- cf[idx, "Pr(>|t|)"]
        sigma_int[i] <- summary(fit_int)$sigma
      }
    }
  }

  # Combined two-part test: 2 df when both steps are informative, 1 df
  # (intensity only) when the two groups contain no zeros, and 1 df
  # (presence/absence only) when the intensity step is not estimable
  chi_pa <- qchisq(p_pa, df = 1, lower.tail = FALSE)
  chi_int <- qchisq(p_int, df = 1, lower.tail = FALSE)
  int_only <- !is.na(chi_int) & no_zeros
  both <- !is.na(chi_int) & !no_zeros
  p_comb <- p_pa
  p_comb[both] <- pchisq(chi_pa[both] + chi_int[both], 2, lower.tail = FALSE)
  p_comb[int_only] <- p_int[int_only]
  chi_pa[int_only] <- 0

  # Direction of the combined effect: sign of the step contributing most to
  # the combined statistic
  dir_comb <- sign(beta_pa)
  int_dominates <- !is.na(chi_int) & chi_int > chi_pa
  dir_comb[int_dominates] <- sign(beta_int[int_dominates])
  discordant <- both & sign(beta_pa) != sign(beta_int)

  list(
    beta_pa = beta_pa,
    p_pa = p_pa,
    beta_int = beta_int,
    se_int = se_int,
    p_int = p_int,
    sigma_int = sigma_int,
    p_comb = p_comb,
    dir_comb = dir_comb,
    discordant = discordant
  )
}

# FDR adjustment across all features and contrasts jointly (as in
# 0_genotype_effects.R)
fdr_matrix <- function(p_mat) {
  p_adj <- p_mat
  p_adj[] <- p.adjust(p_mat, "fdr")
  p_adj
}

# Significance calls coded -1 / 0 / 1
call_matrix <- function(p_adj, direction) {
  out <- ifelse(!is.na(p_adj) & p_adj < alpha, sign(direction), 0)
  out[is.na(out)] <- 0
  out
}

# Long-format table of matrices (features x contrasts)
to_long <- function(mat, value_name) {
  dt <- data.table(Feature = rownames(mat), mat)
  setnames(dt, gsub("Genotype", "", names(dt)))
  dt <- melt(
    dt,
    id.vars = "Feature",
    variable.name = "Genotype",
    value.name = value_name
  )
  dt[, Genotype := as.character(Genotype)]
  dt
}

# Runs the Tobit and two-step analyses for one species and collects results
analyse_species <- function(metabolite_data, design, annotation, species) {
  dat <- prepare_data(metabolite_data, design)
  metabolites <- dat$metabolites
  design_df <- dat$design

  ### Tobit model (identical to 0_genotype_effects.R) ----
  set.seed(1793694715)
  tobit_fit <- Tobit_model2(
    metabolites = metabolites,
    meta_data = design_df,
    formula = ~Genotype,
    ref = "WT",
    test = "two-sided",
    presence_absence = TRUE
  )
  contrasts <- paste0("Genotype", mutants)

  tobit_p_raw <- tobit_fit$p_vals[, contrasts]
  tobit_fallback <- is.na(tobit_p_raw)
  tobit_p <- tobit_p_raw
  tobit_p[tobit_fallback] <- tobit_fit$p_vals_pa[, contrasts][tobit_fallback]
  tobit_p_adj <- fdr_matrix(tobit_p)

  tobit_lfc <- as.matrix(tobit_fit$res[, contrasts])
  tobit_dir <- tobit_lfc
  tobit_dir[tobit_fallback] <- tobit_fit$beta_est_pa[, contrasts][tobit_fallback]
  tobit_call <- call_matrix(tobit_p_adj, tobit_dir)

  ### Two-step model ----
  ts_fit <- two_step_model(metabolites, design_df, ref = "WT")
  ts_p_adj <- fdr_matrix(ts_fit$p_comb)
  ts_call <- call_matrix(ts_p_adj, ts_fit$dir_comb)

  # Separate calls for each step (each adjusted on its own)
  pa_p_adj <- fdr_matrix(ts_fit$p_pa)
  int_p_adj <- fdr_matrix(ts_fit$p_int)
  pa_call <- call_matrix(pa_p_adj, ts_fit$beta_pa)
  int_call <- call_matrix(int_p_adj, ts_fit$beta_int)

  ### Two-step result table (same layout as the Tobit result table) ----
  res <- data.table(Feature = rownames(ts_call), ts_call)
  setnames(res, gsub("Genotype", "", names(res)))
  res <- res[symrk != 0 | ccamk != 0 | nsp1 != 0 | nsp2 != 0]

  feat_tab <- metabolites
  colnames(feat_tab) <- paste(
    design_df$Genotype,
    rownames(design_df),
    sep = "_"
  )
  feat_tab <- data.table(Feature = rownames(feat_tab), feat_tab)
  res <- merge(res, feat_tab, by = "Feature")
  res <- merge(res, annotation, by = "Feature", all.x = TRUE)
  fwrite(
    res,
    paste0("1_tables/", species, "_metabolite_test_results_two_step.csv")
  )

  ### Long-format comparison table ----
  prop_zero <- data.table(
    Feature = rownames(metabolites),
    prop_zero = rowMeans(metabolites == 0)
  )
  pieces <- list(
    to_long(tobit_lfc, "tobit_lfc"),
    to_long(tobit_p, "tobit_p"),
    to_long(tobit_p_adj, "tobit_p_adj"),
    to_long(tobit_fallback * 1, "tobit_fallback_pa"),
    to_long(tobit_call, "tobit_call"),
    to_long(ts_fit$beta_pa, "pa_logodds"),
    to_long(ts_fit$p_pa, "pa_p"),
    to_long(pa_p_adj, "pa_p_adj"),
    to_long(pa_call, "pa_call"),
    to_long(ts_fit$beta_int, "int_lfc"),
    to_long(ts_fit$p_int, "int_p"),
    to_long(int_p_adj, "int_p_adj"),
    to_long(int_call, "int_call"),
    to_long(ts_fit$p_comb, "two_step_p"),
    to_long(ts_p_adj, "two_step_p_adj"),
    to_long(ts_call, "two_step_call"),
    to_long(ts_fit$discordant * 1, "pa_int_discordant")
  )
  comp <- Reduce(
    function(a, b) merge(a, b, by = c("Feature", "Genotype")),
    pieces
  )
  comp <- merge(comp, prop_zero, by = "Feature")
  comp[, Species := species]
  comp[, tobit_fallback_pa := as.logical(tobit_fallback_pa)]
  comp[, pa_int_discordant := as.logical(pa_int_discordant)]
  setcolorder(comp, c("Species", "Feature", "Genotype", "prop_zero"))
  fwrite(
    comp,
    paste0("1_tables/", species, "_tobit_vs_two_step_comparison.csv")
  )

  comp
}

# Running analyses ------------------------------------------------------------
comp_Lj <- analyse_species(
  metabolite_data_Lj,
  design_Lj,
  annotation_Lj,
  "Lotus"
)
comp_Hv <- analyse_species(
  metabolite_data_Hv,
  design_Hv,
  annotation_Hv,
  "Hordeum"
)
comp <- rbind(comp_Lj, comp_Hv)
comp[, Genotype := factor(Genotype, levels = mutants)]
comp[, Species := factor(Species, levels = c("Lotus", "Hordeum"))]

# Comparison ------------------------------------------------------------------
## Numbers of significant features per method ---------------------------------
n_sig <- comp[,
  .(
    n_features = .N,
    tobit = sum(tobit_call != 0),
    tobit_from_pa_fallback = sum(tobit_call != 0 & tobit_fallback_pa),
    two_step = sum(two_step_call != 0),
    two_step_presence = sum(pa_call != 0),
    two_step_intensity = sum(int_call != 0),
    both_tobit_two_step = sum(tobit_call != 0 & two_step_call != 0),
    tobit_only = sum(tobit_call != 0 & two_step_call == 0),
    two_step_only = sum(tobit_call == 0 & two_step_call != 0),
    opposite_direction = sum(tobit_call * two_step_call == -1),
    agreement = mean(tobit_call == two_step_call),
    kappa = {
      # Cohen's kappa for the -1/0/1 calls
      lv <- c(-1, 0, 1)
      tab <- table(factor(tobit_call, lv), factor(two_step_call, lv))
      po <- sum(diag(tab)) / sum(tab)
      pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
      (po - pe) / (1 - pe)
    }
  ),
  by = .(Species, Genotype)
]
setorder(n_sig, Species, Genotype)
fwrite(n_sig, "1_tables/tobit_vs_two_step_summary.csv")
print(n_sig)

## Cross-tabulation of calls --------------------------------------------------
crosstab <- comp[,
  .N,
  by = .(Species, Genotype, tobit_call, two_step_call)
]
setorder(crosstab, Species, Genotype, tobit_call, two_step_call)
fwrite(crosstab, "1_tables/tobit_vs_two_step_crosstab.csv")

## Agreement of effect estimates and p-values ---------------------------------
# Tobit LFC vs. LFC from the intensity step for contrasts where both are
# estimated from the continuous part (i.e. no presence/absence fallback)
est_agreement <- comp[
  !tobit_fallback_pa & !is.na(int_lfc),
  .(
    n = .N,
    cor_lfc_pearson = cor(tobit_lfc, int_lfc),
    median_tobit_minus_int = median(tobit_lfc - int_lfc),
    cor_logp_spearman = cor(
      -log10(tobit_p),
      -log10(two_step_p),
      method = "spearman"
    )
  ),
  by = .(Species, Genotype)
]
setorder(est_agreement, Species, Genotype)
fwrite(est_agreement, "1_tables/tobit_vs_two_step_estimates.csv")
print(est_agreement)

# Same comparison, split by amount of zeros in the feature
comp[,
  zero_class := cut(
    prop_zero,
    breaks = c(-Inf, 0, 0.25, 0.5, 1),
    labels = c("no zeros", "<= 25% zeros", "25-50% zeros", "> 50% zeros")
  )
]
est_by_zero <- comp[
  !tobit_fallback_pa & !is.na(int_lfc),
  .(
    n = .N,
    cor_lfc_pearson = cor(tobit_lfc, int_lfc),
    median_abs_diff_lfc = median(abs(tobit_lfc - int_lfc)),
    n_sig_tobit = sum(tobit_call != 0),
    n_sig_two_step = sum(two_step_call != 0),
    n_sig_both = sum(tobit_call != 0 & two_step_call != 0)
  ),
  by = .(Species, zero_class)
]
setorder(est_by_zero, Species, zero_class)
fwrite(est_by_zero, "1_tables/tobit_vs_two_step_by_zeros.csv")
print(est_by_zero)

## Figures --------------------------------------------------------------------
# Log-fold-changes: Tobit vs. intensity step
lfc_plot <- ggplot(
  comp[!tobit_fallback_pa & !is.na(int_lfc)],
  aes(x = int_lfc, y = tobit_lfc, colour = prop_zero)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.8, alpha = 0.7) +
  facet_grid(Species ~ Genotype) +
  scale_colour_viridis_c(name = "Proportion\nof zeros") +
  labs(
    x = "Log-fold-change, two-step model (intensity step)",
    y = "Log-fold-change, Tobit model"
  ) +
  theme_bw()

ggsave(
  "2_figures/tobit_vs_two_step_lfc.pdf",
  lfc_plot,
  width = 220,
  height = 120,
  units = "mm"
)

# P-values: Tobit vs. combined two-step test
comp[,
  agreement := factor(
    fcase(
      tobit_call != 0 & two_step_call != 0, "Both",
      tobit_call != 0, "Tobit only",
      two_step_call != 0, "Two-step only",
      default = "Neither"
    ),
    levels = c("Both", "Tobit only", "Two-step only", "Neither")
  )
]
p_plot <- ggplot(
  comp,
  aes(x = -log10(two_step_p), y = -log10(tobit_p), colour = agreement)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.8, alpha = 0.7) +
  facet_grid(Species ~ Genotype) +
  scale_colour_manual(
    values = c(
      "Both" = "#1a80bb",
      "Tobit only" = "#a00000",
      "Two-step only" = "#ea801c",
      "Neither" = "grey70"
    ),
    name = "Significant\n(FDR < 0.05)"
  ) +
  labs(
    x = expression(-log[10] ~ "p, two-step model (combined)"),
    y = expression(-log[10] ~ "p, Tobit model")
  ) +
  theme_bw()

ggsave(
  "2_figures/tobit_vs_two_step_pvalues.pdf",
  p_plot,
  width = 220,
  height = 120,
  units = "mm"
)

# Number of significant features per method
n_sig_long <- melt(
  n_sig,
  id.vars = c("Species", "Genotype"),
  measure.vars = c("tobit", "two_step", "two_step_presence", "two_step_intensity"),
  variable.name = "Method",
  value.name = "n"
)
n_sig_long[,
  Method := factor(
    Method,
    levels = c("tobit", "two_step", "two_step_presence", "two_step_intensity"),
    labels = c(
      "Tobit",
      "Two-step (combined)",
      "Two-step (presence)",
      "Two-step (intensity)"
    )
  )
]
n_sig_plot <- ggplot(n_sig_long, aes(x = Genotype, y = n, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~Species, scales = "free_y") +
  scale_fill_manual(values = c("#a00000", "#1a80bb", "#8cc5e3", "#b8b8b8")) +
  labs(y = "Number of significant features (FDR < 0.05)", x = NULL) +
  theme_bw()

ggsave(
  "2_figures/tobit_vs_two_step_n_sig.pdf",
  n_sig_plot,
  width = 200,
  height = 90,
  units = "mm"
)
