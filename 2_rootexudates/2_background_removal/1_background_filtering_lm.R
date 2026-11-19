# Seup ------------------------------------------------------------------------
library(data.table)

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading Lotus data ----------------------------------------------------------

metabolite_data_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_RootEx_mzmine_featurelist.csv"
)

design_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_RootEx_metadata.txt",
  drop = c(2, 4:7)
)

dir_name <- paste(
  "..",
  "1_data",
  "1_Lotus",
  "canopus_structure_summary.tsv",
  sep = "/"
)
annotation_Lj <- fread(dir_name, drop = 21)

# Keeping only columns containing the peak area
cols <- colnames(metabolite_data_Lj)
cols_to_keep <- cols[grepl(":area", cols)]
metabolite_data_Lj <- metabolite_data_Lj[, c("id", cols_to_keep), with = FALSE]

# # Replacing NAs with zeros
# metabolite_data_Lj <- metabolite_data_Lj[,
#   lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
#   id
# ]

# Constructing samples numbers from column names
sample_numbers <- unlist(lapply(
  strsplit(colnames(metabolite_data_Lj)[-1], "_"),
  function(x) x[1]
))
sample_numbers <- gsub("datafile:", "", sample_numbers)

## Fixing typo in feature table
sample_numbers[22] <- "53"
colnames(metabolite_data_Lj)[-1] <- sample_numbers

# Setting feature names in feature table
setnames(metabolite_data_Lj, "id", "Feature")
metabolite_data_Lj[, Feature := paste0("Feature", Feature)]

# Setting sample names in feature tables
colnames(metabolite_data_Lj)[-1] <- paste0(
  "Sample", colnames(metabolite_data_Lj)[-1]
)

# Setting sample names in metadata
design_Lj[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Lj <- design_Lj[ Sample_ID %in% colnames(metabolite_data_Lj)]

# Removing samples not in metadata from feature table
metabolite_data_Lj <- metabolite_data_Lj[,
  c("Feature", design_Lj$Sample_ID),
  with = FALSE
]

# Setting feature names in annotation table
setnames(annotation_Lj, "mappingFeatureId", "Feature")
setcolorder(annotation_Lj, "Feature")
annotation_Lj[, Feature := paste0("Feature", Feature)]

# Loading Hordeum data --------------------------------------------------------
metabolite_data_Hv <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_RootEx_mzmine_featurelist.csv"
)

design_Hv <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_RootEx_metadata.txt",
  drop = c(2, 4:7)
)

dir_name <- paste(
  "..",
  "1_data",
  "2_Hordeum",
  "canopus_structure_summary.tsv",
  sep = "/"
)
annotation_Hv <- fread(dir_name, drop = 21)

# Keeping only columns containing the peak area
cols <- colnames(metabolite_data_Hv)
cols_to_keep <- cols[grepl(":area", cols)]
metabolite_data_Hv <- metabolite_data_Hv[, c("id", cols_to_keep), with = FALSE]

# Replacing NAs with zeros
# metabolite_data_Hv <- metabolite_data_Hv[,
#   lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
#   id
# ]

# Constructing samples numbers from column names
sample_numbers <- unlist(lapply(
  strsplit(colnames(metabolite_data_Hv)[-1], "_"),
  function(x) x[4]
))
colnames(metabolite_data_Hv)[-1] <- sample_numbers

# Setting feature names in feature table
setnames(metabolite_data_Hv, "id", "Feature")
metabolite_data_Hv[,Feature:= paste0("Feature", Feature)]

# Setting sample names in feature tables
colnames(metabolite_data_Hv)[-1] <- paste0(
  "Sample", colnames(metabolite_data_Hv)[-1]
)

# Setting sample names in metadata
design_Hv[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Hv <- design_Hv[ Sample_ID %in% colnames(metabolite_data_Hv)]

# Removing samples not in metadata from feature table
metabolite_data_Hv <- metabolite_data_Hv[,
  c("Feature", design_Hv$Sample_ID),
  with = FALSE
]

# Setting feature names in annotation table
setnames(annotation_Hv, "mappingFeatureId", "Feature")
setcolorder(annotation_Hv, "Feature")
annotation_Hv[, Feature := paste0("Feature", Feature)]

# Subset the data --------------------------------------------------------------

# Keeping only Lotus features present more than 10% of replicates
present_features_idx <- rowMeans(!is.na(metabolite_data_Lj[,-1])) >= 0.1
present_features <- metabolite_data_Lj$Feature[present_features_idx]
metabolite_data_Lj <- metabolite_data_Lj[Feature %in% present_features]

# # Keeping only Lotus features with annotations
# metabolite_data_Lj <- metabolite_data_Lj[Feature %in% annotation_Lj$Feature]

# Keeping only Hordeum features present more than 10% of replicates
present_features_idx <- rowMeans(!is.na(metabolite_data_Lj[,-1])) >= 0.1
present_features <- metabolite_data_Hv$Feature[present_features_idx]
metabolite_data_Hv <- metabolite_data_Hv[Feature %in% present_features]

# # Keeping only Hordeum features with annotations
# metabolite_data_Hv <- metabolite_data_Hv[Feature %in% annotation_Hv$Feature]

# Filter based on sand control samples -----------------------------------------
# Check which features has statistically significantly higher intensity in at
# least one plant genotype compared to sand control samples, and keep only those
# features.

# Lotus
design_Lj[,
  Genotype := factor(
    Genotype,
    levels = c("control", "WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]
design_Lj <- design_Lj[order(Genotype)]

lm_over_features <- function(metabolite_data, design) {
  metabolite_trnsp <- transpose(
    metabolite_data,
    make.names = "Feature",
    keep.names = "Sample_ID"
  )

  empty_matrix <- function(metabolite_data) {
    p_val_mat <- matrix(NA, nrow = nrow(metabolite_data), ncol = 5)
    colnames(p_val_mat) <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
    rownames(p_val_mat) <- metabolite_data$Feature
    return(p_val_mat)
  }

  p_val_mat <- empty_matrix(metabolite_data)
  logfc_mat <- empty_matrix(metabolite_data)
  p_val_mat_pres <- empty_matrix(metabolite_data)
  logfc_mat_pres <- empty_matrix(metabolite_data)

  for (i in 1:(nrow(metabolite_data))) {
    cat(i, "\r")
    dt <- merge(
      design,
      metabolite_trnsp[, c(1, i + 1), with = FALSE],
      by = "Sample_ID"
    )
    colnames(dt)[3] <- "Intensity"
    dt[, present := is.na(Intensity)]

    # Presence-absence based test
    ll <- logistf(present~Genotype, data = dt)
    captured <- capture.output(llres <- summary(ll))
    p_val_mat_pres[i,] <- llres$prob[-1]
    logfc_mat_pres[i,] <- llres$coefficients[-1]

    # Intensity-based test
    present <- dt[, .(pres_prob = any(!is.na(Intensity))), Genotype]
    present <- present[order(Genotype)]
    # if (present$pres_prob[1] == FALSE) {
    #   p_val_mat[i, ] <- ifelse(present$pres_prob[-1], 0, 1)
    #   logfc_mat[i, ] <- ifelse(present$pres_prob[-1], 0, Inf)
    # } else if (sum(present$pres_prob) > 1) {
    if (present$pres_prob[1] == TRUE & sum(present$pres_prob) > 1) {
      l <- lm(log(Intensity) ~ Genotype, data = dt)
      lres <- coef(summary(l))

      p_vals <- lres[-1, "Pr(>|t|)"]
      names(p_vals) <- rownames(lres)[-1]
      names(p_vals) <- gsub("Genotype", "", names(p_vals))
      p_val_mat[i, names(p_vals)] <- p_vals

      logfc <- lres[-1, "Estimate"]
      names(logfc) <- rownames(lres)[-1]
      names(logfc) <- gsub("Genotype", "", names(p_vals))
      logfc_mat[i, names(logfc)] <- logfc
    }
  }
  
  out <- list(
    p_val = p_val_mat,
    logfc = logfc_mat,
    p_val_pres = p_val_mat_pres,
    logfc_pres = logfc_mat_pres
  )
  return(out)
}

lm_out <- lm_over_features(metabolite_data_Lj, design_Lj)
p_adj_mat_Lj <- apply(lm_out$p_val, 2, p.adjust, method = "fdr")
p_adj_mat_pres_Lj <- apply(lm_out$p_val_pres, 2, p.adjust, method = "fdr")

keep_mat_lm <- p_adj_mat_Lj < 0.05 & lm_out$logfc > 0
keep_mat_lm[is.na(keep_mat_lm)] <- FALSE
keep_mat_pres <- p_adj_mat_pres_Lj < 0.05 & lm_out$logfc_pres < 0

keep_mat <- keep_mat_pres|keep_mat_lm
# feature_keep_Lj_idx <- apply(p_adj_mat_Lj, 1, function(x) any(x < 0.05, na.rm = TRUE))
feature_keep_Lj_idx <- apply(keep_mat, 1, any , na.rm = TRUE)
feature_keep_Lj <- rownames(p_adj_mat_Lj)[feature_keep_Lj_idx]

# Hordeum
design_Hv[,
  Genotype := factor(
    Genotype,
    levels = c("control", "WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]
design_Hv <- design_Hv[order(Genotype)]

lm_out <- lm_over_features(metabolite_data_Hv, design_Hv)
p_adj_mat_Hv <- apply(lm_out$p_val, 2, p.adjust, method = "fdr")
p_adj_mat_pres_Hv <- apply(lm_out$p_val_pres, 2, p.adjust, method = "fdr")

keep_mat_lm <- p_adj_mat_Hv < 0.05 & lm_out$logfc > 0
keep_mat_lm[is.na(keep_mat_lm)] <- FALSE
keep_mat_pres <- p_adj_mat_pres_Hv < 0.05 & lm_out$logfc_pres < 0

keep_mat <- keep_mat_pres|keep_mat_lm
# feature_keep_Hv_idx <- apply(p_adj_mat_Hv, 1, function(x) any(x < 0.05, na.rm = TRUE))
feature_keep_Hv_idx <- apply(keep_mat, 1, any , na.rm = TRUE)
feature_keep_Hv <- rownames(p_adj_mat_Hv)[feature_keep_Hv_idx]

# p_val_mat_Hv <- lm_over_features(metabolite_data_Hv, design_Hv)
# p_adj_mat_Hv <- apply(p_val_mat_Hv, 2, p.adjust, method = "fdr")
# feature_keep_Hv_idx <- apply(p_adj_mat_Hv, 1, function(x) any(x < 0.05, na.rm = TRUE))
# feature_keep_Hv <- rownames(p_adj_mat_Hv)[feature_keep_Hv_idx]

# Removing background features -------------------------------------------------
# Lotus
## Removing control samples from feature tables
samples_keep <- design_Lj[Genotype != "control", Sample_ID]
metabolite_data_Lj <- metabolite_data_Lj[,
  c("Feature", samples_keep),
  with = FALSE
]

## Keeping only features with higher intensity in at least one plant genotype 
## than in control samples
metabolite_data_Lj <- metabolite_data_Lj[Feature %in% feature_keep_Lj]

# Hordeum
## Removing control samples from feature tables
samples_keep <- design_Hv[Genotype != "control", Sample_ID]
metabolite_data_Hv <- metabolite_data_Hv[,
  c("Feature", samples_keep),
  with = FALSE
]

## Keeping only features with higher intensity in at least one plant genotype 
## than in control samples
metabolite_data_Hv <- metabolite_data_Hv[Feature %in% feature_keep_Hv]

# Saving feature tables for use in other scripts
fwrite(metabolite_data_Lj, file = "1_tables/feature_table_Lotus_filtered.csv")
fwrite(metabolite_data_Hv, file = "1_tables/feature_table_Hordeum_filtered.csv")

# Saving feature tables for supplementary tables
metabolite_data_Lj[, Feature := gsub("Feature", "", Feature)]
metabolite_data_Hv[, Feature := gsub("Feature", "", Feature)]
setnames(metabolite_data_Lj, "Feature", "id")
setnames(metabolite_data_Hv, "Feature", "id")
colnames(metabolite_data_Lj) <- gsub("Sample", "", colnames(metabolite_data_Lj))
colnames(metabolite_data_Hv) <- gsub("Sample", "", colnames(metabolite_data_Hv))

fwrite(
  metabolite_data_Lj,
  file = "../10_suppl_tables/suppl_tab_filtered_feature_table_lotus.csv"
)
fwrite(
  metabolite_data_Hv,
  file = "../10_suppl_tables/suppl_tab_filtered_feature_table_Hordeum.csv"
)
