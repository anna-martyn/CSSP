# Seup -------------------------------------------------------------------------
library(data.table)
library(logistf)

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Loading Lotus data ----------------------------------------------------------
# Lotus feature table filtered for background features
dir_name <- paste(
  "..",
  "2_background_removal",
  "1_tables",
  "feature_table_Lotus_filtered.csv",
  sep = "/"
)
metabolite_data_Lj <- fread(dir_name)

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

# Setting sample names in metadata
design_Lj[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Lj <- design_Lj[ Sample_ID %in% colnames(metabolite_data_Lj)]

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

for(i in 2:ncol(metabolite_data_Lj)){
  metabolite_data_Lj[[i]] <- ifelse(metabolite_data_Lj[[i]] == 0, NA, metabolite_data_Lj[[i]])
}

## Loading Hordeum data --------------------------------------------------------
dir_name <- paste(
  "..",
  "2_background_removal",
  "1_tables",
  "feature_table_Hordeum_filtered.csv",
  sep = "/"
)
metabolite_data_Hv <- fread(dir_name)
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

# Setting sample names in metadata
design_Hv[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Hv <- design_Hv[ Sample_ID %in% colnames(metabolite_data_Hv)]

# Setting feature names in annotation table
setnames(annotation_Hv, "mappingFeatureId", "Feature")
setcolorder(annotation_Hv, "Feature")
annotation_Hv[, Feature := paste0("Feature", Feature)]

for(i in 2:ncol(metabolite_data_Hv)){
  metabolite_data_Hv[[i]] <- ifelse(metabolite_data_Hv[[i]] == 0, NA, metabolite_data_Hv[[i]])
}

# Genotype effects ------------------------------------------------------------
lm_over_features <- function(metabolite_data, design) {
  metabolite_trnsp <- transpose(
    metabolite_data,
    make.names = "Feature",
    keep.names = "Sample_ID"
  )

  empty_matrix <- function(metabolite_data) {
    p_val_mat <- matrix(NA, nrow = nrow(metabolite_data), ncol = 4)
    colnames(p_val_mat) <- c("symrk", "ccamk", "nsp1", "nsp2")
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

## Lotus ----------------------------------------------------------------------
# Setting factor levels
design_Lj[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

lm_res_Lj <- lm_over_features(metabolite_data_Lj, design_Lj)
p_adj_Lj <- apply(lm_res_Lj$p_val, 2, p.adjust, method = "fdr")
p_adj_pres_Lj <- apply(lm_res_Lj$p_val_pres, 2, p.adjust, method = "fdr")

# p_adj_Lj[highlighted_features_Lj,]
# p_adj_pres_Lj[highlighted_features_Lj,]

# Setting up p-value matrix
p_adj_Lj <- data.table(
  Feature = rownames(test_symrk_Lj),
  symrk_p_adj = test_symrk_Lj$p_adj,
  ccamk_p_adj = test_ccamk_Lj$p_adj,
  nsp1_p_adj = test_nsp1_Lj$p_adj,
  nsp2_p_adj = test_nsp2_Lj$p_adj
)

# Setting up log-fold-change matrix
lfc_Lj <- data.table(Feature = rownames(tobit_fit_Lj$res), tobit_fit_Lj$res[,2:5])
colnames(lfc_Lj) <- gsub("Genotype", "Lfc_", colnames(lfc_Lj))

# Saving above results
fwrite(p_adj_Lj, "1_tables/p_adj_Lj.csv")
fwrite(lfc_Lj, "1_tables/lfc_Lj.csv")

# Saving results
res_Lj <- data.table(
  Feature = rownames(test_symrk_Lj),
  ifelse(tobit_fit_Lj$res[,2:5] > 0, 1, -1)
)
colnames(res_Lj) <- gsub("Genotype", "", colnames(res_Lj))

res_Lj[, ":="(
  symrk = ifelse(
    p_adj_Lj$symrk_p_adj < 0.05 , symrk , 0
  ),
  ccamk = ifelse(
    p_adj_Lj$ccamk_p_adj < 0.05 , ccamk , 0
  ),
  nsp1 = ifelse(
    p_adj_Lj$nsp1_p_adj < 0.05 , nsp1 , 0
  ),
  nsp2 = ifelse(
    p_adj_Lj$nsp2_p_adj < 0.05 , nsp2 , 0
  )
)]

res_Lj <- res_Lj[symrk != 0 | ccamk != 0 | nsp1 != 0 | nsp2 != 0]

Feat_tab_Lj <- metabolites_Lj
colnames(Feat_tab_Lj) <- paste(
  design_Lj_df$Genotype,
  rownames(design_Lj_df),
  sep = "_"
)
Feat_tab_Lj <- data.table(Feature = rownames(Feat_tab_Lj), Feat_tab_Lj)

res_Lj <- merge(res_Lj, Feat_tab_Lj, by = "Feature")
res_Lj <- merge(res_Lj, annotation_Lj, by = "Feature", all.x = TRUE)

# ## Removing duplicate features in annotation table
# n_anno_feature <- table(res_Lj$Feature)
# duplicated_features <- names(n_anno_feature)[n_anno_feature > 1]

# ## Identifying annotation with maximal probabilty for each duplicated feature
# idx_remove_lst <- list()
# for(i in 1:length(duplicated_features)){
#   idx_remove <- which(res_Lj$Feature == duplicated_features[i])
#   idx_keep <- which.max(res_Lj$`ClassyFire#class Probability`[idx_remove])
#   idx_remove <- idx_remove[-idx_keep]
#   idx_remove_lst[[i]] <- idx_remove
# }
# idx_remove <- unlist(idx_remove_lst)
# res_Lj <- res_Lj[-idx_remove]

fwrite(res_Lj, "../10_suppl_tables/Lotus_metabolite_test_results_tobit.csv")
fwrite(res_Lj, "1_tables/Lotus_metabolite_test_results_tobit.csv")

## Hordeum --------------------------------------------------------------------
# Setting factor levels
design_Hv[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

lm_res_Hv <- lm_over_features(metabolite_data_Hv, design_Hv)
p_adj_Hv <- apply(lm_res_Hv$p_val, 2, p.adjust, method = "fdr")
p_adj_pres_Hv <- apply(lm_res_Hv$p_val_pres, 2, p.adjust, method = "fdr")

# p_adj_Hv[highlighted_features_Hv,]
# p_adj_pres_Hv[highlighted_features_Hv,]

# Converting data tables to dataframes for the Tobit regression function
metabolites_Hv <- data.frame(
  metabolite_data_Hv[, -1],
  row.names = metabolite_data_Hv$Feature
)
design_Hv_df <- data.frame(design_Hv[, -1], row.names = design_Hv$Sample_ID)

# Fitting Tobit regression model with Genotype as covariate (WT as 
# reference level)
tobit_fit_Hv <- Tobit_model(
  metabolites = metabolites_Hv,
  meta_data = design_Hv_df,
  formula = ~Genotype
)

# Carrying out likelihood-ratio tests for non-reference level of genotype 
# variable
test_symrk_Hv <- Tobit_LRT(
  tobit_fit_Hv,
  var = "Genotypesymrk",
  p_adjust = "fdr"
)
test_ccamk_Hv <- Tobit_LRT(
  tobit_fit_Hv,
  var = "Genotypeccamk",
  p_adjust = "fdr"
)
test_nsp1_Hv <- Tobit_LRT(tobit_fit_Hv, var = "Genotypensp1", p_adjust = "fdr")
test_nsp2_Hv <- Tobit_LRT(tobit_fit_Hv, var = "Genotypensp2", p_adjust = "fdr")

# Setting up p-value matrix
p_adj_Hv <- data.table(
  Feature = rownames(test_symrk_Hv),
  symrk_p_adj = test_symrk_Hv$p_adj,
  ccamk_p_adj = test_ccamk_Hv$p_adj,
  nsp1_p_adj = test_nsp1_Hv$p_adj,
  nsp2_p_adj = test_nsp2_Hv$p_adj
)

# Setting up log-fold-change matrix
lfc_Hv <- data.table(Feature = rownames(tobit_fit_Hv$res), tobit_fit_Hv$res[,2:5])
colnames(lfc_Hv) <- gsub("Genotype", "Lfc_", colnames(lfc_Hv))

# Saving above results
fwrite(p_adj_Hv, "1_tables/p_adj_Hv.csv")
fwrite(lfc_Hv, "1_tables/lfc_Hv.csv")

# Saving results
res_Hv <- data.table(
  Feature = rownames(test_symrk_Hv),
  ifelse(tobit_fit_Hv$res[,2:5] > 0, 1, -1)
)
colnames(res_Hv) <- gsub("Genotype", "", colnames(res_Hv))

res_Hv[, ":="(
  symrk = ifelse(
    p_adj_Hv$symrk_p_adj < 0.05 , symrk , 0
  ),
  ccamk = ifelse(
    p_adj_Hv$ccamk_p_adj < 0.05 , ccamk , 0
  ),
  nsp1 = ifelse(
    p_adj_Hv$nsp1_p_adj < 0.05 , nsp1 , 0
  ),
  nsp2 = ifelse(
    p_adj_Hv$nsp2_p_adj < 0.05 , nsp2 , 0
  )
)]

res_Hv <- res_Hv[symrk != 0 | ccamk != 0 | nsp1 != 0 | nsp2 != 0]

Feat_tab_Hv <- metabolites_Hv
colnames(Feat_tab_Hv) <- paste(
  design_Hv_df$Genotype,
  rownames(design_Hv_df),
  sep = "_"
)
Feat_tab_Hv <- data.table(Feature = rownames(Feat_tab_Hv), Feat_tab_Hv)

res_Hv <- merge(res_Hv, Feat_tab_Hv, by = "Feature")
res_Hv <- merge(res_Hv, annotation_Hv, by = "Feature", all.x = TRUE)

fwrite(res_Hv, "../10_suppl_tables/Hordeum_metabolite_test_results_tobit.csv")
fwrite(res_Hv, "1_tables/Hordeum_metabolite_test_results_tobit.csv")
