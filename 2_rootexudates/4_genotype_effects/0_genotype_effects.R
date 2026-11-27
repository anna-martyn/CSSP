# Seup ------------------------------------------------------------------------
library(data.table)
library(VGAM)
library(ggplot2)

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading Tobit regression function
source("Tobit_function.R")
genotypes <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")

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
  drop = 4:7
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

# Genotype effects ------------------------------------------------------------
## Lotus ----------------------------------------------------------------------
# Setting factor levels
design_Lj[,
  Genotype := factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
]

# Converting data tables to dataframes for the Tobit regression function
metabolites_Lj <- data.frame(
  metabolite_data_Lj[, -1],
  row.names = metabolite_data_Lj$Feature
)
design_Lj_df <- data.frame(design_Lj[, -1], row.names = design_Lj$Sample_ID)

feat_keep_idx <- rowSums(metabolites_Lj != 0) >= 3
feat_keep <- rownames(metabolites_Lj)[feat_keep_idx]
metabolites_Lj <- metabolites_Lj[feat_keep, ]

# Fitting Tobit regression model with Genotype as covariate (WT as 
# reference level)
tobit_fit_Lj <- Tobit_model2(
  metabolites = metabolites_Lj,
  meta_data = design_Lj_df,
  formula = ~Genotype,
  ref = "WT",
  test = "two-sided"
)

# Identifying non-significant features still of interest, i.e. features 
# completely absent from a mutant, but present in the WT or vice-versa
no_test <- matrix(NA, nrow = nrow(metabolites_Lj), ncol = 5)
rownames(no_test) <- rownames(metabolites_Lj)
colnames(no_test) <- genotypes
for (i in 1:5) {
  idx <- design_Lj[Genotype == genotypes[i], Sample_ID]
  no_test[,i] <- apply(metabolites_Lj[,idx] == 0, 1, all)
}
no_test_dt <- data.table(
  Feature = rownames(no_test),
  no_test[, 1] - no_test[, -1]
)
no_test_dt_Lj <- no_test_dt[symrk != 0 | ccamk != 0 | nsp1 != 0| nsp2 != 0]
# fwrite(no_test_dt, file = "1_tables/non_sig_features_of_interest_Lotus.csv")

# p-values
p_vals_Lj <- tobit_fit_Lj$p_vals
p_vals_Lj[is.na(p_vals_Lj)] <- tobit_fit_Lj$p_vals_pa[is.na(p_vals_Lj)]

p_adj_Lj <- matrix(p.adjust(p_vals_Lj[,-1], "fdr"), ncol = 4)
rownames(p_adj_Lj) <- rownames(p_vals_Lj)
colnames(p_adj_Lj) <- colnames(p_vals_Lj)[-1]

# Residual plot
set.seed(1793694715)
all_res <- unlist(tobit_fit_Lj$resid)
dt <- data.table(
  resid = c(all_res, rnorm(length(all_res))),
  Density = rep(c("Tobit residuals", "Simulated\nstd. normal"), each = length(all_res))
)

resid_plot_Lj <- ggplot(dt, aes(x = resid, fill = Density, colour = Density)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  theme_bw()+
  scale_fill_manual(values = c("#1a80bb", "#a00000"))+
  scale_color_manual(values = c("#1a80bb", "#a00000"))+
  labs(x = "Generalised residuals", y = "Probability density")+
  NULL

# Setting up p-value matrix
p_adj_Lj <- data.table(Feature = rownames(p_adj_Lj),  p_adj_Lj)

# Setting up log-fold-change matrix
lfc_Lj_temp <- tobit_fit_Lj$res[,2:5]
lfc_Lj_temp[is.na(lfc_Lj_temp)] <- tobit_fit_Lj$beta_est_pa[,-1][is.na(lfc_Lj_temp)]
lfc_Lj <- data.table(Feature = rownames(lfc_Lj_temp), lfc_Lj_temp)
colnames(lfc_Lj) <- gsub("Genotype", "Lfc_", colnames(lfc_Lj))

# Saving above results
fwrite(p_adj_Lj, "1_tables/p_adj_Lj.csv")
fwrite(lfc_Lj, "1_tables/lfc_Lj.csv")

# Saving results
res_Lj <- data.table(
  Feature = rownames(lfc_Lj_temp),
  ifelse(lfc_Lj_temp > 0, 1, -1)
)
colnames(res_Lj) <- gsub("Genotype", "", colnames(res_Lj))

res_Lj[, ":="(
  symrk = ifelse(
    p_adj_Lj$Genotypesymrk < 0.05 , symrk , 0
  ),
  ccamk = ifelse(
    p_adj_Lj$Genotypeccamk < 0.05 , ccamk , 0
  ),
  nsp1 = ifelse(
    p_adj_Lj$Genotypensp1 < 0.05 , nsp1 , 0
  ),
  nsp2 = ifelse(
    p_adj_Lj$Genotypensp2 < 0.05 , nsp2 , 0
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

# Converting data tables to dataframes for the Tobit regression function
metabolites_Hv <- data.frame(
  metabolite_data_Hv[, -1],
  row.names = metabolite_data_Hv$Feature
)
design_Hv_df <- data.frame(design_Hv[, -1], row.names = design_Hv$Sample_ID)

feat_keep_idx <- rowSums(metabolites_Hv != 0) >= 3
feat_keep <- rownames(metabolites_Hv)[feat_keep_idx]
metabolites_Hv <- metabolites_Hv[feat_keep, ]

# Fitting Tobit regression model with Genotype as covariate (WT as 
# reference level)
tobit_fit_Hv <- Tobit_model2(
  metabolites = metabolites_Hv,
  meta_data = design_Hv_df,
  formula = ~Genotype,
  ref = "WT",
  test = "two-sided",
  presence_absence = TRUE
)

# Identifying non-significant features still of interest, i.e. features 
# completely absent from a mutant, but present in the WT or vice-versa
no_test <- matrix(NA, nrow = nrow(metabolites_Hv), ncol = 5)
rownames(no_test) <- rownames(metabolites_Hv)
colnames(no_test) <- genotypes
for (i in 1:5) {
  idx <- design_Hv[Genotype == genotypes[i], Sample_ID]
  no_test[,i] <- apply(metabolites_Hv[,idx] == 0, 1, all)
}
no_test_dt <- data.table(
  Feature = rownames(no_test),
  no_test[, 1] - no_test[, -1]
)
no_test_dt_Hv <- no_test_dt[symrk != 0 | ccamk != 0 | nsp1 != 0| nsp2 != 0]
# fwrite(no_test_dt, file = "1_tables/non_sig_features_of_interest_Hordeum.csv")

# p-values
p_vals_Hv <- tobit_fit_Hv$p_vals
p_vals_Hv[is.na(p_vals_Hv)] <- tobit_fit_Hv$p_vals_pa[is.na(p_vals_Hv)]
# p_adj_Hv <- apply(p_vals_Hv[,-1], 2, p.adjust, method = "fdr")

p_adj_Hv <- matrix(p.adjust(p_vals_Hv[,-1], "fdr"), ncol = 4)
rownames(p_adj_Hv) <- rownames(p_vals_Hv)
colnames(p_adj_Hv) <- colnames(p_vals_Hv)[-1]

# Residual plot
set.seed(1793694715)
all_res <- unlist(tobit_fit_Hv$resid)
dt <- data.table(
  resid = c(all_res, rnorm(length(all_res))),
  Density = rep(c("Tobit residuals", "Simulated\nstd. normal"), each = length(all_res))
)

resid_plot_Hv <- ggplot(dt, aes(x = resid, fill = Density, colour = Density)) +
  geom_density(alpha = 0.4, linewidth = 0.5) +
  theme_bw()+
  scale_fill_manual(values = c("#1a80bb", "#a00000"))+
  scale_color_manual(values = c("#1a80bb", "#a00000"))+
  theme()+
  labs(x = "Generalised residuals", y = "Probability density")+
  NULL

# Setting up p-value matrix
p_adj_Hv <- data.table(Feature = rownames(p_adj_Hv), p_adj_Hv)

# Setting up log-fold-change matrix
lfc_Hv <- data.table(Feature = rownames(tobit_fit_Hv$res), tobit_fit_Hv$res[,2:5])
colnames(lfc_Hv) <- gsub("Genotype", "Lfc_", colnames(lfc_Hv))

# Setting up log-fold-change matrix
lfc_Hv_temp <- tobit_fit_Hv$res[,2:5]
lfc_Hv_temp[is.na(lfc_Hv_temp)] <- tobit_fit_Hv$beta_est_pa[,-1][is.na(lfc_Hv_temp)]
lfc_Hv <- data.table(Feature = rownames(lfc_Hv_temp), lfc_Hv_temp)
colnames(lfc_Hv) <- gsub("Genotype", "Lfc_", colnames(lfc_Hv))

# Saving above results
fwrite(p_adj_Hv, "1_tables/p_adj_Hv.csv")
fwrite(lfc_Hv, "1_tables/lfc_Hv.csv")

# Saving results
res_Hv <- data.table(
  Feature = rownames(lfc_Hv_temp),
  ifelse(lfc_Hv_temp > 0, 1, -1)
)
colnames(res_Hv) <- gsub("Genotype", "", colnames(res_Hv))

res_Hv[, ":="(
  symrk = ifelse(
    p_adj_Hv$Genotypesymrk < 0.05 , symrk , 0
  ),
  ccamk = ifelse(
    p_adj_Hv$Genotypeccamk < 0.05 , ccamk , 0
  ),
  nsp1 = ifelse(
    p_adj_Hv$Genotypensp1 < 0.05 , nsp1 , 0
  ),
  nsp2 = ifelse(
    p_adj_Hv$Genotypensp2 < 0.05 , nsp2 , 0
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

# Residual suppl. figure
resid_plot <- ggpubr::ggarrange(
  resid_plot_Lj,
  resid_plot_Hv,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  "../9_final_figures/Suppl_Figure_X_tobit_resid.pdf",
  resid_plot,
  width = 160,
  height = 90,
  units = "mm"
)

# Zero vs. non-zero tables
openxlsx::write.xlsx(
  list("Lotus" = no_test_dt_Lj, "Hordeum" = no_test_dt_Hv),
  file = "../10_suppl_tables/suppl_tab_metabolites_zero_vs_nonzero.xlsx"
)
