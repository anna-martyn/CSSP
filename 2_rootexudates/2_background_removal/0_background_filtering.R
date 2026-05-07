# Seup ------------------------------------------------------------------------
library(data.table)

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading Tobit regression function
source("Tobit_function.R")

# Loading Lotus data ----------------------------------------------------------

metabolite_data_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_RootEx_Apr26_stdAUnew_featurelist.csv"
)

design_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_rootex_metadata.txt",
  drop = c(2, 4:6)
)

dir_name <- paste(
  "..",
  "1_data",
  "1_Lotus",
  "LotusCSSP_RootEx_Apr26_stdAUnew_canopus_structure_summary.tsv",
  sep = "/"
)
annotation_Lj <- fread(dir_name, drop = 21)

# Keeping only columns containing the peak area
cols <- colnames(metabolite_data_Lj)
cols_to_keep <- cols[grepl(":area", cols)]
metabolite_data_Lj <- metabolite_data_Lj[, c("id", cols_to_keep), with = FALSE]

# Replacing NAs with zeros
metabolite_data_Lj <- metabolite_data_Lj[,
  lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
  id
]

# Constructing samples numbers from column names
sample_numbers <- unlist(lapply(
  strsplit(colnames(metabolite_data_Lj)[-1], "_"),
  function(x) x[1]
))
sample_numbers <- gsub("datafile:", "", sample_numbers)
## Fixing typo in feature table
sample_numbers[26] <- "53"
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
  "../1_data/2_Hordeum/HordeumCSSP_RootEx_Apr26_stdAUnew_featurelist.csv"
)

design_Hv <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_rootex_metadata.txt",
  drop = c(2, 4:7)
)

dir_name <- paste(
  "..",
  "1_data",
  "2_Hordeum",
  "HordeumCSSP_RootEx_Apr26_stdAUnew_canopus_structure_summary.tsv",
  sep = "/"
)
annotation_Hv <- fread(dir_name, drop = 21)

# Keeping only columns containing the peak area
cols <- colnames(metabolite_data_Hv)
cols_to_keep <- cols[grepl(":area", cols)]
metabolite_data_Hv <- metabolite_data_Hv[, c("id", cols_to_keep), with = FALSE]

# Replacing NAs with zeros
metabolite_data_Hv <- metabolite_data_Hv[,
  lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
  id
]

# Constructing samples numbers from column names
sample_numbers <- unlist(lapply(
  strsplit(colnames(metabolite_data_Hv)[-1], "_"),
  function(x) x[4]
))
sample_numbers <- gsub("datafile:", "", sample_numbers)
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
present_features_idx <- rowMeans(metabolite_data_Lj[,-1] != 0) >= 0.1
present_features <- metabolite_data_Lj$Feature[present_features_idx]
metabolite_data_Lj <- metabolite_data_Lj[Feature %in% present_features]

# Keeping only Lotus features with annotations
metabolite_data_Lj <- metabolite_data_Lj[Feature %in% annotation_Lj$Feature]

# Keeping only Hordeum features present more than 10% of replicates
present_features_idx <- rowMeans(metabolite_data_Hv[,-1] != 0) >= 0.1
present_features <- metabolite_data_Hv$Feature[present_features_idx]
metabolite_data_Hv <- metabolite_data_Hv[Feature %in% present_features]

# Keeping only Hordeum features with annotations
metabolite_data_Hv <- metabolite_data_Hv[Feature %in% annotation_Hv$Feature]

# Filter based on sand control samples -----------------------------------------
# Check which features has statistically significantly higher intensity in at
# least one plant genotype compared to sand control samples, and keep only those
# features.

# Lotus
## Converting data tables to dataframes for the Tobit regression function
metabolites_Lj <- data.frame(
  metabolite_data_Lj[, -1],
  row.names = metabolite_data_Lj$Feature
)
design_Lj_df <- data.frame(design_Lj[, -1], row.names = design_Lj$Sample_ID)
design_Lj_df$Genotype <- factor(
  design_Lj_df$Genotype,
  levels = c("control", "WT", "symrk", "ccamk", "nsp1", "nsp2")
)

## Fitting Tobit regression model with Genotype as covariate (control as 
## reference level)
sand_filter_Lj <- Tobit_model(metabolites_Lj, design_Lj_df, formula = ~Genotype)

## Carrying out likelihood-ratio tests for non-reference level of genotype 
## variable
test_wt_Lj <- Tobit_LRT(sand_filter_Lj, var = "GenotypeWT", p_adjust = "fdr")
test_symrk_Lj <- Tobit_LRT(sand_filter_Lj, var = "Genotypesymrk", p_adjust = "fdr")
test_ccamk_Lj <- Tobit_LRT(sand_filter_Lj, var = "Genotypeccamk", p_adjust = "fdr")
test_nsp1_Lj <- Tobit_LRT(sand_filter_Lj, var = "Genotypensp1", p_adjust = "fdr")
test_nsp2_Lj <- Tobit_LRT(sand_filter_Lj, var = "Genotypensp2", p_adjust = "fdr")

## Filling dataframe with p-values from likelihood-ratio tests
p_val_mat_Lj <- data.frame(
  WT = test_wt_Lj$p_adj,
  symrk = test_symrk_Lj$p_adj,
  ccamk = test_ccamk_Lj$p_adj,
  nsp1 = test_nsp1_Lj$p_adj,
  nsp2 = test_nsp2_Lj$p_adj
)

## Filling dataframe with log-fold-changes from Tobit regression
lfc_mat_Lj <- data.frame(
  WT = test_wt_Lj$GenotypeWT,
  symrk = test_symrk_Lj$Genotypesymrk,
  ccamk = test_ccamk_Lj$Genotypeccamk,
  nsp1 = test_nsp1_Lj$Genotypensp1,
  nsp2 = test_nsp2_Lj$Genotypensp2
)

## Matrix indicating which features are present with statistically significantly 
## higher intensity than the control
feature_keep_mat_Lj <- p_val_mat_Lj < 0.05 & lfc_mat_Lj > 0
rownames(feature_keep_mat_Lj) <- rownames(test_wt_Lj)

feature_keep_Lj_idx <- apply(feature_keep_mat_Lj, 1, any)
feature_keep_Lj <- rownames(feature_keep_mat_Lj)[feature_keep_Lj_idx]

# Hordeum
## Converting data tables to dataframes for the Tobit regression function
metabolites_Hv <- data.frame(
  metabolite_data_Hv[, -1],
  row.names = metabolite_data_Hv$Feature
)
design_Hv_df <- data.frame(design_Hv[, -1], row.names = design_Hv$Sample_ID)
design_Hv_df$Genotype <- factor(
  design_Hv_df$Genotype,
  levels = c("control", "WT", "symrk", "ccamk", "nsp1", "nsp2")
)

## Fitting Tobit regression model with Genotype as covariate (control samples 
## as reference level)
sand_filter_Hv <- Tobit_model(metabolites_Hv, design_Hv_df, formula = ~Genotype)

## Carrying out likelihood-ratio tests for non-reference level of genotype 
## variable
test_wt_Hv <- Tobit_LRT(sand_filter_Hv, var = "GenotypeWT", p_adjust = "fdr")
test_symrk_Hv <- Tobit_LRT(sand_filter_Hv, var = "Genotypesymrk", p_adjust = "fdr")
test_ccamk_Hv <- Tobit_LRT(sand_filter_Hv, var = "Genotypeccamk", p_adjust = "fdr")
test_nsp1_Hv <- Tobit_LRT(sand_filter_Hv, var = "Genotypensp1", p_adjust = "fdr")
test_nsp2_Hv <- Tobit_LRT(sand_filter_Hv, var = "Genotypensp2", p_adjust = "fdr")

## Filling dataframe with p-values from likelihood-ratio tests
p_val_mat_Hv <- data.frame(
  WT = test_wt_Hv$p_adj,
  symrk = test_symrk_Hv$p_adj,
  ccamk = test_ccamk_Hv$p_adj,
  nsp1 = test_nsp1_Hv$p_adj,
  nsp2 = test_nsp2_Hv$p_adj
)

## Filling dataframe with log-fold-changes from Tobit regression
lfc_mat_Hv <- data.frame(
  WT = test_wt_Hv$GenotypeWT,
  symrk = test_symrk_Hv$Genotypesymrk,
  ccamk = test_ccamk_Hv$Genotypeccamk,
  nsp1 = test_nsp1_Hv$Genotypensp1,
  nsp2 = test_nsp2_Hv$Genotypensp2
)

## Matrix indicating which features are present with statistically significantly 
## higher intensity than the control
feature_keep_mat_Hv <- p_val_mat_Hv < 0.05 & lfc_mat_Hv > 0
rownames(feature_keep_mat_Hv) <- rownames(test_wt_Hv)

feature_keep_Hv_idx <- apply(feature_keep_mat_Hv, 1, any)
feature_keep_Hv <- rownames(feature_keep_mat_Hv)[feature_keep_Hv_idx]

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
