# Seup ------------------------------------------------------------------------
library(data.table)
library(VGAM)
library(logistf)

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading Tobit regression function
source("Tobit_function.R")

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
metabolite_data_Hv <- metabolite_data_Hv[,
  lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
  id
]

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

all(colnames(metabolites_Lj) == rownames(design_Lj_df))
p_vals_pa <- matrix(NA, ncol = 5, nrow = nrow(metabolites_Lj))
for (i in 1:nrow(metabolites_Lj)) {
  cat(i, "\r")
  df <- data.frame(
    design_Lj_df,
    intensity = ifelse(as.numeric(metabolites_Lj[i, ]) == 0, 0, 1)
  )
  lg <- logistf(intensity~Genotype, data = df)
  p_vals_pa[i,] <- lg$prob[-1]
}

p_adj_pa <- apply(p_vals_pa, 2, p.adjust, method = "fdr")
keep_idx <- apply(p_adj_pa < 0.05, 1, any)
feature_keep_pa <- rownames(metabolites_Lj)[keep_idx]

feature_other <- setdiff(rownames(metabolites_Lj), feature_keep_pa)

# cidx <- design_Lj[Genotype == "control", Sample_ID]
# cidx <- which(colnames(metabolites_Lj) %in% cidx)
# ct_all0 <- apply(metabolites_Lj[,cidx], 1, function(x) all(x == 0))
# zcidx <- which(ct_all0)

# pres_met <- apply(metabolites_Lj[,-cidx], 1, function(x) mean(x != 0) > 0.1)
# feature_keep_pa <- names(ct_all0)[ct_all0 & pres_met]

set.seed(1793694715)
tobit_fit_Lj <- Tobit_model2(
  metabolites = metabolites_Lj[feature_other,],
  meta_data = design_Lj_df,
  formula = ~Genotype,
  ref = "control",
  test = "one-sided",
  presence_absence = FALSE
)

p_val <- tobit_fit_Lj$p_vals
p_val[is.na(p_val)] <- 1
p_adj <- apply(p_val[,-1], 2, p.adjust, method = "fdr")

feature_keep_tr <- apply(p_adj < 0.05, 1, any, na.rm = T)
feature_keep_tr <- rownames(p_adj)[feature_keep_tr]

feature_keep_Lj <- union(feature_keep_pa, feature_keep_tr)

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

all(colnames(metabolites_Hv) == rownames(design_Hv_df))
p_vals_pa <- matrix(NA, ncol = 5, nrow = nrow(metabolites_Hv))
for (i in 1:nrow(metabolites_Hv)) {
  cat(i, "\r")
  df <- data.frame(
    design_Hv_df,
    intensity = ifelse(as.numeric(metabolites_Hv[i, ]) == 0, 0, 1)
  )
  lg <- logistf(intensity~Genotype, data = df)
  p_vals_pa[i,] <- lg$prob[-1]
}

p_adj_pa <- apply(p_vals_pa, 2, p.adjust, method = "fdr")
keep_idx <- apply(p_adj_pa < 0.05, 1, any)
feature_keep_pa <- rownames(metabolites_Hv)[keep_idx]
feature_other <- setdiff(rownames(metabolites_Hv), feature_keep_pa)

# cidx <- design_Hv[Genotype == "control", Sample_ID]
# cidx <- which(colnames(metabolites_Hv) %in% cidx)
# ct_all0 <- apply(metabolites_Hv[,cidx], 1, function(x) all(x == 0))
# zcidx <- which(ct_all0)

# pres_met <- apply(metabolites_Hv[,-cidx], 1, function(x) mean(x != 0) > 0.1)
# feature_keep_pa <- names(ct_all0)[ct_all0 & pres_met]

set.seed(1793694715)
tobit_fit_Hv <- Tobit_model2(
  metabolites = metabolites_Hv[feature_other,],
  meta_data = design_Hv_df,
  formula = ~Genotype,
  ref = "control",
  test = "one-sided",
  presence_absence = FALSE
)

p_val <- tobit_fit_Hv$p_vals
p_val[is.na(p_val)] <- 1
p_adj <- apply(p_val[,-1], 2, p.adjust, method = "fdr")

feature_keep_tr <- apply(p_adj[,-1] < 0.05, 1, any, na.rm = T)
feature_keep_tr <- rownames(p_adj)[feature_keep_tr]
feature_keep_Hv <- union(feature_keep_pa, feature_keep_tr)

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

# Removing metabolites with less than three non-zero observations
non_zero_obs <- rowSums(metabolite_data_Lj != 0) >= 3
metabolite_data_Lj <- metabolite_data_Lj[non_zero_obs,]

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

# Removing metabolites with less than three non-zero observations
non_zero_obs <- rowSums(metabolite_data_Hv != 0) >= 3
metabolite_data_Hv <- metabolite_data_Hv[non_zero_obs,]

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
