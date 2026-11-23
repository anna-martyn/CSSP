# Seup ------------------------------------------------------------------------
# Loading packages
pkg <- c("data.table", "codacore", "tensorflow", "caret")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
## Lotus
lotus_asv_table <- fread(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
colnames(lotus_asv_table)[1] <- "ASVid"

lotus_design <- fread(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", drop = c(5,7,8)
)
### Removing Nodule samples from metadata
lotus_design <- lotus_design[Compartment != "Nodules"]
soils <- c("NPK", "PK", "UF")
lotus_design[, Soil := factor(Soil, levels = soils)]

lotus_taxonomy <- fread("../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv")

## Hordeum
hordeum_asv_table <- fread(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv"
)
colnames(hordeum_asv_table)[1] <- "ASVid"

hordeum_design <- fread(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", drop = c(2,3,8)
)
### Removing soil samples from metadata
hordeum_design <- hordeum_design[Compartment != "Soil"]
hordeum_design[, Soil := factor(Soil, levels = soils)]

hordeum_taxonomy <- fread(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv"
)

# Cleaning up taxonomy
taxa_clean <- function(taxonomy){
  taxa_levels <- c(
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  )
  taxonomy[,c(taxa_levels):= tstrsplit(Taxon, "; ", fill = "u__Unknown")]
  taxonomy[,c(taxa_levels):= lapply(.SD, substr, 4, 1000), .SDcols = taxa_levels]
  taxonomy[Kingdom == "ssigned", Kingdom:= "Unassigned"]
  taxonomy[,Taxon:=NULL]
  setcolorder(taxonomy, c("Feature ID", taxa_levels, "Confidence"))
  return(taxonomy)
}

lotus_taxonomy <- taxa_clean(lotus_taxonomy)
hordeum_taxonomy <- taxa_clean(hordeum_taxonomy)

# Removing non-bacterial ASVs from 
non_bac_asv <- lotus_taxonomy[Kingdom != "Bacteria", "Feature ID"][[1]]
lotus_asv_table <- lotus_asv_table[!(ASVid %in% non_bac_asv)]

# Setting up list for easy retrieval
asv_table <- list(Lotus = lotus_asv_table, Hordeum = hordeum_asv_table)
design <- rbind(lotus_design, hordeum_design)
taxonomy <- rbind(lotus_taxonomy, hordeum_taxonomy)
setnames(taxonomy, "Feature ID", "ASVid")

# Setting up RA tables
convert_to_RA <- function(asv_table){
  data.table(
    ASVid = asv_table$ASVid,
    t(t(asv_table[,-1])/colSums(asv_table[,-1]))
  )
}
asv_table_RA <- lapply(asv_table, convert_to_RA)

# Functions for nested classification model -----------------------------------
Nested_model <- function(asv_table, y, seed, lambda, overlap = T){
  # Setting up NPK vs non-NPK binary variable
  NPK_bin_var <- ifelse(y == "NPK", "NPK", "Non-NPK")
  NPK_bin_var <- relevel(factor(NPK_bin_var), ref = "NPK")
  
  # Tranposing ASV table
  asv_table_transp <- transpose(
    l = asv_table,
    make.names = "ASVid",
    keep.names = "SampleID"
  )

  # Training codacore model for binary classification between NPK and non-NPK soils
  set.seed(seed)
  tf$random$set_seed(seed)
  Model1 <- codacore(
    x = data.frame(asv_table_transp[,-1], row.names = asv_table_transp$SampleID),
    y = NPK_bin_var,
    objective = "binary classification",
    logRatioType = "SLR",
    lambda = lambda,
    overlap = overlap
  )
  
  # Removing NPK samples
  Non_NPK_idx <- which(y != "NPK")
  y_non_NPK <- y[Non_NPK_idx]
  asv_table_non_NPK_transp <- asv_table_transp[Non_NPK_idx]
  
  # Training codacore model for binary classification between PK and UF soils
  Model2 <- codacore(
    x = data.frame(
      asv_table_non_NPK_transp[,-1],
      row.names = asv_table_non_NPK_transp$SampleID
    ),
    y = y_non_NPK,
    objective = "binary classification",
    logRatioType = "SLR",
    lambda = lambda,
    overlap = overlap
  )
  
  return( list(Model1 = Model1, Model2 = Model2) )
}

# Function that returns predtictions based on a nested model fit obtained from 
# the Nested_model function
Nested_model_pred <- function(Model, test_data, thresh = 0.5){
  # Transformation that converts predicted log-odds into probabilities
  inverse_logit <- function(x) exp(x)/(1+exp(x))
  
  # Transposing test data
  test_data_transp <- transpose(
    test_data, keep.names = "SampleID", make.names = "ASVid"
  )

  # Obtain predictions between NPK and non-NPK soils from Model 1
  pred <- predict(
    Model$Model1,
    data.frame(
      test_data_transp[,-1],
      row.names = test_data_transp$SampleID
    ),
    logits = F
  )
  # Converting log-odds into probabilities and obtaining predicted soil-type
  prob <- inverse_logit(pred)
  prob[is.nan(prob)] <- 1
  pred <- ifelse(prob < thresh, "NPK", "non-NPK")
  
  # Keeping only samples from test data predicted to be from non-NPK soils
  lower_idx <- which(pred == "non-NPK")
  test_data_lower_transp <- test_data_transp[lower_idx]
  # Obtaining prediction from Model 2, PK vs UF
  pred_lower <- predict(
    Model$Model2,
    data.frame(
      test_data_lower_transp[,-1],
      row.names = test_data_lower_transp$SampleID
    ),
    logits = F
  )
  # Converting log-odds into probabilities and obtaining predicted soil-type
  prob <- inverse_logit(pred_lower)
  prob[is.nan(prob)] <- 1
  pred_lower <- ifelse(prob < thresh, "PK", "UF")
  
  # Returning full vector of trinary predictions
  pred[lower_idx] <- pred_lower
  return(pred)
}

# Run the prediction analysis --------------------------------------------------
Opt <- expand.grid(
  Host = c("Lotus", "Hordeum"),
  Compartment = c("Root", "Rhizosphere")
)

design_pred_list <- list()
for(i in 1:nrow(Opt)){
  run <- paste(as.matrix(Opt)[i,], collapse = "_")
  # cat(run, "\r")
  
  # Setting up the current host-compartment combination for the current iteration
  current_host <- Opt$Host[i]
  current_comp <- Opt$Compartment[i]

  # Subsetting ASV table and metadata to only contain samples from current 
  # plant-compartment combination
  asv_table_sub <- asv_table[[current_host]]
  design_sub <- design[
    Plant == current_host & Compartment == current_comp
  ]
  sample_sub <- design_sub$SampleID
  asv_table_sub <- asv_table_sub[,c("ASVid", sample_sub), with = F]
  
  # Filtering out ASVs present in less than 10% of samples
  pres_idx <- which(rowMeans(asv_table_sub[,-1] != 0) >= 0.1)
  asv_table_sub <- asv_table_sub[pres_idx]
  
  # Leave-one-genotype-out cross validation
  gt <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  pred_list <- list()
  for (j in 1:5) {
    cat(paste0(run, "_", j, "      "), "\r")
    
    test_samples <- design_sub$SampleID[design_sub$Genotype == gt[j]]
    training_samples <- setdiff(design_sub$SampleID, test_samples)

    design_train <- design_sub[SampleID %in% training_samples]
    design_test <- design_sub[SampleID %in% test_samples]

    x_train <- asv_table_sub[, c("ASVid", design_train$SampleID), with = F]
    x_test <- asv_table_sub[, c("ASVid", design_test$SampleID), with = F]

    y_train <- design_train$Soil
    y_test <- design_test$Soil

    tf$random$set_seed(1700294030)
    cc <- Nested_model(
      asv_table = x_train,
      y = as.character(y_train),
      seed = 1700299435,
      lambda = 1
    )
    
    pred <- Nested_model_pred(cc, x_test) |> factor(levels = soils)
    design_test[, Soil_pred := pred]
    pred_list[[j]] <- design_test
  }

  design_pred_list <- c(design_pred_list, pred_list)
}

design_pred <- rbindlist(design_pred_list)

pred_res <- design_pred[,
  .(
    Accuracy = mean(Soil == Soil_pred),
    Accuracy_CI_lower = binom.test(sum(Soil == Soil_pred), .N)$conf.int[1],
    Accuracy_CI_upper = binom.test(sum(Soil == Soil_pred), .N)$conf.int[2]
  ),
  .(Plant, Compartment)
]
pred_res <- pred_res[c(3,1,4,2)]

cm <- list()
for (i in 1:4) {
  pred_sub <- design_pred[Compartment == Opt$Compartment[i] & Plant == Opt$Host[i]]
  cm[[i]] <- confusionMatrix(pred_sub$Soil_pred, pred_sub$Soil)
}
names(cm) <- paste("Conf. mat", Opt[,1], Opt[,2])

cm_dt <- lapply(cm, function(x){
  mat <- x$table |> as.numeric() |> matrix(ncol = 3) |> data.frame()
  dt <- data.table(Predicted = soils, mat)
  colnames(dt)[-1] <- paste("Obs", soils, sep = "_")
  return(dt)
})

pred_res_soil <- design_pred[,
  .(
    Accuracy = mean(Soil == Soil_pred),
    Accuracy_CI_lower = binom.test(sum(Soil == Soil_pred), .N)$conf.int[1],
    Accuracy_CI_upper = binom.test(sum(Soil == Soil_pred), .N)$conf.int[2]
  ),
  .(Plant, Compartment, Soil)
]

setnames(design_pred, "Soil_pred", "Predicted soil type")

res_to_save <- list(
  "Full results" = design_pred,
  "Overall accuracy" = pred_res,
  "Soil specific accuracy" = pred_res_soil
)
res_to_save <- c(res_to_save, cm_dt)

# Saving results
write.xlsx(
  res_to_save,
  file = "../3_suppl_tables/Suppl_Table_X_pred_leave_one_genotype_out.xlsx"
)

