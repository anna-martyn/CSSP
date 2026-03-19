# Seup ------------------------------------------------------------------------
# Loading packages
pkg <- c(
  "data.table", "magrittr", "ggplot2", "ggh4x", "codacore", 
  "tensorflow", "ggpubr", "gridExtra", "cowplot", "ggtext"
)
for(pk in pkg){
  library(pk, character.only = T)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
## Lotus
lotus_asv_table <- fread(
  "../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
colnames(lotus_asv_table)[1] <- "ASVid"

lotus_design <- fread(
  "../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", drop = c(5,7,8)
)
### Removing Nodule samples from metadata
lotus_design <- lotus_design[Compartment != "Nodules"]

lotus_taxonomy <- fread("../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv")

## Hordeum
hordeum_asv_table <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv"
)
colnames(hordeum_asv_table)[1] <- "ASVid"

hordeum_design <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", drop = c(2,3,8)
)
### Removing soil samples from metadata
hordeum_design <- hordeum_design[Compartment != "Soil"]

hordeum_taxonomy <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv"
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

Opt <- expand.grid(
  Host = c("Lotus", "Hordeum"),
  Compartment = c("Root", "Rhizosphere")
)

# Functions for nested classification model -----------------------------------

# The following function distinguishes between 'highest order' and 'lower order' 
# levels. 'Highest order' refers to the level that the first predicts for, while 
# the lower levels refers to the levels the subsequent model predicts between. In 
# practice, we want the first model to distinguish between NPK and non-NPK soils, 
# so NPK is the 'highest order level'. The second model distinguishes between PK 
# and UF soils, so those are referred to as 'lower order levels'.
Nested_model <- function(asv_table, y, highest_order, seed, lambda, overlap = T){
  # Setting up binary variables that checks if a samples belongs to the highest 
  # order level or not.
  Highest_order_resp <- ifelse(
    y == highest_order, highest_order, paste("Not", highest_order)
  )
  Highest_order_resp <- relevel(factor(Highest_order_resp), ref = highest_order)
  
  # Tranposing ASV table
  asv_table_transp <- transpose(
    l = asv_table,
    make.names = "ASVid",
    keep.names = "SampleID"
  )

  # Training codacore model for binary classification between highest order 
  # level or not (in Practice: NPK vs. non-NPK in practice)
  set.seed(seed)
  tf$random$set_seed(seed)
  cc1 <- codacore(
    x = data.frame(asv_table_transp[,-1], row.names = asv_table_transp$SampleID),
    y = Highest_order_resp,
    objective = "binary classification",
    logRatioType = "SLR",
    lambda = lambda,
    overlap = overlap
  )
  
  # Removing samples that do not belong to the highest order level 
  # (in practice: remove NPK samples)
  Not_highest_idx <- which(y != highest_order)
  y2 <- y[Not_highest_idx]
  asv_table_lower_transp <- asv_table_transp[Not_highest_idx]
  
  # Training codacore model for binary classification between the lower order levels 
  # (In practice: Predict if a samples is from PK or UF soil)
  cc2 <- codacore(
    x = data.frame(
      asv_table_lower_transp[,-1],
      row.names = asv_table_lower_transp$SampleID
    ),
    y = y2,
    objective = "binary classification",
    logRatioType = "SLR",
    lambda = lambda,
    overlap = overlap
  )
  
  return( list(Highest = cc1, Lowest = cc2) )
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

  # Obtain predictions from highest order model 
  # (In practice: NPK vs. non-NPK prediction)
  pred <- predict(
    Model$Highest,
    data.frame(test_data_transp[,-1], row.names = test_data_transp$SampleID),
    logits = F
  )
  # Converting log-odds into probabilities and obtaining predicted soil-type
  prob <- inverse_logit(pred)
  prob[is.nan(prob)] <- 1
  pred <- ifelse(prob < thresh, "NPK", "Not NPK")
  
  # Keeping only samples from test data predicted to be from non-NPK soil
  lower_idx <- which(pred == "Not NPK")
  test_data_lower_transp <- test_data_transp[lower_idx]
  # Obtaining prediction from Model 2, PK vs UF
  pred_lower <- predict(
    Model$Lowest,
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
Ratio_summary_list <- list()
Res_list <- list()
Ratio_amount_list <- list()
tax_summary_list <- list()
for(i in 1:nrow(Opt)){
  run <- paste(as.matrix(Opt)[i,], collapse = "_")
  Current_plant <- as.matrix(Opt)[i,1]
  cat(run, "\r")

  asv_table_sub <- asv_table[[Opt$Host[i]]]
  design_sub <- design[
    Plant == Opt$Host[i] & Compartment == Opt$Compartment[i]
  ]
  sample_sub <- design_sub$SampleID
  asv_table_sub <- asv_table_sub[,c("ASVid", sample_sub), with = F]
  
  # Filter out ASVs present in less than 10% of samples
  pres_idx <- which(rowMeans(asv_table_sub[,-1] != 0) >= 0.1)
  asv_table_sub <- asv_table_sub[pres_idx]

  # Splitting into test and training data
  set.seed(1700294030)
  test_samples <- design_sub[
    ,.(SampleID = sample(SampleID, 2, replace = FALSE)), by = .(Genotype, Soil)
  ]$SampleID
  training_samples <- setdiff(design_sub$SampleID, test_samples)

  x_train <- asv_table_sub[,c("ASVid", training_samples), with = F]
  x_test <- asv_table_sub[,c("ASVid", test_samples), with = F]

  y_train <- design_sub[SampleID %in% training_samples, Soil]
  y_test <- design_sub[SampleID %in% test_samples, Soil]

  # Fitting prediction model
  tf$random$set_seed(1700294030)
  cc <- Nested_model(
    OTU_table = x_train,
    y = y_train,
    highest_order = "NPK",
    seed = 1700299435,
    lambda = 1
  )

  pred <- Nested_model_pred(cc, x_test)
  eval <- mean(y_test == pred)
  
  # Summarising results
  dt <- data.table(SampleID = names(pred), Obs = y_test, Pred = pred)
  dt <- dt[order(Obs)]
  dt[,":="(Host = Opt$Host[i], Compartment = Opt$Compartment[i])]
  dt <- merge(dt, design_sub[,c(1,3)], by = "SampleID")
  Res_list[[i]] <- dt

  # Extracting ASVs used in predictions
  num_vec_NPK <- lapply(
    cc$Highest$ensemble, function(x) x_train$ASVid[x$hard$numerator]
  )
  denom_vec_NPK <- lapply(
    cc$Highest$ensemble, function(x) x_train$ASVid[x$hard$denominator]
  )

  num_vec_PK <- lapply(
    cc$Lowest$ensemble, function(x) x_train$ASVid[x$hard$numerator]
  )
  denom_vec_PK <- lapply(
    cc$Lowest$ensemble, function(x) x_train$ASVid[x$hard$denominator]
  )

  # Coutning the ASVs in the ratios for table C
  Mat_temp <- matrix(
    c(
      unlist(lapply(num_vec_NPK, length)),
      unlist(lapply(denom_vec_NPK, length))
    ),
    ncol = 2
  )
  # NPK_ASVs_in_Ratios <- apply(Mat_temp, 1, function(x) paste(x[1], x[2], sep = "/"))
  NPK_ASVs_in_Ratios <- rowSums(Mat_temp)
  NPK_ratios <- paste0(
    length(num_vec_NPK), " (",
    paste(NPK_ASVs_in_Ratios, collapse = ", "), ")"
  )

  Mat_temp <- matrix(
    c(
      unlist(lapply(num_vec_PK, length)),
      unlist(lapply(denom_vec_PK, length))
    ),
    ncol = 2
  )
  # PK_ASVs_in_Ratios <- apply(Mat_temp, 1, function(x) paste(x[1], x[2], sep = "/"))
  PK_ASVs_in_Ratios <- rowSums(Mat_temp)
  PK_ratios <- paste0(
    length(num_vec_PK), " (",
    paste(PK_ASVs_in_Ratios, collapse = ", "), ")"
  )

  Ratio_amount_list[[i]] <- data.table(
    Host = Opt$Host[i],
    Compartment = Opt$Compartment[i],
    "NPK ratios" = NPK_ratios,
    "PK ratios" = PK_ratios,
    Accuracy = eval
  )

  # Constructing supplementary table
  Sum_num_NPK <- lapply(
    1:length(num_vec_NPK), function(i) data.table(
      ASVid = num_vec_NPK[[i]],
      Role = paste0("Numerator", i),
      Prediction = "NPK vs. non-NPK"
    )
  )
  Sum_denom_NPK <- lapply(
    1:length(denom_vec_NPK), function(i) data.table(
      ASVid = denom_vec_NPK[[i]],
      Role = paste0("Denominator", i),
      Prediction = "NPK vs. non-NPK"
    )
  )

  Sum_num_PK <- lapply(
    1:length(num_vec_PK), function(i) data.table(
      ASVid = num_vec_PK[[i]],
      Role = paste0("Numerator", i),
      Prediction = "PK vs. UF"
    )
  )
  Sum_denom_PK <- lapply(
    1:length(denom_vec_PK), function(i) data.table(
      ASVid = denom_vec_PK[[i]],
      Role = paste0("Denominator", i),
      Prediction = "PK vs. UF"
    )
  )

  lst <- c(Sum_num_NPK, Sum_denom_NPK, Sum_num_PK, Sum_denom_PK)
  Ratio_summary_temp <- rbindlist(lst)
  Ratio_summary_temp[,":="(Compartment = Opt$Compartment[i], Host = Opt$Host[i])]
  pred_ASVs <- unique(Ratio_summary_temp$ASVid)
  
  asv_table_RA_sub <- asv_table_RA[[Opt$Host[i]]]
  asv_table_RA_sub <- asv_table_RA_sub[ASVid %in% pred_ASVs]
  asv_table_RA_sub <- transpose(
    asv_table_RA_sub, keep.names = "SampleID", make.names = "ASVid"
  )
  asv_table_RA_sub <- merge(design_sub, asv_table_RA_sub, by = "SampleID")
  mean_RAs <- asv_table_RA_sub[,lapply(.SD, mean), Soil, .SDcols = pred_ASVs]
  mean_RAs <- transpose(mean_RAs, make.names = "Soil", keep.names = "ASVid")
  colnames(mean_RAs)[-1] <- paste(colnames(mean_RAs)[-1], "RA", sep = "_")

  Ratio_summary_temp <- merge(Ratio_summary_temp, mean_RAs)
  Ratio_summary_list[[i]] <- Ratio_summary_temp

  # Summarising RA and taxonomic information for ASVs used for estimation
  ASV_pred <- unique(Ratio_summary_temp$ASVid)
  asv_table_RA_sub <- asv_table_RA[[Opt$Host[i]]]
  asv_table_RA_sub <- asv_table_RA_sub[ASVid %in% ASV_pred]
  asv_table_RA_sub <- asv_table_RA_sub[,c("ASVid", sample_sub), with = F]
  asv_table_RA_sub <- merge(taxonomy[,c(1, 5)], asv_table_RA_sub, by = "ASVid")
  Order_table_RA_sub <- asv_table_RA_sub[,lapply(.SD, mean), Order, .SDcols = sample_sub]
  Order_table_RA_sub_T <- transpose(
    Order_table_RA_sub, make.names = "Order", keep.names = "SampleID"
  )
  Order_table_RA_sub_T <- merge(design_sub, Order_table_RA_sub_T, by = "SampleID")
  Order_info <- melt(
    Order_table_RA_sub_T,
    id.vars = 1:5,
    variable.name = "Order",
    value.name = "RA"
  )
  tax_summary_list[[i]] <- Order_info[,.(RA = mean(RA)), .(Plant, Compartment, Soil, Order)]
}

tax_summary <- rbindlist(tax_summary_list)
tax_summary[,.(RA = sum(RA)), .(Plant, Compartment, Soil)]
fwrite(tax_summary, "Pred_taxonomic_composition.csv")

Res <- rbindlist(Res_list)
Ratio_summary <- rbindlist(Ratio_summary_list)
Ratio_summary <- merge(Ratio_summary, taxonomy, by = "ASVid")
Ratio_summary <- Ratio_summary[order(Host, Compartment, Prediction, Role)]
fwrite(Ratio_summary, file = "LotusHordeum_Askov_prediction_ratios_summary.csv")

Ratio_amount <- rbindlist(Ratio_amount_list)
Ratio_amount[,":="(
  Accuracy = paste0(round(Accuracy*100, 1), "%"),
  Compartment = relevel(Compartment, ref = "Rhizosphere")
)]
Ratio_amount <- Ratio_amount[order(Host, Compartment)]
fwrite(Ratio_amount, "Pred_accuracy_summary.csv")

Res[,.(Accuracy = mean(Obs == Pred)), .(Host, Compartment)]
fwrite(Res, "Prediction_results.csv")
