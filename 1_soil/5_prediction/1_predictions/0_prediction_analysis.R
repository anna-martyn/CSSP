# Seup ------------------------------------------------------------------------
# Loading packages
pkg <- c("data.table", "codacore","tensorflow")
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
ratio_summary_list <- list()
res_list <- list()
ratio_amount_list <- list()
tax_summary_list <- list()
for(i in 1:nrow(Opt)){
  run <- paste(as.matrix(Opt)[i,], collapse = "_")
  cat(run, "\r")
  
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
    asv_table = x_train,
    y = y_train,
    seed = 1700299435,
    lambda = 1
  )

  pred <- Nested_model_pred(cc, x_test)
  eval <- mean(y_test == pred)
  
  # Summarising results
  pred_res <- data.table(SampleID = names(pred), Obs = y_test, Pred = pred)
  pred_res <- pred_res[order(Obs)]
  pred_res[,":="(Host = current_host, Compartment = current_comp)]
  pred_res <- merge(
    x = pred_res,
    y = design_sub[,c("SampleID","Genotype")],
    by = "SampleID"
  )
  res_list[[i]] <- pred_res

  # Extracting ASVs used in predictions
  ## Numerators model 1
  num_model1 <- lapply(
    cc$Model1$ensemble, function(x) x_train$ASVid[x$hard$numerator]
  )
  ## Denominators model 1
  denom_model1 <- lapply(
    cc$Model1$ensemble, function(x) x_train$ASVid[x$hard$denominator]
  )
  ## Numerators model 2
  num_model2 <- lapply(
    cc$Model2$ensemble, function(x) x_train$ASVid[x$hard$numerator]
  )
  ## Denominators model 2
  denom_model2 <- lapply(
    cc$Model2$ensemble, function(x) x_train$ASVid[x$hard$denominator]
  )

  # Counting the ASVs in ratios for table C
  count_mat_model1 <- matrix(
    c(
      unlist(lapply(num_model1, length)),
      unlist(lapply(denom_model1, length))
    ),
    ncol = 2
  )
  # asv_in_ratios_model1 <- apply(count_mat_model1, 1, function(x) paste(x[1], x[2], sep = "/"))
  asv_in_ratios_model1 <- rowSums(count_mat_model1)
  ratios_model1 <- paste0(
    length(num_model1), " (",
    paste(asv_in_ratios_model1, collapse = ", "), ")"
  )

  count_mat_model2 <- matrix(
    c(
      unlist(lapply(num_model2, length)),
      unlist(lapply(denom_model2, length))
    ),
    ncol = 2
  )
  # asv_in_ratios_model2 <- apply(count_mat_model2, 1, function(x) paste(x[1], x[2], sep = "/"))
  asv_in_ratios_model2 <- rowSums(count_mat_model2)
  ratios_model2 <- paste0(
    length(num_model2), " (",
    paste(asv_in_ratios_model2, collapse = ", "), ")"
  )

  ratio_amount_list[[i]] <- data.table(
    Host = current_host,
    Compartment = current_comp,
    "NPK ratios" = ratios_model1,
    "PK ratios" = ratios_model2,
    Accuracy = eval
  )

  # Constructing supplementary table
  # Summarising all ratios in both models
  ## Numerators in model 1
  model1_num_summary <- lapply(
    1:length(num_model1), function(i) data.table(
      ASVid = num_model1[[i]],
      Role = paste0("Numerator", i),
      Prediction = "NPK vs. non-NPK"
    )
  )
  ## Denominators in model 1
  model1_denom_summary <- lapply(
    1:length(denom_model1), function(i) data.table(
      ASVid = denom_model1[[i]],
      Role = paste0("Denominator", i),
      Prediction = "NPK vs. non-NPK"
    )
  )
  ## Numerators in model 2
  model2_num_summary <- lapply(
    1:length(num_model2), function(i) data.table(
      ASVid = num_model2[[i]],
      Role = paste0("Numerator", i),
      Prediction = "PK vs. UF"
    )
  )
  ## Denominators in model 2
  model2_denom_summary <- lapply(
    1:length(denom_model2), function(i) data.table(
      ASVid = denom_model2[[i]],
      Role = paste0("Denominator", i),
      Prediction = "PK vs. UF"
    )
  )
  # Collecting information on both numerators and denominators into a single data table
  both_models_ratios_list <- c(
    model1_num_summary, model1_denom_summary,
    model2_num_summary, model2_denom_summary
  )
  both_models_ratios <- rbindlist(both_models_ratios_list)
  both_models_ratios[,":="(Compartment = current_comp, Host = current_host)]
  
  # Vector containing ASVs involvled in at least one log-ratio used for predictions
  pred_asv <- unique(both_models_ratios$ASVid)
  
  # Subsetting ASV table with RAs to keep only samples from the current host and 
  # ASVs involved in log-ratios used for predictions
  asv_table_RA_sub <- asv_table_RA[[current_host]]
  asv_table_RA_sub <- asv_table_RA_sub[ASVid %in% pred_asv]
  
  # Transposing ASV table
  asv_table_RA_sub <- transpose(
    asv_table_RA_sub, keep.names = "SampleID", make.names = "ASVid"
  )
  
  # Adding information from metadata to merged ASV table
  # Only samples from the current compartment are kept, since design_sub only 
  # contains these samples, and only common samples are kept in merge, unless
  # all = TRUE is set
  asv_table_RA_sub <- merge(design_sub, asv_table_RA_sub, by = "SampleID")
  # Taking the mean RA for each prediction ASV within soils and across genotypes
  mean_RAs <- asv_table_RA_sub[,lapply(.SD, mean), Soil, .SDcols = pred_asv]
  mean_RAs <- transpose(mean_RAs, make.names = "Soil", keep.names = "ASVid")
  colnames(mean_RAs)[-1] <- paste(colnames(mean_RAs)[-1], "RA", sep = "_")

  # Merging ratio summaries with RA means
  both_models_ratios <- merge(both_models_ratios, mean_RAs)
  ratio_summary_list[[i]] <- both_models_ratios

  # Summarising RA and taxonomic information for ASVs used for prediction
  ## Subsetting ASV table with RAs, keeping only samples from current 
  ## host-compartment combination and ASVs used for prediction
  asv_table_RA_sub <- asv_table_RA[[current_host]]
  asv_table_RA_sub <- asv_table_RA_sub[ASVid %in% pred_asv]
  asv_table_RA_sub <- asv_table_RA_sub[,c("ASVid", sample_sub), with = F]
  
  ## Adding order information to ASV table
  asv_table_RA_sub <- merge(
    x = taxonomy[,c("ASVid", "Order")],
    y = asv_table_RA_sub,
    by = "ASVid"
  )
  
  ## Aggregating RAs at the order level
  orders_RA <- asv_table_RA_sub[,lapply(.SD, sum), Order, .SDcols = sample_sub]
  
  ## Tranposing order information and merging it with metadata
  orders_RA_transp <- transpose(
    orders_RA, make.names = "Order", keep.names = "SampleID"
  )
  orders_RA_transp <- merge(design_sub, orders_RA_transp, by = "SampleID")

  ## Converting to long form and saving information on orders
  order_info <- melt(
    orders_RA_transp,
    id.vars = 1:5,
    variable.name = "Order",
    value.name = "RA"
  )
  tax_summary_list[[i]] <- order_info[,.(RA = mean(RA)), .(Plant, Compartment, Soil, Order)]
}

# Collecting taxonomic information on ASVs used for predictions in all 
# plant-compartment combinations
tax_summary <- rbindlist(tax_summary_list)
tax_summary[,.(RA = sum(RA)), .(Plant, Compartment, Soil)]
fwrite(tax_summary, "1_tables/Pred_taxonomic_composition.csv")

# Collecting information on which ASVs are used for predictions in all 
# plant-compartment combinations
res <- rbindlist(res_list)
ratio_summary <- rbindlist(ratio_summary_list)
ratio_summary <- merge(ratio_summary, taxonomy, by = "ASVid")
ratio_summary <- ratio_summary[order(Host, Compartment, Prediction, Role)]
fwrite(
  ratio_summary,
  file = "1_tables/LotusHordeum_Askov_prediction_ratios_summary.csv"
)

# Collecting information on the amount of ASVs used for predictions in all 
# plant-compartment combinations and their accuracies
ratio_amount <- rbindlist(ratio_amount_list)
ratio_amount[,":="(
  Accuracy = paste0(round(Accuracy*100, 1), "%"),
  Compartment = relevel(Compartment, ref = "Rhizosphere")
)]
ratio_amount <- ratio_amount[order(Host, Compartment)]
fwrite(ratio_amount, "1_tables/Pred_accuracy_summary.csv")

# Collecting information on the observed vs. predicted soils in all
# samples in the test data
fwrite(res, "1_tables/Prediction_results.csv")

res[,.(Accuracy = mean(Obs == Pred)), .(Host, Compartment)]
