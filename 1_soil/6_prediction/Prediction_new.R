# Load packages and set colours --------------------------

pkg <- c(
  "data.table", "magrittr", "ggplot2", "ggh4x", "codacore", 
  "tensorflow", "ggpubr", "gridExtra", "cowplot", "ggtext"
)
for(pk in pkg){
  library(pk, character.only = T)
}

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define colours for the genotypes, soils, and bacterial orders of interest.
cols <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

colors <- c(NPK = "#6F944F", PK = "#B2563C", UF = "#3C7D82")

order_colors <- data.frame(
  group = c(
    "Burkholderiales", "Caulobacterales", "Flavobacteriales", "Micrococcales",
    "Pseudomonadales", "Rhizobiales", "Streptomycetales", "Sphingomonadales",
    "Pseudonocardiales", "Unknown", "Other"
  ),
  colors = c(
    "#645394", "#8e3563", "#05294a", "#44AAAA",
    "#88CCAA", "#fdbb6b", "#fed5a4", "lightyellow",
    "#95bb72", "grey", "lightgrey"
  )
)

# Load the input data --------------------------------------
## Lotus
ASV_table_Lotus <- fread("../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv")
colnames(ASV_table_Lotus)[1] <- "ASVid"

meta_data_Lotus <- fread(
  "../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", drop = c(5,7,8)
)
meta_data_Lotus <- meta_data_Lotus[Compartment != "Nodules"]

taxonomy_Lotus <- fread("../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv")

## Hordeum
ASV_table_Hordeum <- fread("../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv")
colnames(ASV_table_Hordeum)[1] <- "ASVid"

meta_data_Hordeum <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", drop = c(2,3,8)
)
meta_data_Hordeum <- meta_data_Hordeum[Compartment != "Soil"]

taxonomy_Hordeum <- fread("../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv")

# Clean up the taxonomy file layouts
taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxonomy_Lotus[,c(taxa_levels):= tstrsplit(Taxon, "; ", fill = "u__Unknown")]
taxonomy_Lotus[,c(taxa_levels):= lapply(.SD, substr, 4, 1000), .SDcols = taxa_levels]
taxonomy_Lotus[Kingdom == "ssigned", Kingdom:= "Unassigned"]
taxonomy_Lotus[,Taxon:=NULL]
setcolorder(taxonomy_Lotus, c("Feature ID", taxa_levels, "Confidence"))

taxonomy_Hordeum[,c(taxa_levels):= tstrsplit(Taxon, "; ", fill = "u__Unknown")]
taxonomy_Hordeum[,c(taxa_levels):= lapply(.SD, substr, 4, 1000), .SDcols = taxa_levels]
taxonomy_Hordeum[,Taxon:=NULL]
setcolorder(taxonomy_Hordeum, c("Feature ID", taxa_levels, "Confidence"))

# Set up list for easy retrieval later
ASV_table <- list(Lotus = ASV_table_Lotus, Hordeum = ASV_table_Hordeum)
meta_data <- rbind(meta_data_Lotus, meta_data_Hordeum)
taxonomy <- rbind(taxonomy_Lotus, taxonomy_Hordeum)
setnames(taxonomy, "Feature ID", "ASVid")

# Set up RA tables
ASV_table_Lotus_RA <- data.table(
  ASVid = ASV_table_Lotus$ASVid,
  t(t(ASV_table_Lotus[,-1])/colSums(ASV_table_Lotus[,-1]))
)

ASV_table_Hordeum_RA <- data.table(
  ASVid = ASV_table_Hordeum$ASVid,
  t(t(ASV_table_Hordeum[,-1])/colSums(ASV_table_Hordeum[,-1]))
)

ASV_table_RA <- list(Lotus = ASV_table_Lotus_RA, Hordeum = ASV_table_Hordeum_RA)

Opt <- expand.grid(Host = c("Lotus", "Hordeum"), Compartment = c("Root", "Rhizosphere"))

# Functions for the nested classification model --------------------------------
Nested_model <- function(OTU_table, y, highest_order, seed, lambda, overlap = T){
  Highest_order_resp <- ifelse(y == highest_order, highest_order, paste("Not", highest_order))
  Highest_order_resp <- relevel(factor(Highest_order_resp), ref = highest_order)
  
  OTU_table_T <- transpose(OTU_table, make.names = "ASVid", keep.names = "SampleID")

  set.seed(seed)
  tf$random$set_seed(seed)
  cc1 <- codacore(
    x = data.frame(OTU_table_T[,-1], row.names = OTU_table_T$SampleID),
    y = Highest_order_resp,
    objective = "binary classification",
    logRatioType = "SLR",
    lambda = lambda,
    overlap = overlap
  )
  
  Not_highest_idx <- which(y != highest_order)
  y2 <- y[Not_highest_idx]
  OTU_table_lower_T <- OTU_table_T[Not_highest_idx]
  
  cc2 <- codacore(
    x = data.frame(OTU_table_lower_T[,-1], row.names = OTU_table_lower_T$SampleID),
    y = y2,
    objective = "binary classification",
    logRatioType = "SLR",
    lambda = lambda,
    overlap = overlap
  )
  
  return( list(Highest = cc1, Lowest = cc2) )
}

inverse_logit <- function(x) exp(x)/(1+exp(x))
Nested_model_pred <- function(Model, eval_data, thresh = 0.5){
  eval_data_T <- transpose(eval_data, keep.names = "SampleID", make.names = "ASVid")

  pred <- predict(
    Model$Highest,
    data.frame(eval_data_T[,-1], row.names = eval_data_T$SampleID),
    logits = F
  )
  prob <- inverse_logit(pred)
  prob[is.nan(prob)] <- 1
  pred <- ifelse(prob < thresh, "NPK", "Not NPK")
  
  lower_idx <- which(pred == "Not NPK")
  x_lower <- eval_data_T[lower_idx]
  pred2 <- predict(
    Model$Lowest,
    data.frame(x_lower[,-1], row.names = x_lower$SampleID),
    logits = F
  )
  prob2 <- inverse_logit(pred2)
  prob2[is.nan(prob2)] <- 1
  pred2 <- ifelse(prob2 < thresh, "PK", "UF")
  
  pred[lower_idx] <- pred2
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

  ASV_table_sub <- ASV_table[[Opt$Host[i]]]
  meta_data_sub <- meta_data[
    Plant == Opt$Host[i] & Compartment == Opt$Compartment[i]
  ]
  sample_sub <- meta_data_sub$SampleID
  ASV_table_sub <- ASV_table_sub[,c("ASVid", sample_sub), with = F]
  
  # Filter out ASVs present in less than 10% of samples
  pres_idx <- which(rowMeans(ASV_table_sub[,-1] != 0) >= 0.1)
  ASV_table_sub <- ASV_table_sub[pres_idx]

  # Splitting into test and training data
  set.seed(1700294030)
  test_samples <- meta_data_sub[
    ,.(SampleID = sample(SampleID, 2, replace = FALSE)), by = .(Genotype, Soil)
  ]$SampleID
  training_samples <- setdiff(meta_data_sub$SampleID, test_samples)

  x_train <- ASV_table_sub[,c("ASVid", training_samples), with = F]
  x_test <- ASV_table_sub[,c("ASVid", test_samples), with = F]

  y_train <- meta_data_sub[SampleID %in% training_samples, Soil]
  y_test <- meta_data_sub[SampleID %in% test_samples, Soil]

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
  dt <- merge(dt, meta_data_sub[,c(1,3)], by = "SampleID")
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
  NPK_ASVs_in_Ratios <- apply(Mat_temp, 1, function(x) paste(x[1], x[2], sep = "/"))
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
  PK_ASVs_in_Ratios <- apply(Mat_temp, 1, function(x) paste(x[1], x[2], sep = "/"))
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
  Ratio_summary_list[[i]] <- Ratio_summary_temp

  # Summarising RA and taxonomic information for ASVs used for estimation
  ASV_pred <- unique(Ratio_summary_temp$ASVid)
  ASV_table_RA_sub <- ASV_table_RA[[Opt$Host[i]]]
  ASV_table_RA_sub <- ASV_table_RA_sub[ASVid %in% ASV_pred]
  ASV_table_RA_sub <- ASV_table_RA_sub[,c("ASVid", sample_sub), with = F]
  ASV_table_RA_sub <- merge(taxonomy[,c(1, 5)], ASV_table_RA_sub, by = "ASVid")
  Order_table_RA_sub <- ASV_table_RA_sub[,lapply(.SD, mean), Order, .SDcols = sample_sub]
  Order_table_RA_sub_T <- transpose(
    Order_table_RA_sub, make.names = "Order", keep.names = "SampleID"
  )
  Order_table_RA_sub_T <- merge(meta_data_sub, Order_table_RA_sub_T, by = "SampleID")
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

Res <- rbindlist(Res_list)
Ratio_summary <- rbindlist(Ratio_summary_list)
Ratio_summary <- merge(Ratio_summary, taxonomy, by = "ASVid")
fwrite(Ratio_summary, file = "LotusHordeum_Askov_prediction_ratios_summary.csv")

Ratio_amount <- rbindlist(Ratio_amount_list)

Res[,.(Accuracy = mean(Obs == Pred)), .(Host, Compartment)]

# Visualising results ------------------------------------------------

# Prediction results
Res[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Obs = factor(Obs, levels = c("NPK", "PK", "UF")),
  Pred = factor(Pred, levels = c("UF", "PK", "NPK"))
)]
Res[,Prediction := Obs == Pred]
genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

ggplot(data = Res) +
  geom_count(aes(x = Obs, y = Pred), color = "lightgrey")+
  scale_size_continuous(range=c(1.5,15)) +
  geom_jitter(data = Res[Prediction == T],
              aes(x = Obs, y = Pred, fill = Genotype),
              position = position_jitter(width = 0.35, height = 0.35, seed = 1),
              shape = 21, stroke = 0.25) +
  geom_jitter(data = Res[Prediction == F],
              aes(x = Obs, y = Pred, fill = Genotype),
              position = position_jitter(width = 0.15, height = 0.15, seed = 1),
              shape = 21, stroke = 0.25) +
  scale_shape_manual(values = c(21, 21))+
  theme_bw() +
  labs(x = "Observed", y = "Predicted")+
  guides(size = "none", fill = guide_legend(override.aes = list(size=3))) +
  facet_grid(Host ~ Compartment) +
  scale_fill_manual(values = cols, labels = genotype_labels_legend) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -8),
        strip.background = element_rect(colour = NA),
        axis.title.y = element_text(size = 8, family = "Helvetica"),
        axis.title.x = element_text(size = 8, family = "Helvetica"),
        axis.text.y = element_text(size = 8, family = "Helvetica",
                                   colour = "black"),
        axis.text.x = element_text(size = 8, family = "Helvetica",
                                   colour = "black"),
        legend.text = element_markdown(size = 8, family = "Helvetica",
                                       colour = "black"),
        legend.title = element_text(size = 8, family = "Helvetica",
                                    colour = "black"),
        strip.text = element_text(size = 8, family = "Helvetica",
                                  face = "bold"),
        legend.key.size = unit(5, "mm"))+
  NULL -> g1; g1

# Taxonomic composition
tax_summary[!(Order %in% order_colors$group), Order:= "Other"]
ggplot(data = tax_summary, aes(x = Soil, y = RA, fill = Order))+
  geom_bar(stat = "identity", position = "stack", linewidth = 0.1) +
  facet_wrap2(vars(Plant, Compartment), strip = strip_nested(), nrow = 1)+
  scale_fill_manual(values = order_colors$colors, breaks = order_colors$group,
                    name = "Bacterial orders")+
  labs(x = NULL, y = "Cumulative Mean Relative Abundance")+
  guides(fill = guide_legend(nrow = 4, title.position = "top"))+
  scale_y_continuous(expand = c(.001, .001), limits = c(0, 0.2))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.box.margin = margin(0, 0, 0, -15),
        legend.margin = margin(t = -5),
        strip.background = element_rect(colour = NA),
        legend.key.size = unit(0.25, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'),
        legend.justification = c(0.75, 0),
        axis.title.y = element_text(size = 8, family = "Helvetica"),
        axis.text.y = element_text(size = 8, family = "Helvetica",
                                   colour = "black"),
        axis.text.x = element_text(size = 8, family = "Helvetica",
                                   colour = "black"),
        legend.text = element_text(size = 8, family = "Helvetica"),
        legend.title = element_text(size = 8, family = "Helvetica"),
        strip.text = element_text(size = 8, family = "Helvetica",
                                  face = "bold"))+
  NULL -> g2; g2
