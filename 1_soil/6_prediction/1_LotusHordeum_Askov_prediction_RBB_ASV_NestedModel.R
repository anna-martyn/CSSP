# Clean up.
options(warn=-1)
rm(list=ls())

# Load packages and set colours ------------------------------------------------

# Load all required packages.
pkg <- c("data.table", "magrittr", "ggplot2", "RColorBrewer", "ggh4x", "dplyr", "tidyr",
         "codacore", "tensorflow", "ggpubr", "gridExtra", "cowplot", "ggtext")
for(pk in pkg){
  library(pk, character.only = T)
}

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define colours for the genotypes, soils, and bacterial orders of interest.
cols <- c("WT" = "#A9C289",
          "symrk" = "#FEDA8B",
          "ccamk" = "#FDB366",
          "nsp1" = "#C0E4EF",
          "nsp2" = "#6EA6CD")

colors <- c(NPK = "#6F944F", PK = "#B2563C", UF = "#3C7D82")

order_colors <- data.frame(group=c("Burkholderiales","Caulobacterales",
                                   "Flavobacteriales", "Micrococcales",
                                   "Pseudomonadales", "Rhizobiales",
                                   "Streptomycetales", "Unknown"),
                           colors=c("#645394", "#8e3563", 
                                    "#05294a", "#44AAAA",
                                    "#95bb72", "#fdbb6b",
                                    "#fed5a4", "grey") )

# Load the input data ----------------------------------------------------------
## Lotus.
ASV_table_Lotus <- fread("../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv")
colnames(ASV_table_Lotus)[1] <- "ASVid"

meta_data_Lotus <- fread(
  "../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", drop = c(5,7,8)
)
# meta_data_Lotus <- fread("../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt")
meta_data_Lotus <- meta_data_Lotus[Compartment != "Nodules"]

taxonomy_Lotus <- fread("../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv", sep="\t", header=TRUE, fill=TRUE)

## Hordeum.
ASV_table_Hordeum <- fread("../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv")
colnames(ASV_table_Hordeum)[1] <- "ASVid"

meta_data_Hordeum <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", drop = c(2,3,8)
)
# meta_data_Hordeum <- fread("../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt")
meta_data_Hordeum <- meta_data_Hordeum[Compartment != "Soil"]

taxonomy_Hordeum <- fread("../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv", sep="\t", header=TRUE, fill=TRUE)

# Clean up the taxonomy file layouts.
rename_tax <- function(tax_table){
  colnames(tax_table)[colnames(tax_table) == "Feature ID"] <- "ASVid"
  tax_table %>%
    separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
             sep = "; ", fill = "right") %>%
    mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
    replace(is.na(.), "Unknown") %>%
    select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)
}

taxonomy_Lotus <- rename_tax(taxonomy_Lotus)
taxonomy_Hordeum <- rename_tax(taxonomy_Hordeum)

# taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# taxonomy_Lotus[,c(taxa_levels):= tstrsplit(Taxon, "; ", fill = "u__Unknown")]
# taxonomy_Lotus[,c(taxa_levels):= lapply(.SD, substr, 4, 1000), .SDcols = taxa_levels]
# taxonomy_Lotus[Kingdom == "ssigned", Kingdom:= "Unassigned"]
# taxonomy_Lotus[,Taxon:=NULL]
# setcolorder(taxonomy_Lotus, c("Feature ID", taxa_levels, "Confidence"))

# taxonomy_Hordeum[,c(taxa_levels):= tstrsplit(Taxon, "; ", fill = "u__Unknown")]
# taxonomy_Hordeum[,c(taxa_levels):= lapply(.SD, substr, 4, 1000), .SDcols = taxa_levels]
# taxonomy_Hordeum[,Taxon:=NULL]
# setcolorder(taxonomy_Hordeum, c("Feature ID", taxa_levels, "Confidence"))


# Combine both datasets.
ASV_table <- list(Lotus = ASV_table_Lotus,
                  Hordeum = ASV_table_Hordeum)

meta_data <- list(Lotus = meta_data_Lotus,
                  Hordeum = meta_data_Hordeum)

taxonomy <- list(Lotus = taxonomy_Lotus,
                 Hordeum = taxonomy_Hordeum)

ASV_table_t <- lapply(ASV_table,
                      function(x) transpose(x,
                                            keep.names = "SampleID",
                                            make.names = "ASVid")
                      )

Full_data <- lapply(1:2,
                    function(i) merge(meta_data[[i]],
                                      ASV_table_t[[i]],
                                      by = "SampleID")
                    )
names(Full_data) <- c("Lotus", "Hordeum")

Opt <- expand.grid(Host = c("Lotus", "Hordeum"),
                   Compartment = c("Root", "Rhizosphere"))

# Set up tables for results ----------------------------------------------------
Res <- data.table(Host = NA,
                  Compartment = NA,
                  Accuracy = NA,
                  RA = NA,
                  PK_ratios = NA,
                  NPK_ratios = NA)[-1]

res2 <- data.table(SampleID = NA, Obs = NA, Pred = NA,
                   Host = NA, Compartment = NA)[-1]

Barplot_data <- data.table(Soil = NA,
                           Order = NA,
                           Abundance = NA,
                           Plant = NA,
                           Compartment = NA)[-1]

SLR_dt <- data.table(Soil = NA,
                     SLR = NA,
                     Pred_type = NA,
                     Plant = NA,
                     Compartment = NA)[-1]

ratio_summary <- data.table(
  ASV_ID = NA,
  Role = NA,
  Prediction = NA,
  Compartment = NA,
  Host = NA
)

# Functions for the nested classification model --------------------------------
Nested_model_three_categories <- function(OTU_table, y, highest_order,
                                          seed, lambda, overlap = T){
  Highest_order_resp <- ifelse(y == highest_order, 
                               highest_order, 
                               paste("Not", highest_order))
  Highest_order_resp <- relevel(factor(Highest_order_resp),
                                ref = highest_order)
  
  set.seed(seed)
  tf$random$set_seed(seed)
  cc1 <- codacore(x = OTU_table,
                  y = Highest_order_resp,
                  objective = "binary classification",
                  logRatioType = "SLR",
                  lambda = lambda,
                  overlap = overlap)
  
  Not_highest_idx <- y != highest_order
  y2 <- y[Not_highest_idx]
  OTU_table_lower <- OTU_table[Not_highest_idx,]
  
  cc2 <- codacore(x = OTU_table_lower,
                  y = y2,
                  objective = "binary classification",
                  logRatioType = "SLR",
                  lambda = lambda,
                  overlap = overlap)
  
  return( list(Highest = cc1, Lowest = cc2) )
}

Nested_model_three_categories_predict <- function(Model, eval_data){
  pred <- predict(Model$Highest, eval_data, logits = F)
  prob <- boot::inv.logit(pred)
  pred <- ifelse(prob<0.5, "NPK", "Not NPK")
  
  lower_idx <- which(pred == "Not NPK")
  x_lower <- eval_data[lower_idx,]
  pred2 <- predict(Model$Lowest, x_lower, logits = F)
  prob2 <- boot::inv.logit(pred2)
  pred2 <- ifelse(prob2<0.5, "PK", "UF")
  
  pred[lower_idx] <- pred2
  return(pred)
}

# Run the prediction analysis --------------------------------------------------
for(i in 1:nrow(Opt)){
  run <- paste(as.matrix(Opt)[i,], collapse = "_")
  Current_plant <- as.matrix(Opt)[i,1]
  cat(run, "\r")
  
  data_subset <- Full_data[[Current_plant]][Compartment == Opt$Compartment[i]]
  
  OTU_present <- apply(data_subset[,-(1:5)], 2, function(x) mean(x != 0) >= 0.1)
  OTU_keep <- colnames(data_subset)[-(1:5)][OTU_present]
  
  set.seed(1700294030)
  temp <- setDT(data_subset)[,sample(SampleID, 2, replace = FALSE),
                             by = .(Genotype, Soil)]
  test_data_samples <- temp$V1
  train_data_samples <- data_subset[!(SampleID %in% temp$V1 ), SampleID]
  
  train_data <- data_subset[SampleID %in% train_data_samples]
  test_data <- data_subset[SampleID %in% test_data_samples]
  
  x_train <- data.frame(train_data[,..OTU_keep],
                        row.names = train_data$SampleID)
  x_test <- data.frame(test_data[,..OTU_keep],
                       row.names = test_data$SampleID)
  
  colnames(x_train) <- gsub("X", "", colnames(x_train))
  colnames(x_test) <- gsub("X", "", colnames(x_test))
  
  y_train <- train_data$Soil
  y_test <- test_data$Soil
  
  tf$random$set_seed(1700294030)
  cc <- Nested_model_three_categories(OTU_table = x_train,
                                      y = y_train,
                                      highest_order = "NPK",
                                      seed = 1700299435,
                                      lambda = 1)
  
  pred <- Nested_model_three_categories_predict(cc, x_test)
  eval <- mean(y_test == pred)
  
  dt <- data.table(SampleID = names(pred), Obs = y_test, Pred = pred)
  dt <- dt[order(Obs)]
  dt[,":="(Host = Opt$Host[i], Compartment = Opt$Compartment[i])]
  res2 <- rbind(res2, dt)
  
  num_vec_NPK <- lapply(cc$Highest$ensemble,
                        function(x) colnames(x_train)[x$hard$numerator])
  denom_vec_NPK <- lapply(cc$Highest$ensemble,
                          function(x) colnames(x_train)[x$hard$denominator])
  
  num_NPK_dt <- data.table(
    ASV_ID = unlist(num_vec_NPK),
    Role = paste0(
      "Numerator",
      rep(1:length(num_vec_NPK), lapply(num_vec_NPK, length))
    ),
    Prediction = "NPK vs. non-NPK",
    Compartment = as.character(Opt$Compartment[i]),
    Host = as.character(Opt$Host[i])
  )
  
  denom_NPK_dt <- data.table(
    ASV_ID = unlist(denom_vec_NPK),
    Role = paste0(
      "Denominator",
      rep(1:length(denom_vec_NPK), lapply(denom_vec_NPK, length))
    ),
    Prediction = "NPK vs. non-NPK",
    Compartment = as.character(Opt$Compartment[i]),
    Host = as.character(Opt$Host[i])
  )
  
  num_NPK <- lapply(num_vec_NPK, function(x) paste(x, collapse = "+"))
  denom_NPK <- lapply(denom_vec_NPK, function(x) paste(x, collapse = "+"))
  
  # Remove duplicate ratios for NPK.
  ratios_names_NPK <- paste(unlist(num_NPK), unlist(denom_NPK), sep = "/")
  NPK_idx <- match(unique(ratios_names_NPK), ratios_names_NPK)

  num_vec_NPK <- num_vec_NPK[NPK_idx]
  denom_vec_NPK <- denom_vec_NPK[NPK_idx]
  num_NPK <- num_NPK[NPK_idx]
  denom_NPK <- denom_NPK[NPK_idx]
  
  num_vec_PK <- lapply(cc$Lowest$ensemble,
                       function(x) colnames(x_train)[x$hard$numerator])
  denom_vec_PK <- lapply(cc$Lowest$ensemble,
                         function(x) colnames(x_train)[x$hard$denominator])
  
  num_PK_dt <- data.table(
    ASV_ID = unlist(num_vec_PK),
    Role = paste0(
      "Numerator",
      rep(1:length(num_vec_PK), lapply(num_vec_PK, length))
    ),
    Prediction = "PK vs. UF",
    Compartment = as.character(Opt$Compartment[i]),
    Host = as.character(Opt$Host[i])
  )
  
  denom_PK_dt <- data.table(
    ASV_ID = unlist(denom_vec_PK),
    Role = paste0(
      "Denominator",
      rep(1:length(denom_vec_PK), lapply(denom_vec_PK, length))
    ),
    Prediction = "PK vs. UF",
    Compartment = as.character(Opt$Compartment[i]),
    Host = as.character(Opt$Host[i])
  )
  
  ratio_summary <- rbind(
    ratio_summary, num_NPK_dt, denom_NPK_dt, num_PK_dt, denom_PK_dt
  )
  
  num_PK <- lapply(num_vec_PK, function(x) paste(x, collapse = "+"))
  denom_PK <- lapply(denom_vec_PK, function(x) paste(x, collapse = "+"))
  
  # Remove duplicate ratios for PK.
  ratios_names_PK <- paste(unlist(num_PK), unlist(denom_PK), sep = "/")
  PK_idx <- match(unique(ratios_names_PK), ratios_names_PK)

  num_vec_PK <- num_vec_PK[PK_idx]
  denom_vec_PK <- denom_vec_PK[PK_idx]
  num_PK <- num_PK[PK_idx]
  denom_PK <- denom_PK[PK_idx]
  
  r1 <- length(num_NPK)
  r2 <- length(num_PK)
  
  all_pred_ASV <- unique(c(unlist(num_vec_NPK), unlist(denom_vec_NPK),
                           unlist(num_vec_PK), unlist(denom_vec_PK)))
  # pred_ASV_tax <- taxonomy[[Current_plant]][Feature %in% all_pred_ASV]
  
  Taxonomy_df <- data.frame(as.data.frame(taxonomy[[Current_plant]][,-1]),
                            row.names = taxonomy[[Current_plant]]$ASV)
  
  full_data <- rbind(x_train, x_test)
  Sample_info <- meta_data[[Current_plant]][SampleID %in% rownames(full_data),
                                            c("SampleID", "Soil"), with = F]
  Sample_info <- Sample_info[match(rownames(full_data), SampleID)]
  
  Used_ASV_info <- cbind(
    Sample_info[,-1], full_data[all_pred_ASV]/rowSums(full_data)
  )
  
  num1 <- rowSums(as.data.frame(full_data[,num_vec_NPK[[1]]]))
  num1[num1 == 0] <- 1 
  denom1 <- rowSums(as.data.frame(full_data[,denom_vec_NPK[[1]]]))
  denom1[denom1 == 0] <- 1
  SLR <- log(num1/denom1)
  
  num1 <- rowSums(as.data.frame(full_data[,num_vec_PK[[1]]]))
  num1[num1 == 0] <- 1 
  denom1 <- rowSums(as.data.frame(full_data[,denom_vec_PK[[1]]]))
  denom1[denom1 == 0] <- 1
  SLR_PK <- log(num1/denom1)
  
  SLR_dt_temp <- data.table(Soil = rep(Sample_info$Soil, 2),
                            SLR = c(SLR, SLR_PK),
                            Pred_type = rep(c("NPK", "PK"), each = length(SLR)),
                            Plant = Opt$Host[i],
                            Compartment = Opt$Compartment[i])
  
  SLR_dt <- rbind(SLR_dt, SLR_dt_temp)
  
  Used_ASV_info2 <- melt(Used_ASV_info,
                         value.name = "Abundance",
                         variable.name = "ASV")
  
  rbind(
    Used_ASV_info2[Soil == "NPK",.(Abundance = mean(Abundance)), .(Soil, ASV)],
    Used_ASV_info2[Soil == "PK",.(Abundance = mean(Abundance)), .(Soil, ASV)],
    Used_ASV_info2[Soil == "UF",.(Abundance = mean(Abundance)), .(Soil,ASV)]
  ) -> Used_ASV_info2
  
  Used_ASV_info2[,Order:=Taxonomy_df[ASV,]$Order]
  rbind(
    Used_ASV_info2[Soil == "NPK",.(Abundance = sum(Abundance)), .(Soil, Order)],
    Used_ASV_info2[Soil == "PK",.(Abundance = sum(Abundance)), .(Soil, Order)],
    Used_ASV_info2[Soil == "UF",.(Abundance = sum(Abundance)), .(Soil,Order)]
  ) -> Used_ASV_info2
  
  Used_ASV_info2[,":="(Plant = Opt$Host[i], Compartment = Opt$Compartment[i])]
  
  Barplot_data <- rbind(Barplot_data, Used_ASV_info2)
  
  acc_RA <- mean(rowSums(full_data[,all_pred_ASV])/rowSums(full_data))
  
  max11 <- lapply(cc$Highest$ensemble, function(x) sum(x$hard$denominator)) %>%
    unlist() %>% max()
  max12 <- lapply(cc$Lowest$ensemble, function(x) sum(x$hard$denominator)) %>% 
    unlist() %>% max()
  max1 <- max(c(max11, max12))
  max21 <- lapply(cc$Highest$ensemble, function(x) sum(x$hard$numerator)) %>%
    unlist() %>% max()
  max22 <- lapply(cc$Lowest$ensemble, function(x) sum(x$hard$numerator)) %>% 
    unlist() %>% max()
  max2 <- max(c(max21, max22))
  
  N_ratios_PK <- length(num_PK)
  N_ASVs_in_num <- unlist(lapply(num_vec_PK, length))
  N_ASVs_in_denom <- unlist(lapply(denom_vec_PK, length))
  N_ASVs_in_ratios <- N_ASVs_in_num + N_ASVs_in_denom
  N_ASVs_in_ratios_PK_as_char <- paste(N_ASVs_in_ratios, collapse = ", ")
  PK_ratios_info <- paste0(N_ratios_PK, " (", N_ASVs_in_ratios_PK_as_char, ")")
  
  N_ratios_NPK <- length(num_NPK)
  N_ASVs_in_num <- unlist(lapply(num_vec_NPK, length))
  N_ASVs_in_denom <- unlist(lapply(denom_vec_NPK, length))
  N_ASVs_in_ratios <- N_ASVs_in_num + N_ASVs_in_denom
  N_ASVs_in_ratios_NPK_as_char <- paste(N_ASVs_in_ratios, collapse = ", ")
  NPK_ratios_info <- paste0(N_ratios_NPK, " (",
                            N_ASVs_in_ratios_NPK_as_char, ")")
  
  dt <- data.table(Host = Opt$Host[i],
                   Compartment = Opt$Compartment[i],
                   Accuracy = round(eval*100,1),
                   RA = round(acc_RA*100,1),
                   PK_ratios = PK_ratios_info,
                   NPK_ratios = NPK_ratios_info)
  
  Res <- rbind(Res, dt)
}

full_taxonomy <- rbind(taxonomy_Hordeum, taxonomy_Lotus)
colnames(full_taxonomy)[1] <- "ASV_ID"
ratio_summary <- merge(
  ratio_summary, full_taxonomy, by = "ASV_ID", all.x = T
)
ratio_summary <- ratio_summary[-1]
ratio_summary <- ratio_summary[order(Host, Compartment, Prediction, Role)]

fwrite(ratio_summary, "LotusHordeum_Askov_prediction_ratios_summary.csv")

# Visualizing results ----------------------------------------------------------
meta_data_full <- rbind(meta_data_Hordeum, meta_data_Lotus)
meta_data_full <- data.frame(meta_data_full[,-1],
                             row.names = meta_data_full$SampleID)
res2[,Genotype:=meta_data_full[res2$SampleID,"Genotype"]]

res2[,":="(Host = factor(Host, levels = c("Lotus", "Hordeum")),
           Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
           Obs = factor(Obs, levels = c("NPK", "PK", "UF")),
           Pred = factor(Pred, levels = c("UF", "PK", "NPK")))]

res2[,Prediction := Obs == Pred]

genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

res2[,Genotype:=factor(Genotype, 
                       levels = c("WT", "symrk", "ccamk", "nsp1","nsp2"))]

ggplot(data = res2) +
  geom_count(aes(x = Obs, y = Pred), color = "lightgrey")+
  scale_size_continuous(range=c(1.5,15)) +
  geom_jitter(data = res2[Prediction == T],
              aes(x = Obs, y = Pred, fill = Genotype),
              position = position_jitter(width = 0.35, height = 0.35, seed = 1),
              shape = 21, stroke = 0.25) +
  geom_jitter(data = res2[Prediction == F],
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

Barplot_data[,m:=paste(Plant, Compartment, sep = " \n")]
Barplot_data[,m:=factor(m, levels = c("Lotus \nRhizosphere",
                                      "Hordeum \nRhizosphere",
                                      "Lotus \nRoot",
                                      "Hordeum \nRoot"))]
SLR_dt[,m:=paste(Plant, Compartment)]
SLR_dt[Pred_type == "NPK", Pred_type:="NPK vs non-NPK"]
SLR_dt[Pred_type == "PK", Pred_type:="PK vs UF"]

reorder <- c("Burkholderiales", "Caulobacterales",
             "Flavobacteriales", "Pseudomonadales",
             "Rhizobiales", "Streptomycetales",
             "Micrococcales", "Unknown")
order_colors <- order_colors[match(reorder, order_colors$group),]
Barplot_data[,Compartment:=factor(Compartment, 
                                  levels = c("Rhizosphere", "Root"))]
ggplot(data = Barplot_data, aes(x = Soil, y = Abundance, fill = Order))+
  geom_bar(stat = "identity", position = "stack", linewidth = 0.1) +
  facet_wrap2(vars(Plant, Compartment), strip = strip_nested(), nrow = 1)+
  scale_fill_manual(values = order_colors$colors, breaks = order_colors$group,
                    name = "Bacterial orders")+
  labs(x = NULL, y = "Cumulative Mean Relative Abundance")+
  guides(fill = guide_legend(nrow = 3, title.position = "top"))+
  scale_y_continuous(expand = c(.001, .001), limits = c(0, 0.44))+
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

ggplot() + theme_void() -> blank

R <- Res[,c(1,2,6,5,3)]
R[,Accuracy:=paste0(Accuracy, "%")]
R[Host == "Barley", Host:="Hordeum"]
colnames(R) <- gsub("_", " ", colnames(R))
R <- R[c(3:4,1:2)]

# IMPORTANT! 
# Run Flowchart.R BEFORE running the code below to produce the correct figure!
# Correct order of events:
#  1. run line 1-535 of this script
#  2, Run flowchart.R
#  3. run the rest of this scripts.

tg <- tableGrob(R, theme = ttheme_default(base_size = 8), rows = NULL)

gg1 <- plot_grid(grob, tg,
                 labels = c("A", "C"),
                 rel_heights = c(0.8, 0.2),
                 ncol = 1)
gg2 <- ggarrange(g1, g2, labels = c("B", "D"), ncol = 1, heights = c(0.5, 0.5))
gg <- ggarrange(gg1, gg2, ncol = 2)

ggsave(filename = "Figure3_Askov_prediction.pdf", plot = gg,
       width = 210, height = 200, units = "mm")

fwrite(Res, "LotusHordeum_Askov_prediction_summary.csv")
