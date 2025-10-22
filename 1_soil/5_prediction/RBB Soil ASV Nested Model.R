# Load packages and set colours ------------------------------------------------
pkg <- c("data.table", "magrittr", "ggplot2", "RColorBrewer", "ggh4x",
         "codacore", "tensorflow", "ggpubr", "gridExtra", "cowplot")
for(pk in pkg){
  library(pk, character.only = T)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
                           colors=c("#645394", "#AA4488", 
                                    "#ffeeef", "#44AAAA",
                                    "#88CCAA", "#fdbb6b",
                                    "#ffd7b5", "grey") )

# Loading the data -------------------------------------------------------------
ASV_table_lotus <- fread("../1_data/1_Lotus/feature-table.tsv")
colnames(ASV_table_lotus)[1] <- "ASV_ID"
meta_data_dir <- paste("..", "1_data", "1_Lotus", 
                       "Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt",
                       sep = "/")
meta_data_lotus <- fread(meta_data_dir, drop = c(5,7:10))
taxonomy_lotus <- fread("../1_data/1_Lotus/taxonomy_clean.csv")

meta_data_lotus <- meta_data_lotus[Plant != "Soil"]
meta_data_lotus <- meta_data_lotus[Compartment != "Nodules"]
meta_data_lotus[Compartment == "Endosphere/Rhizoplane", Compartment:="Root"]
colnames(meta_data_lotus)[1] <- "Sample_ID"

ASV_tab_dir <- paste("..", "1_data", "2_Barley",
                     "BarleyCSSP_Askov_reseq_ASVtable_10_4.tsv", sep = "/")
ASV_table_barley <- fread(ASV_tab_dir)
colnames(ASV_table_barley)[1] <- "ASV_ID"
meta_data_dir <- paste("..", "1_data", "2_Barley", 
                       "BarleyCSSP_Askov_reseq_metadata.txt",
                       sep = "/")
meta_data_barley <- fread(meta_data_dir, drop = c(2:3,8))
taxonomy_barley <- fread("../1_data/2_Barley/taxonomy_clean.csv")
meta_data_barley <- meta_data_barley[Plant != "Soil"]
meta_data_barley[,Compartment:=fcase(Compartment == "rhizo", "Rhizosphere",
                                     Compartment == "endo", "Root")]

ASV_table <- list(Lotus = ASV_table_lotus,
                  Barley = ASV_table_barley)

meta_data <- list(Lotus = meta_data_lotus,
                  Barley = meta_data_barley)

taxonomy <- list(Lotus = taxonomy_lotus,
                 Barley = taxonomy_barley)

ASV_table_t <- lapply(ASV_table,
                      function(x) transpose(x,
                                            keep.names = "Sample_ID",
                                            make.names = "ASV_ID")
                      )

Full_data <- lapply(1:2,
                    function(i) merge(meta_data[[i]],
                                      ASV_table_t[[i]],
                                      by = "Sample_ID")
                    )
names(Full_data) <- c("Lotus", "Barley")

Opt <- expand.grid(Host = c("Lotus", "Barley"),
                   Compartment = c("Root", "Rhizosphere"))

# Settings up tables for results -----------------------------------------------
Res <- data.table(Host = NA,
                  Compartment = NA,
                  Accuracy = NA,
                  RA = NA,
                  PK_ratios = NA,
                  NPK_ratios = NA)[-1]

res2 <- data.table(Sample_ID = NA, Obs = NA, Pred = NA,
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

# Running the prediction analysis ----------------------------------------------
for(i in 1:nrow(Opt)){
  run <- paste(as.matrix(Opt)[i,], collapse = "_")
  Current_plant <- as.matrix(Opt)[i,1]
  cat(run, "\r")
  
  data_subset <- Full_data[[Current_plant]][Compartment == Opt$Compartment[i]]
  
  OTU_present <- apply(data_subset[,-(1:5)], 2, function(x) mean(x != 0) >= 0.1)
  OTU_keep <- colnames(data_subset)[-(1:5)][OTU_present]
  
  set.seed(1700294030)
  temp <- setDT(data_subset)[,sample(Sample_ID, 2, replace = FALSE),
                             by = .(Genotype, Soil)]
  test_data_samples <- temp$V1
  train_data_samples <- data_subset[!(Sample_ID %in% temp$V1 ), Sample_ID]
  
  train_data <- data_subset[Sample_ID %in% train_data_samples]
  test_data <- data_subset[Sample_ID %in% test_data_samples]
  
  x_train <- data.frame(train_data[,..OTU_keep],
                        row.names = train_data$Sample_ID)
  x_test <- data.frame(test_data[,..OTU_keep],
                       row.names = test_data$Sample_ID)
  
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
  
  dt <- data.table(Sample_ID = names(pred), Obs = y_test, Pred = pred)
  dt <- dt[order(Obs)]
  dt[,":="(Host = Opt$Host[i], Compartment = Opt$Compartment[i])]
  res2 <- rbind(res2, dt)
  
  num_vec_NPK <- lapply(cc$Highest$ensemble,
                        function(x) colnames(x_train)[x$hard$numerator])
  denom_vec_NPK <- lapply(cc$Highest$ensemble,
                          function(x) colnames(x_train)[x$hard$denominator])
  
  num_NPK <- lapply(num_vec_NPK, function(x) paste(x, collapse = "+"))
  denom_NPK <- lapply(denom_vec_NPK, function(x) paste(x, collapse = "+"))
  
  # Removing duplicate ratios for NPK
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
  
  num_PK <- lapply(num_vec_PK, function(x) paste(x, collapse = "+"))
  denom_PK <- lapply(denom_vec_PK, function(x) paste(x, collapse = "+"))
  
  # Removing duplicate ratios for PK
  ratios_names_PK <- paste(unlist(num_PK), unlist(denom_PK), sep = "/")
  PK_idx <- match(unique(ratios_names_PK), ratios_names_PK)

  num_vec_PK <- num_vec_PK[PK_idx]
  denom_vec_PK <- denom_vec_PK[PK_idx]
  num_PK <- num_PK[PK_idx]
  denom_PK <- denom_PK[PK_idx]
  
  r1 <- length(num_NPK)
  r2 <- length(num_PK)
  
  Ratio_information1 <- data.table("Ratio for NPK prediction" = c( rbind(paste0("Numerator", 1:r1),
                                                   paste0("Denominator", 1:r1)) ),
                                  ASV = c( rbind(unlist(num_NPK), unlist(denom_NPK)) ))
  
  Ratio_information2 <- data.table("Ratio for NPK prediction" = c( rbind(paste0("Numerator", 1:r2),
                                                    paste0("Denominator", 1:r2)) ),
                                   ASV = c( rbind(unlist(num_PK), unlist(denom_PK)) ))
  
  Ratio_information2 <- rbind(
    data.table("Ratio for NPK prediction" = "Ratio for PK prediction", ASV = ""),
    Ratio_information2
  )
  
  Ratio_information <- rbind(Ratio_information1, Ratio_information2)
  
  dir <- paste("Ratios/Ratio_information_", run, ".csv", sep = "")
  fwrite(Ratio_information, dir)
  
  all_pred_ASV <- unique(c(unlist(num_vec_NPK), unlist(denom_vec_NPK),
                           unlist(num_vec_PK), unlist(denom_vec_PK)))
  pred_ASV_tax <- taxonomy[[Current_plant]][Feature %in% all_pred_ASV]
  
  dir2 <- paste("Taxonomy/Pred_ASV_", run, ".csv", sep = "")
  fwrite(pred_ASV_tax, dir2)
  
  Taxonomy_df <- data.frame(as.data.frame(taxonomy[[Current_plant]][,-1]),
                            row.names = taxonomy[[Current_plant]]$ASV)
  
  full_data <- rbind(x_train, x_test)
  Sample_info <- meta_data[[Current_plant]][Sample_ID %in% rownames(full_data),
                                            c("Sample_ID", "Soil"), with = F]
  Sample_info <- Sample_info[match(rownames(full_data), Sample_ID)]
  
  Used_ASV_info <- cbind(Sample_info[,-1], full_data[all_pred_ASV]/rowSums(full_data))
  
  num1 <- rowSums(as.data.frame(full_data[,num_vec_NPK[[1]]])); num1[num1 == 0] <- 1 
  denom1 <- rowSums(as.data.frame(full_data[,denom_vec_NPK[[1]]])); denom1[denom1 == 0] <- 1
  SLR <- log(num1/denom1)
  
  num1 <- rowSums(as.data.frame(full_data[,num_vec_PK[[1]]])); num1[num1 == 0] <- 1 
  denom1 <- rowSums(as.data.frame(full_data[,denom_vec_PK[[1]]])); denom1[denom1 == 0] <- 1
  SLR_PK <- log(num1/denom1)
  
  SLR_dt_temp <- data.table(Soil = rep(Sample_info$Soil, 2),
                            SLR = c(SLR, SLR_PK),
                            Pred_type = rep(c("NPK", "PK"), each = length(SLR)),
                            Plant = Opt$Host[i],
                            Compartment = Opt$Compartment[i])
  
  SLR_dt <- rbind(SLR_dt, SLR_dt_temp)
  
  Used_ASV_info2 <- melt(Used_ASV_info, value.name = "Abundance", variable.name = "ASV")
  
  rbind(Used_ASV_info2[Soil == "NPK",.(Abundance = mean(Abundance)), .(Soil, ASV)],
        Used_ASV_info2[Soil == "PK",.(Abundance = mean(Abundance)), .(Soil, ASV)],
        Used_ASV_info2[Soil == "UF",.(Abundance = mean(Abundance)), .(Soil,ASV)]) -> Used_ASV_info2
  
  Used_ASV_info2[,Order:=Taxonomy_df[ASV,]$Order]
  rbind(Used_ASV_info2[Soil == "NPK",.(Abundance = sum(Abundance)), .(Soil, Order)],
        Used_ASV_info2[Soil == "PK",.(Abundance = sum(Abundance)), .(Soil, Order)],
        Used_ASV_info2[Soil == "UF",.(Abundance = sum(Abundance)), .(Soil,Order)]) -> Used_ASV_info2
  
  Used_ASV_info2[,":="(Plant = Opt$Host[i], Compartment = Opt$Compartment[i])]
  
  Barplot_data <- rbind(Barplot_data, Used_ASV_info2)
  
  acc_RA <- mean(rowSums(full_data[,all_pred_ASV])/rowSums(full_data))
  
  max11 <- lapply(cc$Highest$ensemble, function(x) sum(x$hard$denominator)) %>% unlist() %>% max()
  max12 <- lapply(cc$Lowest$ensemble, function(x) sum(x$hard$denominator)) %>% unlist() %>% max()
  max1 <- max(c(max11, max12))
  max21 <- lapply(cc$Highest$ensemble, function(x) sum(x$hard$numerator)) %>% unlist() %>% max()
  max22 <- lapply(cc$Lowest$ensemble, function(x) sum(x$hard$numerator)) %>% unlist() %>% max()
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

# Visualizing results ----------------------------------------------------------
meta_data_full <- rbind(meta_data_barley, meta_data_lotus)
meta_data_full <- data.frame(meta_data_full[,-1],
                             row.names = meta_data_full$Sample_ID)
res2[,Genotype:=meta_data_full[res2$Sample_ID,"Genotype"]]

res2[Host == "Barley", Host:="Hordeum"]
res2[,":="(Host = factor(Host, levels = c("Lotus", "Hordeum")),
           Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
           Obs = factor(Obs, levels = c("NPK", "PK", "UF")),
           Pred = factor(Pred, levels = c("UF", "PK", "NPK")))]

res2[,Prediction := Obs == Pred]
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
  # guides(fill = guide_legend(nrow = 2))+
  facet_grid(Host ~ Compartment) +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -8),
        strip.background = element_rect(colour = NA),
        axis.title.y = element_text(size = 8, family = "Helvetica"),
        axis.title.x = element_text(size = 8, family = "Helvetica"),
        axis.text.y = element_text(size = 8, family = "Helvetica",
                                   colour = "black"),
        axis.text.x = element_text(size = 8, family = "Helvetica",
                                   colour = "black"),
        legend.text = element_text(size = 8, family = "Helvetica"),
        legend.title = element_text(size = 8, family = "Helvetica"),
        strip.text = element_text(size = 8, family = "Helvetica",
                                  face = "bold"),
        legend.key.size = unit(5, "mm"))+
  NULL -> g1; g1

Barplot_data[Plant == "Barley", Plant:="Hordeum"]
Barplot_data[,m:=paste(Plant, Compartment, sep = " \n")]
Barplot_data[,m:=factor(m, levels = c("Lotus \nRhizosphere",
                                      "Hordeum \nRhizosphere",
                                      "Lotus \nRoot",
                                      "Hordeum \nRoot"))]
SLR_dt[,m:=paste(Plant, Compartment)]
SLR_dt[Pred_type == "NPK", Pred_type:="NPK vs non-NPK"]
SLR_dt[Pred_type == "PK", Pred_type:="PK vs UF"]
# SLR_dt[,Pred_type:=paste(Pred_type, "Prediction")]

# Barplot_data$Order <- factor(Barplot_data$Order,
#                              levels = )
# Barplot_data <- Barplot_data[order(Order)]
reorder <- c("Burkholderiales", "Caulobacterales",
             "Flavobacteriales", "Pseudomonadales",
             "Rhizobiales", "Streptomycetales",
             "Micrococcales", "Unknown")
order_colors <- order_colors[match(reorder, order_colors$group),]
ggplot(data = Barplot_data, aes(x = Soil, y = Abundance, fill = Order))+
  geom_bar(stat = "identity", position = "stack", linewidth = 0.1) +
  # facet_wrap(~Plant+Compartment, nrow = 1)+
  # facet_wrap(~m, nrow = 1)+
  facet_wrap2(vars(Plant, Compartment), strip = strip_nested(), nrow = 1)+
  scale_fill_manual(values = order_colors$colors, breaks = order_colors$group,
                    name = "Bacterial orders")+
  # ylim(0, 1.05)+
  labs(x = NULL, y = "Cumulative Mean Relative Abundance")+
  # guides(fill = guide_legend(nrow = 3, title.position = "top"))+
  guides(fill = guide_legend(nrow = 3, title.position = "top"))+
  scale_y_continuous(expand = c(.001, .001), limits = c(0, 0.44))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.box.margin = margin(0, 0, 0, -15),
        legend.margin = margin(t = -5),
        strip.background = element_rect(colour = NA),
        legend.key.size = unit(0.25, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'),
        # panel.background = element_rect(fill = "white"),
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

# ggplot(data = SLR_dt, aes(x = Soil, y = SLR, fill = Soil))+
#   geom_boxplot()+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   facet_grid(Pred_type~m)+
#   guides(fill = "none")+
#   scale_fill_manual(values = colors, breaks = names(colors)) +
#   labs(x = NULL, y = "Most informative log-ratio")+
#   NULL -> g3

ggplot() + theme_void() -> blank

R <- Res[,c(1,2,6,5,3)]
R[,Accuracy:=paste0(Accuracy, "%")]
R[Host == "Barley", Host:="Hordeum"]
colnames(R) <- gsub("_", " ", colnames(R))
R <- R[c(3:4,1:2)]
# setcolorder(R, c("Host", "Compartment", "PK ratios", "NPK ratios", "Accuracy"))

# gg1 <- ggarrange(g1, g2, labels = c("A", "B"))
# gg2 <- ggarrange(g3, blank, widths = c(0.6, 0.4), labels = c("C", "D"))
# gg1 <- ggarrange(blank, g1, labels = c("A", "B"), widths = c(0.4, 0.6))
# gg2 <- ggarrange(g2, tableGrob(R), labels = c("C", "D"))
# gg <- ggarrange(gg1, gg2, nrow = 2)

# source("Flowchart.R")
# fc <- readRDS("Flowchart.rds")

# IMPORTANT! 
# Run Flowchart.R to produce the correct figure!

tg <- tableGrob(R, theme = ttheme_default(base_size = 8), rows = NULL)

gg1 <- plot_grid(grob, tg,
                 labels = c("A", "C"),
                 rel_heights = c(0.8, 0.2),
                 ncol = 1)
gg2 <- ggarrange(g1, g2, labels = c("B", "D"), ncol = 1, heights = c(0.5, 0.5))
gg <- ggarrange(gg1, gg2, ncol = 2)

ggsave(filename = "Prediction_plot_with_flowchart.pdf", plot = gg,
       width = 210, height = 200, units = "mm")

fwrite(Res, "Prediction_summary.csv")
