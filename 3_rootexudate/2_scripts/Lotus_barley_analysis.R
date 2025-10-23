# Code by Ib Thorsgaard Jensen
# Loading pacakges -------------------------------------------------------------
pkg <- c("data.table", "magrittr", "ggplot2", "ggfortify", "multcompView",
         "ggpubr","RColorBrewer", "ComplexHeatmap", "colorRamp2", "gridtext",
         "cowplot")
for(pk in pkg){
  library(pk, character.only = T)
}

#Changes working directory to the path of this file. Only works in Rstudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
source("Tobit_metabolites.R")

cols <- c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD", "#cecece")
gt <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
names(cols) <- c(gt, "control")

# Loading data -----------------------------------------------------------------
## Lotus ----
metabolite_data_Lj <- fread(
  "../1_data/1_Lotus/241129_CSSP_mutants_neg_GNPS_quant.csv",
  header = T, drop = 2:13
)
metabolite_data_Lj[,V73:=NULL]
sample_numbers <- unlist(
  lapply(strsplit(colnames(metabolite_data_Lj)[-1], "_"), function(x) x[1])
)
colnames(metabolite_data_Lj)[-1] <- sample_numbers

meta_data_Lj <- fread("../1_data/1_Lotus/LotusCSSP_rootex_metadata.txt")
annotation_Lj <- fread(
  "../1_data/1_Lotus/Lotus_Dec24_canopus_structure_summary.tsv"
)

### Changing data format ----
#### metabolite data ----
colnames(metabolite_data_Lj)[1] <- "Feature"
metabolite_data_Lj[,Feature:= paste0("Feature", Feature)]

QC_idx <- grepl("QC", colnames(metabolite_data_Lj))
col_rem <- colnames(metabolite_data_Lj)[QC_idx]
metabolite_data_Lj[,c(col_rem):=NULL]
colnames(metabolite_data_Lj)[-1] <- paste0(
  "Sample", colnames(metabolite_data_Lj)[-1]
)

#### metadata ----
meta_data_Lj[,Sample_ID:=paste0("Sample", Sample_ID)]
meta_data_Lj <- meta_data_Lj[ Sample_ID %in% colnames(metabolite_data_Lj)]

#### Annotations ----
w <- which(colnames(annotation_Lj) == "mappingFeatureId")
colnames(annotation_Lj)[w] <- "Feature"
setcolorder(annotation_Lj, "Feature")
annotation_Lj[,Feature:= paste0("Feature", Feature)]

## Barley ----------------------------------------------------------------------
metabolite_data_Hv <- fread("../1_data/2_Barley/01042025_Barley_neg_quant.csv",
                            header = T, drop = 2:13)
metabolite_data_Hv[,V38:=NULL]
sample_numbers <- unlist(
  lapply(strsplit(colnames(metabolite_data_Hv)[-1], "_"), function(x) x[4])
)
colnames(metabolite_data_Hv)[-1] <- sample_numbers

meta_data_Hv <- fread("../1_data/2_Barley/BarleyCSSP_3diaz_metadata.txt")
annotation_Hv <- fread(
  "../1_data/2_Barley/Barley_Apr25_canopus_structure_summary.tsv"
)

meta_data_Hv <- meta_data_Hv[Sample_ID %in% sample_numbers]
meta_data_Hv[Genotype == "sand", Genotype := "control"]

### Changing data format ----
#### Metabolite data ----
colnames(metabolite_data_Hv)[1] <- "Feature"
metabolite_data_Hv[,Feature:= paste0("Feature", Feature)]
colnames(metabolite_data_Hv)[-1] <- paste0(
  "Sample", colnames(metabolite_data_Hv)[-1]
)

#### Metadata ----
meta_data_Hv[,Sample_ID:=paste0("Sample", Sample_ID)]
meta_data_Hv <- meta_data_Hv[Sample_ID %in% colnames(metabolite_data_Hv)]

#### Annotations ----
w <- which(colnames(annotation_Hv) == "mappingFeatureId")
colnames(annotation_Hv)[w] <- "Feature"
setcolorder(annotation_Hv, "Feature")
annotation_Hv[,Feature:= paste0("Feature", Feature)]

# Subsetting -------------------------------------------------------------------

# Subset only Inoculum = none
sample_keep <- meta_data_Lj[Inoculum == "none", Sample_ID]
meta_data_Lj <- meta_data_Lj[Sample_ID %in% sample_keep]
metabolite_data_Lj <- metabolite_data_Lj[,c("Feature", sample_keep), with = F]

sample_keep <- meta_data_Hv[Inoculum == "none", Sample_ID]
meta_data_Hv <- meta_data_Hv[Sample_ID %in% sample_keep]
metabolite_data_Hv <- metabolite_data_Hv[,c("Feature", sample_keep), with = F]

# Remove features present in too few replicates
is_feature_present <- rowMeans(metabolite_data_Lj[,-1] == 0) <= 0.9
feature_present <- metabolite_data_Lj$Feature[is_feature_present]
metabolite_data_Lj <- metabolite_data_Lj[Feature %in% feature_present]
metabolite_data_Lj <- metabolite_data_Lj[Feature %in% annotation_Lj$Feature]

is_feature_present <- rowMeans(meta_data_Hv[,-1] == 0) <= 0.9
feature_present <- metabolite_data_Hv$Feature[is_feature_present]
metabolite_data_Hv <- metabolite_data_Hv[Feature %in% feature_present]
metabolite_data_Hv <- metabolite_data_Hv[Feature %in% annotation_Hv$Feature]

# PCA plots --------------------------------------------------------------------
feattab_trans_Lj <- transpose(l = metabolite_data_Lj,
                              keep.names = "Sample_ID",
                              make.names = "Feature")

feattab_meta_Lj <- merge(meta_data_Lj[,c("Sample_ID", "Genotype")],
                         feattab_trans_Lj,
                         by = "Sample_ID")

feattab_trans_Hv <- transpose(l = metabolite_data_Hv,
                              keep.names = "Sample_ID",
                              make.names = "Feature")

feattab_meta_Hv <- merge(meta_data_Hv[,c("Sample_ID", "Genotype")],
                         feattab_trans_Hv,
                         by = "Sample_ID")

pLj <- prcomp(feattab_meta_Lj[,-(1:2)], center = T, scale. = T)
var_exp_Lj <- (pLj$sdev^2/sum(pLj$sdev^2))[1:2]*100
var_exp_Lj <- paste0("(", round(var_exp_Lj, 2), "%)")

pHv <- prcomp(feattab_meta_Hv[,-(1:2)], center = T, scale. = T)
var_exp_Hv <- (pHv$sdev^2/sum(pHv$sdev^2))[1:2]*100
var_exp_Hv <- paste0("(", round(var_exp_Hv, 2), "%)")

dt1 <- data.table(pLj$x[,1:2], Genotype = feattab_meta_Lj$Genotype,
                  Host = "Lotus")
dt2 <- data.table(pHv$x[,1:2], Genotype = feattab_meta_Hv$Genotype,
                  Host = "Hordeum")
dt <- rbind(dt1, dt2)
ggplot(data = dt1, mapping = aes(x = PC1, y = PC2, fill = Genotype))+
  geom_point(size = 4, shape = 21)+
  theme_bw()+
  scale_fill_manual(name = "Genotype",
                    breaks = names(cols),
                    values = cols)+
  labs(x = paste("PC1", var_exp_Lj[1]), y = paste("PC2", var_exp_Lj[2]))+
  facet_wrap(~Host)+
  theme(legend.position = "bottom",
        strip.background = element_rect(colour = NA),
        strip.text = element_text(face = "bold"))+
  guides(fill = guide_legend(nrow = 1))+
  NULL -> p1

ggplot(data = dt2, mapping = aes(x = PC1, y = PC2, fill = Genotype))+
  geom_point(size = 4, shape = 21)+
  theme_bw()+
  scale_fill_manual(name = "Genotype",
                    breaks = names(cols),
                    values = cols)+
  labs(x = paste("PC1", var_exp_Hv[1]),
       y = paste("PC2", var_exp_Hv[2]))+
  facet_wrap(~Host)+
  theme(legend.position = "bottom",
        strip.background = element_rect(colour = NA),
        strip.text = element_text(face = "bold"))+
  guides(fill = guide_legend(nrow = 1))+
  NULL -> p2

ggarrange(p1, p2, common.legend = TRUE, legend = "bottom")

# Filter based on sand ---------------------------------------------------------
## Lotus ----
metabolites_Lj <- data.frame(metabolite_data_Lj[,-1],
                             row.names = metabolite_data_Lj$Feature)
meta_data_Lj_df <- data.frame(meta_data_Lj[,-1], 
                              row.names = meta_data_Lj$Sample_ID)
meta_data_Lj_df$Genotype <- factor(meta_data_Lj_df$Genotype,
                                   levels = c("control", "WT", "symrk",
                                              "ccamk", "nsp1", "nsp2"))

sand_filter <- Tobit_model(metabolites_Lj, meta_data_Lj_df,
                           formula = ~Genotype)
test_WT <- Tobit_LRT(sand_filter, var = "GenotypeWT", p_adjust = "fdr")
test_symrk <- Tobit_LRT(sand_filter, var = "Genotypesymrk", p_adjust = "fdr")
test_ccamk <- Tobit_LRT(sand_filter, var = "Genotypeccamk", p_adjust = "fdr")
test_nsp1 <- Tobit_LRT(sand_filter, var = "Genotypensp1", p_adjust = "fdr")
test_nsp2 <- Tobit_LRT(sand_filter, var = "Genotypensp2", p_adjust = "fdr")

p_vals_mat_Lj <- matrix(NA, nrow = nrow(metabolites_Lj), ncol = 5)
colnames(p_vals_mat_Lj) <- names(cols)[-6]
p_vals_mat_Lj[,"WT"] <- test_WT$p_adj
p_vals_mat_Lj[,"symrk"] <- test_symrk$p_adj
p_vals_mat_Lj[,"ccamk"] <- test_ccamk$p_adj
p_vals_mat_Lj[,"nsp1"] <- test_nsp1$p_adj
p_vals_mat_Lj[,"nsp2"] <- test_nsp2$p_adj

feature_keep_Lj <- rownames(test_nsp2)[apply(p_vals_mat_Lj<0.05, 1, any)]

## Barley ----
metabolites_Hv <- data.frame(metabolite_data_Hv[,-1],
                             row.names = metabolite_data_Hv$Feature)
meta_data_Hv_df <- data.frame(meta_data_Hv[,-1], 
                              row.names = meta_data_Hv$Sample_ID)
meta_data_Hv_df$Genotype <- factor(meta_data_Hv_df$Genotype,
                                   levels = c("control", "WT", "symrk",
                                              "ccamk", "nsp1", "nsp2"))

sand_filter <- Tobit_model(metabolites_Hv, meta_data_Hv_df,
                           formula = ~Genotype)
test_WT <- Tobit_LRT(sand_filter, var = "GenotypeWT", p_adjust = "fdr")
test_symrk <- Tobit_LRT(sand_filter, var = "Genotypesymrk", p_adjust = "fdr")
test_ccamk <- Tobit_LRT(sand_filter, var = "Genotypeccamk", p_adjust = "fdr")
test_nsp1 <- Tobit_LRT(sand_filter, var = "Genotypensp1", p_adjust = "fdr")
test_nsp2 <- Tobit_LRT(sand_filter, var = "Genotypensp2", p_adjust = "fdr")

p_vals_mat_Hv <- matrix(NA, nrow = nrow(metabolites_Hv), ncol = 5)
colnames(p_vals_mat_Hv) <- names(cols)[-6]
p_vals_mat_Hv[,"WT"] <- test_WT$p_adj
p_vals_mat_Hv[,"symrk"] <- test_symrk$p_adj
p_vals_mat_Hv[,"ccamk"] <- test_ccamk$p_adj
p_vals_mat_Hv[,"nsp1"] <- test_nsp1$p_adj
p_vals_mat_Hv[,"nsp2"] <- test_nsp2$p_adj

feature_keep_Hv <- rownames(test_nsp2)[apply(p_vals_mat_Hv<0.05, 1, any)]

# Removing non-present features
# Lotus
samples_keep <- meta_data_Lj[Genotype != "control", Sample_ID]
meta_data_Lj_df <- meta_data_Lj_df[samples_keep,]
metabolites_Lj <- metabolites_Lj[feature_keep_Lj, samples_keep]

feature_present_Lj <- rownames(metabolites_Lj)[
  apply(metabolites_Lj, 1, function(x) mean(x!=0)>=0.1)
]
metabolites_Lj <- metabolites_Lj[feature_present_Lj,]

# Barley
samples_keep <- meta_data_Hv[Genotype != "control", Sample_ID]
meta_data_Hv_df <- meta_data_Hv_df[samples_keep,]
metabolites_Hv <- metabolites_Hv[feature_keep_Hv, samples_keep]

feature_present_Hv <- rownames(metabolites_Hv)[
  apply(metabolites_Hv, 1, function(x) mean(x!=0)>=0.1)
]
metabolites_Hv <- metabolites_Hv[feature_present_Hv,]

# PCA plots after filtering ----------------------------------------------------
pLj <- prcomp(t(metabolites_Lj), center = T, scale. = T)
pHv <- prcomp(t(metabolites_Hv), center = T, scale. = T)

var_exp_Lj <- (pLj$sdev^2/sum(pLj$sdev^2))[1:2]*100
var_exp_Lj <- paste0("(", round(var_exp_Lj, 2), "%)")

var_exp_Hv <- (pHv$sdev^2/sum(pHv$sdev^2))[1:2]*100
var_exp_Hv <- paste0("(", round(var_exp_Hv, 2), "%)")

dtLj <- data.table(pLj$x[,1:2], Genotype = meta_data_Lj_df$Genotype,
                   Host = "Lotus")
dtHv <- data.table(pHv$x[,1:2], Genotype = meta_data_Hv_df$Genotype,
                   Host = "Hordeum")

centroids_Lj <- dtLj[,.(PC1_cent = mean(PC1), PC2_cent = mean(PC2)),
                     list(Genotype, Host)]
segments_Lj <- merge(dtLj, centroids_Lj, by=c("Genotype", "Host"))
ggplot(data = dtLj, mapping = aes(x = PC1, y = PC2, colour = Genotype))+
  geom_point(size = 2, stroke = 0.25)+
  geom_segment(data = segments_Lj, aes(x = PC1, y = PC2, xend = PC1_cent,
                                       yend = PC2_cent, color = Genotype),
               alpha = 0.5, show.legend = FALSE)+
  scale_colour_manual(name = "Genotype",
                    breaks = names(cols),
                    values = cols)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(colour = 'black', size = 8, face = "bold"),
        strip.background = element_rect(colour = NA),
        # legend.position = "none",
        legend.position = "bottom",
        # legend.background = element_rect(linewidth = 0.2, colour = 1),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8,
                                   margin = margin(l = -0.1, unit = "pt")),
        legend.title = element_text(size = 8),
        legend.margin = margin(t = 0, r = 5, l = 5),
        legend.key = element_rect(fill = NA),
        legend.key.spacing.y = unit(-0.1, "cm"),
        legend.key.spacing.x = unit(0.001, "cm"),
        plot.title = element_text(size = 8, hjust = 0.5))+
  labs(x = paste("PC1", var_exp_Lj[1]), y = paste("PC2", var_exp_Lj[2]))+
  facet_wrap(~Host)+
  guides(fill = guide_legend(nrow = 1))+
  NULL -> p1

centroids_Hv <- dtHv[,.(PC1_cent = mean(PC1), PC2_cent = mean(PC2)),
                     list(Genotype, Host)]
segments_Hv <- merge(dtHv, centroids_Hv, by=c("Genotype", "Host"))
ggplot(data = dtHv, mapping = aes(x = PC1, y = PC2, colour = Genotype))+
  geom_point(size = 2, stroke = 0.25)+
  scale_colour_manual(name = "Genotype",
                    breaks = names(cols),
                    values = cols)+
  geom_segment(data = segments_Hv, aes(x = PC1, y = PC2, xend = PC1_cent,
                                       yend = PC2_cent, color = Genotype),
               alpha = 0.5, show.legend = FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(colour = 'black', size = 8, face = "bold"),
        strip.background = element_rect(colour = NA),
        legend.position = "none",
        # legend.position = "bottom",
        # legend.background = element_rect(size = 0.2, colour = 1),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8,
                                   margin = margin(l = -0.1, unit = "pt")),
        legend.title = element_text(size = 8),
        legend.margin = margin(t = 0, r = 5, l = 5),
        legend.key = element_rect(fill = NA),
        legend.key.spacing.y = unit(-0.1, "cm"),
        legend.key.spacing.x = unit(0.001, "cm"),
        plot.title = element_text(size = 8, hjust = 0.5))+
  labs(x = paste("PC1", var_exp_Hv[1]),
       y = paste("PC2", var_exp_Hv[2]))+
  facet_wrap(~Host)+
  guides(fill = guide_legend(nrow = 1))+
  NULL -> p2

PCA_legend <- get_plot_component(p1, 'guide-box', return_all = TRUE)
PCA_legend <- PCA_legend[[3]]
p1 <- p1 + guides(fill = "none")

PCA_all <- ggarrange(p1, p2, nrow = 2)
PCA_all + guides(colour = NULL)
# Genotyope effects ------------------------------------------------------------
## Lotus ----
meta_data_Lj_df$Genotype <- droplevels(meta_data_Lj_df$Genotype)
DA_metabolites_Lj <- Tobit_model(metabolites_Lj,
                                 meta_data_Lj_df,
                                 formula = ~Genotype)
test_symrk_Lj <- Tobit_LRT(DA_metabolites_Lj,
                           var = "Genotypesymrk",
                           p_adjust = "fdr")
test_ccamk_Lj <- Tobit_LRT(DA_metabolites_Lj,
                           var = "Genotypeccamk",
                           p_adjust = "fdr")
test_nsp1_Lj <- Tobit_LRT(DA_metabolites_Lj,
                          var = "Genotypensp1",
                          p_adjust = "fdr")
test_nsp2_Lj <- Tobit_LRT(DA_metabolites_Lj,
                          var = "Genotypensp2",
                          p_adjust = "fdr")
p_adj_mat_Lj <- matrix(NA, nrow = nrow(metabolites_Lj), ncol = 4)
colnames(p_adj_mat_Lj) <- gt[-1]
p_adj_mat_Lj[,"symrk"] <- test_symrk_Lj$p_adj
p_adj_mat_Lj[,"ccamk"] <- test_ccamk_Lj$p_adj
p_adj_mat_Lj[,"nsp1"] <- test_nsp1_Lj$p_adj
p_adj_mat_Lj[,"nsp2"] <- test_nsp2_Lj$p_adj

### Savings results ----
DA_Lj <- matrix(NA, nrow = nrow(metabolites_Lj), ncol = 4)
colnames(DA_Lj) <- gt[-1]
rownames(DA_Lj) <- rownames(metabolites_Lj)
DA_Lj[,"symrk"] <- test_symrk_Lj$p_adj < 0.05
DA_Lj[,"ccamk"] <- test_ccamk_Lj$p_adj < 0.05
DA_Lj[,"nsp1"] <- test_nsp1_Lj$p_adj < 0.05
DA_Lj[,"nsp2"] <- test_nsp2_Lj$p_adj < 0.05
feature_sig <- rownames(DA_Lj)[rowSums(DA_Lj)>0]
DA_Lj_full <- DA_Lj
DA_Lj <- DA_Lj[feature_sig,]

DA_metabolites_Lj$res
res_Lj <- ifelse(DA_metabolites_Lj$res[,2:5] > 0, 1, -1)
res_Lj <- res_Lj[rownames(DA_Lj),]
res_Lj[!DA_Lj] <- 0
res_Lj <- data.table(Feature = rownames(res_Lj), res_Lj)
colnames(res_Lj) <- gsub("Genotype", "", colnames(res_Lj))

Feat_tab_Lj <- metabolites_Lj
colnames(Feat_tab_Lj) <- paste(meta_data_Lj_df$Genotype,
                               rownames(meta_data_Lj_df), sep = "_")
Feat_tab_Lj <- data.table(Feature = rownames(Feat_tab_Lj), Feat_tab_Lj)

res_Lj <- merge(res_Lj, Feat_tab_Lj, by = "Feature")
res_Lj <- merge(res_Lj, annotation_Lj, by = "Feature", all.x = T)

fwrite(res_Lj, "../4_tables/Metabolite_test_results_tobit_Lj.csv")

## Barley ----
meta_data_Hv_df$Genotype <- droplevels(meta_data_Hv_df$Genotype)
DA_metabolites_Hv <- Tobit_model(metabolites_Hv,
                                 meta_data_Hv_df,
                                 formula = ~Genotype)
test_symrk_Hv <- Tobit_LRT(DA_metabolites_Hv,
                           var = "Genotypesymrk",
                           p_adjust = "fdr")
test_ccamk_Hv <- Tobit_LRT(DA_metabolites_Hv,
                           var = "Genotypeccamk",
                           p_adjust = "fdr")
test_nsp1_Hv <- Tobit_LRT(DA_metabolites_Hv,
                          var = "Genotypensp1",
                          p_adjust = "fdr")
test_nsp2_Hv <- Tobit_LRT(DA_metabolites_Hv, 
                          var = "Genotypensp2",
                          p_adjust = "fdr")
p_adj_mat_Hv <- matrix(NA, nrow = nrow(metabolites_Hv), ncol = 4)
colnames(p_adj_mat_Hv) <- gt[-1]
p_adj_mat_Hv[,"symrk"] <- test_symrk_Hv$p_adj
p_adj_mat_Hv[,"ccamk"] <- test_ccamk_Hv$p_adj
p_adj_mat_Hv[,"nsp1"] <- test_nsp1_Hv$p_adj
p_adj_mat_Hv[,"nsp2"] <- test_nsp2_Hv$p_adj

### Savings results ----
DA_Hv <- matrix(NA, nrow = nrow(metabolites_Hv), ncol = 4)
colnames(DA_Hv) <- gt[-1]
rownames(DA_Hv) <- rownames(metabolites_Hv)
DA_Hv[,"symrk"] <- test_symrk_Hv$p_adj < 0.05
DA_Hv[,"ccamk"] <- test_ccamk_Hv$p_adj < 0.05
DA_Hv[,"nsp1"] <- test_nsp1_Hv$p_adj < 0.05
DA_Hv[,"nsp2"] <- test_nsp2_Hv$p_adj < 0.05
feature_sig <- rownames(DA_Hv)[rowSums(DA_Hv)>0]
DA_Hv_full <- DA_Hv
DA_Hv <- DA_Hv[feature_sig,]

DA_metabolites_Hv$res
res_Hv <- ifelse(DA_metabolites_Hv$res[,2:5] > 0, 1, -1)
res_Hv <- res_Hv[rownames(DA_Hv),]
res_Hv[!DA_Hv] <- 0
res_Hv <- data.table(Feature = rownames(res_Hv), res_Hv)
colnames(res_Hv) <- gsub("Genotype", "", colnames(res_Hv))

Feat_tab_Hv <- metabolites_Hv
colnames(Feat_tab_Hv) <- paste(meta_data_Hv_df$Genotype, 
                               rownames(meta_data_Hv_df), sep = "_")
Feat_tab_Hv <- data.table(Feature = rownames(Feat_tab_Hv), Feat_tab_Hv)

res_Hv <- merge(res_Hv, Feat_tab_Hv, by = "Feature")
res_Hv <- merge(res_Hv, annotation_Hv, by = "Feature", all.x = T)

fwrite(res_Hv, "../4_tables/Metabolite_test_results_tobit_Hv.csv")

# Volcano plots ----------------------------------------------------------------
p_vals_dt <- rbind(
  data.table(Feature = rownames(DA_metabolites_Lj$res),
             p_adj_mat_Lj, Host = "Lotus"),
  data.table(Feature = rownames(DA_metabolites_Hv$res),
             p_adj_mat_Hv, Host = "Hordeum")
)

p_vals_dt <- melt(p_vals_dt,
                  id.vars = c("Feature", "Host"),
                  variable.name = "Genotype",
                  value.name = "p_adj")

logFC_dt <- rbind(
  data.table(Feature = rownames(DA_metabolites_Lj$res),
             DA_metabolites_Lj$res[,2:5], Host = "Lotus"),
  data.table(Feature = rownames(DA_metabolites_Hv$res),
             DA_metabolites_Hv$res[,2:5], Host = "Hordeum")
)

logFC_dt <- melt(logFC_dt,
                 id.vars = c("Feature", "Host"),
                 variable.name = "Genotype",
                 value.name = "logFC")

logFC_dt[,Genotype:=gsub("Genotype", "", Genotype)]
logFC_dt[,logFC:=logFC/log(2)]
res_table <- merge(x = logFC_dt, y = p_vals_dt,
                   by = c("Feature", "Genotype", "Host"))
res_table[,Sig:=p_adj<0.05]
res_table[,DA:=fcase(Sig == T & logFC > 0, "Enriched",
                     Sig == T & logFC < 0, "Depleted",
                     default = "NS")]
res_table[,Genotype:=factor(Genotype,
                            levels = c("symrk", "ccamk", "nsp1", "nsp2"))]

text_data <- res_table[,.(N_DA=sum(Sig)), list(Host, Genotype)]

res_table[,Host:=factor(Host, levels = c("Lotus", "Hordeum"))]
ggplot(data = res_table, aes(x = logFC, y = -log10(p_adj), colour = DA))+
  geom_point(size = 0.5) +
  xlab(expression("log"[2]*"FC vs WT")) +
  ylab(expression("-log"[10]*"p-value (adjuted for FDR)")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_light() +
  scale_color_manual(values = c("Enriched" = "#902121",
                                "Depleted" = "darkblue",
                                "NS" = "gray"),
                     name = "DEM") +
  facet_grid(factor(Host, levels = c("Lotus", "Hordeum"))~Genotype)+
  geom_label(data = text_data, aes(x = -25, y = 8.5, label = N_DA),
             colour = "black", fill = "grey", alpha = 0.2, size = 8/.pt)+
  theme(
    panel.border = element_rect(color = "black", size = 0.5),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, hjust = 0.5),
    strip.text = element_text(colour = 'black', size = 8, face = "bold"),
    # legend.position = "none",
    legend.position = "bottom",
    # legend.background = element_rect(size = 0.2, colour = 1),
    strip.background =element_rect(fill="lightgrey"),
    # legend.margin = margin(t = 2, r = 2, l = 2, b = 0),
    legend.key = element_rect(fill = NA),
    legend.key.spacing.x = unit(5, "pt"),
    legend.key.size = unit(5, "pt"),
    legend.box.spacing = unit(5, "pt"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, hjust = 0.5, margin = margin(r = 5))
  ) +
  xlim(-30,30)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  # guides(colour = "none")+
  NULL -> volcanoes

volcano_legend <- get_plot_component(volcanoes, 'guide-box', return_all = TRUE)
volcano_legend <- volcano_legend[[3]]
volcanoes <- volcanoes + guides(colour = "none")
plot_grid(volcano_legend)
# Boxplots ---------------------------------------------------------------------
## Lotus ----
feat_set <- paste0( "Feature", c(269, 1067, 455, 882, 973,
                                 1047, 945, 1053, 976) )

c("Feature269" = "Coumaric\nacid", "Feature1067" = "Chalcone",
  "Feature455" = "Ferulic acid", "Feature882" = "Liqiritigenin",
  "Feature973" = "Naringenin", "Feature1047" = "BiochaninA/\nOlmelin",
  "Feature945" = "Formononetin", "Feature1053" = "Vestitione",
  "Feature976" = "Vestitol") -> name_change

feat_set <- factor(feat_set, levels = feat_set)

feat_set_dt <- data.table(t(metabolites_Lj[as.character(feat_set),]),
                          Genotype = meta_data_Lj_df$Genotype)

feat_set_dt <- melt(feat_set_dt,
                     variable.name = "Feature",
                     value.name = "Intensity",
                     id.vars = length(feat_set)+1)

test_symrk_Lj[as.character(feat_set),]

r1 <- data.table(Feature = feat_set,
                 test_symrk_Lj[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "symrk")
r2 <- data.table(Feature = feat_set,
                 test_ccamk_Lj[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "ccamk")
r3 <- data.table(Feature = feat_set,
                 test_nsp1_Lj[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "nsp1")
r4 <- data.table(Feature = feat_set,
                 test_nsp2_Lj[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "nsp2")

rr <- rbind(r1, r2, r3, r4)
rr[,text:=ifelse(p_adj < 0.05, "*", NA)]

feat_set_dt2 <- feat_set_dt[,max(Intensity), list(Feature, Genotype)]
rr2 <- merge(rr, feat_set_dt2, by = c("Feature", "Genotype"))

# Manually adjusting position of aserisks
rr2[2, V1:=1000]
rr2[3, V1:=2000]
rr2[21, V1:=14000]
rr2[27, V1:=2000]
rr2[34, V1:=125000]

feat_set_dt[,Feature:=name_change[as.character(Feature)]]
feat_set_dt[,Feature:=factor(Feature, levels = name_change)]
rr2[,Feature:=name_change[as.character(Feature)]]
rr2[,Feature:=factor(Feature, levels = name_change)]

ggplot(feat_set_dt, aes(x=Genotype, y=Intensity, fill=Genotype)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5, outlier.color = "red") +
  facet_wrap(~Feature, scales = "free")+
  geom_label(data = rr2, aes(y = V1, label = text), label.size = NA,
             alpha = 0, size = 20/.pt)+
  scale_fill_manual(values = cols, breaks = names(cols))+
  theme_bw()+
  ggtitle("Lotus")+
  theme(legend.position = "right",
        strip.background = element_rect(colour = NA),
        axis.title.y=element_blank(),
        axis.title = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(colour = 'black', size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        legend.margin = margin(t = 20, unit = "pt"))+
  labs(x = NULL, title = "Lotus")+
  guides(fill = "none")+
  NULL -> feature_box_Lj

## Barley ----
feat_set <- paste0( "Feature", c(2546, 1548, 495, 2970, 3069,
                                 3095, 2306, 2308, 2309) )

c("Feature2546" = "Gibberellin", "Feature1548" = "Abscisic \nacid",
  "Feature495" = "Esculetin", "Feature2970" = "Flavonoid-like\nglycoside",
  "Feature3069" = "Isoorientin", "Feature3095" = "Paeonin C",
  "Feature2306" = "Oxylipin 1", "Feature2308" = "Oxylipin 2",
  "Feature2309" = "Oxylipin 3") -> name_change

feat_set <- factor(feat_set, levels = feat_set)

feat_set_dt <- data.table(t(metabolites_Hv[as.character(feat_set),]),
                          Genotype = meta_data_Hv_df$Genotype)

feat_set_dt <- melt(feat_set_dt,
                    variable.name = "Feature",
                    value.name = "Intensity",
                    id.vars = length(feat_set)+1)

test_symrk_Hv[as.character(feat_set),]

r1 <- data.table(Feature = feat_set,
                 test_symrk_Hv[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "symrk")
r2 <- data.table(Feature = feat_set,
                 test_ccamk_Hv[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "ccamk")
r3 <- data.table(Feature = feat_set,
                 test_nsp1_Hv[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "nsp1")
r4 <- data.table(Feature = feat_set,
                 test_nsp2_Hv[as.character(feat_set), c("p_vals", "p_adj")],
                 Genotype = "nsp2")

rr <- rbind(r1, r2, r3, r4)
rr[,text:=ifelse(p_adj < 0.05, "*", NA)]

feat_set_dt2 <- feat_set_dt[,max(Intensity), list(Feature, Genotype)]
rr2 <- merge(rr, feat_set_dt2, by = c("Feature", "Genotype"))
rr2[Feature == "Feature3329",V1:=V1-1000]
rr2[Feature == "Feature387",V1:=V1-40]
rr2[V1<0, V1:=0]

# Manually adjusting position of aserisks
rr2[3, V1:=100]
rr2[6, V1:=80]
rr2[7, V1:=60]
rr2[21, V1:=950]
rr2[27,V1:=3100]
rr2[31,V1:=7400]
rr2[35,V1:=6400]

data.table(Feature = paste0("Feature", c(2306, 2308, 2309)),
           Genotype = c("nsp1", "nsp2", "nsp2"),
           p_vals = rep(NA, 3),
           p_adj = rep(NA, 3),
           text = rep("", 3),
           V1 = c(4000, 12800, 8100)) -> extend_limits

rr2 <- rbind(rr2, extend_limits)

feat_set_dt[,Feature:=name_change[as.character(Feature)]]
feat_set_dt[,Feature:=factor(Feature, levels = name_change)]
rr2[,Feature:=name_change[as.character(Feature)]]
rr2[,Feature:=factor(Feature, levels = name_change)]

ggplot(feat_set_dt, aes(x=Genotype, y=Intensity, fill=Genotype)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5, outlier.color = "red") +
  facet_wrap(~Feature, scales = "free")+
  geom_label(data = rr2, aes(y = V1, label = text), label.size = NA,
             alpha = 0, size = 20/.pt)+
  scale_fill_manual(values = cols, breaks = names(cols))+
  theme_bw()+
  ggtitle("Hordeum")+
  theme(legend.position = "right",
        strip.background = element_rect(colour = NA),
        axis.title.y=element_blank(),
        axis.title = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(colour = 'black', size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        legend.margin = margin(t = 20, unit = "pt"))+
  labs(x = NULL)+
  guides(fill = "none")+
  NULL -> feature_box_Hv

feature_box <- ggarrange(feature_box_Lj, feature_box_Hv, ncol = 1)

# Bubble plot ------------------------------------------------------------------
## Sum tables ----
### Lotus ----
sum_tab_Lj <- data.table(
  res_Lj[,c("Feature", "symrk", "ccamk", "nsp1", "nsp2",
            "ClassyFire#class", "ClassyFire#class Probability")]
)
colnames(sum_tab_Lj)[6:7] <- c("Class", "Class_prob")
sum_tab_Lj[is.na(Class) | Class == "", Class:="Unknown class"]
sum_tab_Lj[Class_prob<0.6, Class:="Low probability"]
sum_tab_Lj[,":="(Feature = NULL, Class_prob = NULL)]
sum_tab_Lj_pos <- sum_tab_Lj[,lapply(.SD, function(x) sum(x == 1) ), Class]
sum_tab_Lj_neg <- sum_tab_Lj[,lapply(.SD, function(x) sum(x == -1) ), Class]

sum_tab_Lj_pos <- melt(sum_tab_Lj_pos,
                       id.vars = "Class",
                       variable.name = "Genotype",
                       value.name = "Amount_pos")
sum_tab_Lj_neg <- melt(sum_tab_Lj_neg,
                       id.vars = "Class",
                       variable.name = "Genotype",
                       value.name = "Amount_neg")
sum_tab_Lj <- merge(sum_tab_Lj_pos, sum_tab_Lj_neg,
                    by = c("Class", "Genotype"))
sum_tab_Lj[,Host:="Lotus"]

### Barley ----
sum_tab_Hv <- data.table(
  res_Hv[,c("Feature", "symrk", "ccamk", "nsp1", "nsp2",
            "ClassyFire#class", "ClassyFire#class Probability")]
)
colnames(sum_tab_Hv)[6:7] <- c("Class", "Class_prob")
sum_tab_Hv[is.na(Class) | Class == "", Class:="Unknown class"]
sum_tab_Hv[Class_prob<0.6, Class:="Low probability"]
sum_tab_Hv[,":="(Feature = NULL, Class_prob = NULL)]
sum_tab_Hv_pos <- sum_tab_Hv[,lapply(.SD, function(x) sum(x == 1) ), Class]
sum_tab_Hv_neg <- sum_tab_Hv[,lapply(.SD, function(x) sum(x == -1) ), Class]

sum_tab_Hv_pos <- melt(sum_tab_Hv_pos,
                       id.vars = "Class",
                       variable.name = "Genotype",
                       value.name = "Amount_pos")
sum_tab_Hv_neg <- melt(sum_tab_Hv_neg,
                       id.vars = "Class",
                       variable.name = "Genotype",
                       value.name = "Amount_neg")
sum_tab_Hv <- merge(sum_tab_Hv_pos, sum_tab_Hv_neg,
                    by = c("Class", "Genotype"))
sum_tab_Hv[,Host:="Hordeum"]

### Combining ----
sum_tab <- rbind(sum_tab_Lj, sum_tab_Hv)
sum_tab[,Amount:=Amount_pos+Amount_neg]

Total_DEM <- sum_tab[,.(N = sum(Amount)), list(Host, Genotype)]

AA <- sum_tab[,max(Amount), Class]
sp1 <- AA[V1 > 1, Class]
sum_tab <- sum_tab[Class %in% sp1]
sum_tab[,Host:=factor(Host, levels = c("Lotus", "Hordeum"))]
sum_tab <- sum_tab[Class != "Unknown class"]
sum_tab[,Direction:=fcase(
  Amount_pos/Amount_neg > 1.25, "Enriched",
  Amount_pos/Amount_neg < 0.8, "Depleted",
  default = "Similar"
)]
sum_tab <- sum_tab[Amount != 0]

# 'Cheating' by shortening class names. Remove if not desired.
# sum_tab[,Class:=gsub(" and derivatives", "", Class)]
# sum_tab[,Class:=gsub(" and substituted derivatives", "", Class)]
# sum_tab <- sum_tab[!(Class %in% "Allyl-type 1,3-dipolar organic compounds")]
# End of 'cheating'

avg_sum <- sum_tab[,.(avg = mean(Amount)),Class]
avg_sum <- avg_sum[order(avg)]
sum_tab[,Class:=factor(Class, levels = avg_sum$Class)]

dummy <- data.table(Class = "Total DEMs", Genotype = rep(gt[-1], 2),
                    Amount_pos = NA,
                    Amount_neg = NA,
                    Host = rep(c("Lotus", "Hordeum"), each = 2),
                    Amount = NA, 
                    Direction = NA)
sum_tab <- rbind(sum_tab, dummy)
sum_tab[,pos_prob:=Amount_pos/Amount]
Text <- data.table(Class = "Total DEMs", 
                   Host = rep(c("Lotus", "Hordeum"), each = 4 ),
                   Genotype = rep(gt[-1], 2))
Text <- merge(Text, Total_DEM, by = c("Host", "Genotype"))
# Text[,Class:=""]

sum_tab[,Host:=factor(Host, levels = c("Lotus", "Hordeum"))]
ggplot(data = sum_tab, mapping = aes(x = Genotype, y = Class, size = Amount,
                                     fill = pos_prob))+
  geom_label(data = Text, mapping = aes(x = Genotype, y = Class, label = N), 
             size = 3, fill = "white", label.size = NA)+
  geom_point(shape = 21)+
  facet_wrap(~factor(Host, levels = c("Lotus", "Hordeum")))+
  scale_fill_gradient2(midpoint = 0.5, low = "darkblue",
                       mid = "white", high = "#902121",
                       name = "Proportion of enriched DEMs")+
  scale_size_continuous(breaks = c(1, 10, 50, 100))+
  # scale_fill_manual(values = c("Enriched" = "#902121",
  #                              "Depleted" = "darkblue",
  #                              "Similar" = "white"),
  #                   name = "DEMs are mostly") +
  labs(size = "# of DEMs")+
  ggtitle("")+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        strip.background = element_rect(colour = NA),
        panel.border = element_rect(color = "black", size = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8,
                                   colour = "black"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(colour = 'black', size = 8, face = "bold"),
        legend.position = "bottom",
        legend.box="vertical",
        legend.justification = c(1.5, 0),
        # legend.background = element_rect(size = 0.2, colour = 1),
        legend.text = element_text(size = 8,
                                   margin = margin(t = 2, unit = "pt")),
        legend.title = element_text(size = 8,
                                    margin = margin(r = 5)),
        legend.key = element_rect(fill = NA),
        legend.key.spacing.x = unit(5, "pt"),
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.margin = margin(t = 2, r = 8, l = 5, unit = "pt"))+
  # guides(size = "none", fill = "none")+
  # guides(fill = guide_colourbar(barheight = 0.5, label.position = "bottom"))+
  NULL -> Bubble

#  Collect plots ---------------------------------------------------------------
Blank <- ggplot() + theme_void()

combined_legend <- ggarrange(PCA_legend, volcano_legend,
                             )
gg <- ggarrange(PCA_all, volcanoes, widths = c(0.2, 0.8))
gg2 <- plot_grid(gg, combined_legend, nrow = 2, rel_heights = c(0.95, 0.05))
gg3 <- plot_grid(Bubble, feature_box, rel_widths = c(0.54, 0.46))
gg4 <- plot_grid(gg2, gg3, nrow = 2, rel_heights = c(0.35, 0.65))

ggsave("../3_figures/Lotus_barley_metabolites.pdf", gg4,
       width = 210, height = 297, units = "mm")

