# Load packages and data -------------------------------------------------------
pkg <- c("data.table", "magrittr", "ggplot2", "vegan", "Maaslin2", "patchwork",
         "RColorBrewer", "ComplexHeatmap", "colorRamp2", "ggh4x")
for(pk in pkg){
  library(pk, character.only = T)
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Structural_zeros.R")

cols <- c("WT" = "#33a02c",
          "symrk" = "#ff7f00",
          "ccamk" = "#1f78b4",
          "nsp1" = "#e31a1c",
          "nsp2" = "#ffd700")

order_colors <- data.frame(group=c("Burkholderiales","Caulobacterales",
                                   "Chloroflexales", "Flavobacteriales",
                                   "Frankiales","Gaiellales","Gemmatimonadales",
                                   "Micrococcales","Micromonosporales",
                                   "Propionibacteriales", "Pseudomonadales",
                                   "Pseudonocardiales","Rhizobiales",
                                   "Sphingomonadales","Streptomycetales",
                                   "Unknown", "Xanthomonadales","Other"),
                           colors=c("#645394","#AA4488","#CC99BB","#ffeeef",
                                    "#114477","#4477AA","#77AADD","#44AAAA",
                                    "#77CCCC","#117744","#88CCAA","#CDEBC5",
                                    "lightyellow","#fdbb6b","#ffd7b5",
                                    "darkgrey","#ffc0cb","lightgrey"))

ASV_table <- fread(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
colnames(ASV_table)[1] <- "ASV_ID"
ASV_table_full <- data.table(ASV_table)
meta_data <- fread(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", drop = c(5,7,8)
)
meta_data_full <- data.table(meta_data)
taxonomy <- fread("../../1_data/1_Lotus/taxonomy_clean.tsv")

ASV_table_with_metadata <- transpose(ASV_table,
                                     make.names = "ASV_ID",
                                     keep.names = "SampleID")
ASV_table_with_metadata <- merge(meta_data,
                                 ASV_table_with_metadata,
                                 by = "SampleID")

# Setting reference levels for the soil and genotype
ASV_table_with_metadata[,Genotype:=factor(Genotype,
                                          levels = c("WT", "symrk", "ccamk",
                                                     "nsp1", "nsp2"))]

Experimental_conditions <- as.matrix(
  expand.grid(Compartment = c("Root", "Rhizosphere"),
              Soil = c("NPK", "PK", "UF"))
)

# Differential abundance analysis with MaAslin2 --------------------------------
ASV_tables_subsets <- apply(
  Experimental_conditions, 1,
  function(x){
    data_subset <- ASV_table_with_metadata[Compartment == x[1] & Soil == x[2]]
    
    meta_data <- data.frame(data_subset[,2:5],
                            row.names = data_subset$SampleID)
    meta_data$Genotype <- droplevels(meta_data$Genotype)
    
    ASV_table <- t(data_subset[,-(1:5)])
    colnames(ASV_table) <- data_subset$SampleID
    
    Out <- list(ASV_table = ASV_table, meta_data = meta_data)
    return(Out)
  }
)

lapply(ASV_tables_subsets,
       function(X){
         # Run MaAslin2
         M <- Maaslin2(input_data = X$ASV_table,
                       input_metadata = X$meta_data,
                       output = "Maaslin2 files",
                       plot_heatmap = F,
                       fixed_effects = "Genotype",
                       min_prevalence = 0.1)
         res_dt <- data.table(M$results)
         res_dt <- res_dt[order(feature)]
         
         # Run structural zeros function
         X$meta_data$library_size <- colSums(X$ASV_table)
         
         S <- Structural_zeros3(X$ASV_table, X$meta_data,
                                group = "Genotype", ref = "WT",
                                min_reads = 20, min_present_reps = 3)
         
         # Constructing logical matrix to indicate which ASVs are DA
         res_mat <- matrix(0, nrow = nrow(X$ASV_table), ncol = 4)
         rownames(res_mat) <- sort(rownames(X$ASV_table))
         colnames(res_mat) <- c("symrk", "ccamk", "nsp1", "nsp2")
         
         DA_symrk <- res_dt[value == "symrk" & qval < 0.05]$feature
         DA_symrk <- gsub("X", "", DA_symrk)
         res_mat[DA_symrk,"symrk"] <- sign(res_dt[value == "symrk" &
                                                    qval < 0.05]$coef)
         DA_ccamk <- res_dt[value == "ccamk" & qval < 0.05]$feature
         DA_ccamk <- gsub("X", "", DA_ccamk)
         res_mat[DA_ccamk,"ccamk"] <- sign(res_dt[value == "ccamk" &
                                                    qval < 0.05]$coef)
         DA_nsp1 <- res_dt[value == "nsp1" & qval < 0.05]$feature
         DA_nsp1 <- gsub("X", "", DA_nsp1)
         res_mat[DA_nsp1,"nsp1"] <- sign(res_dt[value == "nsp1" &
                                                  qval < 0.05]$coef)
         DA_nsp2 <- res_dt[value == "nsp2" & qval < 0.05]$feature
         DA_nsp2 <- gsub("X", "", DA_nsp2)
         res_mat[DA_nsp2,"nsp2"] <- sign(res_dt[value == "nsp2" &
                                                  qval < 0.05]$coef)
         
         all_struc_zero_DA_ASVS <- unlist(S$struc_zero_DA)
         N_DA_ASV_by_genotype <- unlist(lapply(S$struc_zero_DA, length))
         names(all_struc_zero_DA_ASVS) <- rep(c("symrk", "ccamk",
                                                "nsp1", "nsp2"),
                                              N_DA_ASV_by_genotype)
                                              
         data.table(ASV = all_struc_zero_DA_ASVS,
                    Genotype = names(all_struc_zero_DA_ASVS),
                    DA_by_struc_zero = T) -> struc_zero_res_dt
         
         # Fill in structural zeros results from symrk
         struc_zero_symrk <- setdiff(S$struc_zero_DA$symrk,
                                     res_dt[value == "symrk" &
                                              qval < 0.05]$feature)
         S_symrk <- S$struc_zero_table[struc_zero_symrk, c("WT", "symrk")]
         struc_zero_FC_symrk <- apply(S_symrk, 1, diff)
         res_mat[names(struc_zero_FC_symrk), "symrk"] <- struc_zero_FC_symrk
         
         # Fill in structural zeros results from ccamk
         struc_zero_ccamk <- setdiff(S$struc_zero_DA$ccamk,
                                     res_dt[value == "ccamk" &
                                              qval < 0.05]$feature)
         S_ccamk <- S$struc_zero_table[struc_zero_ccamk, c("WT", "ccamk")]
         struc_zero_FC_ccamk <- apply(S_ccamk, 1, diff)
         res_mat[names(struc_zero_FC_ccamk), "ccamk"] <- struc_zero_FC_ccamk
         
         # Fill in structural zeros results from nsp1
         struc_zero_nsp1 <- setdiff(S$struc_zero_DA$nsp1,
                                     res_dt[value == "nsp1" &
                                              qval < 0.05]$feature)
         S_nsp1 <- S$struc_zero_table[struc_zero_nsp1, c("WT", "nsp1")]
         struc_zero_FC_nsp1 <- apply(S_nsp1, 1, diff)
         res_mat[names(struc_zero_FC_nsp1), "nsp1"] <- struc_zero_FC_nsp1
         
         # Fill in structural zeros results from nsp2
         struc_zero_nsp2 <- setdiff(S$struc_zero_DA$nsp2,
                                     res_dt[value == "nsp2" &
                                              qval < 0.05]$feature)
         S_nsp2 <- S$struc_zero_table[struc_zero_nsp2, c("WT", "nsp2")]
         struc_zero_FC_nsp2 <- apply(S_nsp2, 1, diff)
         res_mat[names(struc_zero_FC_nsp2), "nsp2"] <- struc_zero_FC_nsp2
         
         # Mean RAs across all samples of DA ASVs
         WT_samples <- rownames(X$meta_data)[X$meta_data$Genotype == "WT"]
         lib_size <- colSums(X$ASV_table[,WT_samples])
         apply(res_mat, 2,
               function(y){
                 ASVs <- rownames(res_mat)[y != 0]
                 
                 ASV_subset <- X$ASV_table[ASVs,WT_samples]
                 mean(rowSums(t(ASV_subset)/lib_size))
               }) -> mean_RAs
         
         mean_RA_res <- data.table(Compartment = X$meta_data$Compartment[1],
                                   Soil = X$meta_data$Soil[1],
                                   Genotype = names(mean_RAs),
                                   mean_RAs)
         
         # Summarising results
         res_dt[,":="(DA_by_sig = qval < 0.05, metadata = NULL, stderr = NULL,
                      pval = NULL, N = NULL, N.not.zero = NULL, name = NULL,
                      Compartment = X$meta_data$Compartment[1],
                      Soil = X$meta_data$Soil[1])]
         colnames(res_dt)[1:3] <- c("ASV", "Genotype", "LogFC")
         res_dt <- merge(res_dt, struc_zero_res_dt,
                         by = c("ASV", "Genotype"), all.x = T)
         res_dt[is.na(DA_by_struc_zero), DA_by_struc_zero:=F]
         res_dt[,qval:=NULL]
         res_dt <- res_dt[DA_by_sig|DA_by_struc_zero]
         
         # Summarizing the number of DA ASVs by condition
         Amount <- res_dt[,.(DA_ASVs = sum(DA_by_sig|DA_by_struc_zero)),
                          list(Genotype)]
         Amount[,":="(Compartment = X$meta_data$Compartment[1],
                      Soil = X$meta_data$Soil[1],
                      TOtal_ASVs = nrow(M$results)/4)]
         setcolorder(Amount, c("Genotype", "Compartment", "Soil",
                               "TOtal_ASVs", "DA_ASVs"))
         
         out <- list(res_mat = res_mat, res_dt = res_dt, 
                     Amount = Amount, mean_RA_res = mean_RA_res)
         
         return(out)
       }) -> DAA_results

DA_ASV_amounts <- DAA_results %>% lapply(function(x) x$Amount) %>% rbindlist()
# fwrite(DA_ASV_amounts, "ASV_DA_Overview_Lotus_seperate.csv")

DA_ASV_results <- DAA_results %>% lapply(function(x) x$res_dt) %>% rbindlist()
setcolorder(DA_ASV_results, c("ASV", "Soil", "Compartment", "Genotype"))
DA_ASV_results[,Genotype:=factor(Genotype, levels = c("WT", "symrk", "ccamk",
                                                      "nsp1", "nsp2"))]
DA_ASV_results <- DA_ASV_results[order(Soil, Compartment, Genotype, ASV)]
# fwrite(DA_ASV_results, "DA_results_Lotus_seperate.csv")

mean_RA_res <- DAA_results %>% 
  lapply(function(x) x$mean_RA_res) %>% 
  rbindlist()

# Heatmap ----------------------------------------------------------------------
res_mat_full <- t(Reduce("cbind", lapply(DAA_results, function(x) x$res_mat)))
res_mat_full <- data.table(Genotype = rownames(res_mat_full),
                           res_mat_full)

heatmap_data <- lapply(
  DAA_results,
  function(X){
    dt <- data.table(ASV_ID = rownames(X$res_mat), X$res_mat)
    dtt <- melt(dt, id.vars = 1, variable.name = "Genotype", value.name = "DAA")
    dtt[,":="(Compartment = X$res_dt$Compartment[1],
              Soil = X$res_dt$Soil[1])]
    return(dtt)
  }
)

heatmap_data <- rbindlist(heatmap_data)

# heatmap_data <- melt(res_mat_full, id.vars = 1, variable.name = "ASV_ID",
#                      value.name = "DAA")

heatmap_data[
  ,DAA:=fcase(DAA == 0, "NS",
              DAA == 1, "Enriched",
              DAA == -1, "Depleted")
]

## Removing low-abundance ASVs from heatmap ------------------------------------
k <- ncol(ASV_table_with_metadata)
RA_table <- ASV_table_with_metadata[,6:k]/rowSums(ASV_table_with_metadata[,6:k])
ASV_table_with_metadata[,6:k:=RA_table]
Mean_RA_by_condition <- ASV_table_with_metadata[,lapply(.SD, mean),
                                                by = list(Soil, Compartment, 
                                                          Genotype),
                                                .SDcols = colnames(RA_table)]
High_abn <- apply(Mean_RA_by_condition[,-(1:3)], 2, function(x) any(x > 0.005))
High_abn_ASV <- colnames(RA_table)[High_abn]

htmp_hiabn <- heatmap_data[ASV_ID %in% High_abn_ASV]

# Taxonomy annotation ----------------------------------------------------------
## Taxonomy ----
tax_bar <- taxonomy[Feature %in% unique(htmp_hiabn$ASV_ID)]
# tax_bar <- tax_bar[match(unique(htmp_hiabn$ASV_ID), Feature), c("Feature", "Order")]

colors_orders <- c(
  "Burkholderiales" = "#645394",
  "Caulobacterales" = "#8e3563", 
  "Chloroflexales" = "#CC99BB",
  "Flavobacteriales" = "#05294a",
  "Frankiales" = "#114477",
  "Gaiellales" = "#4477AA", 
  "Gemmatimonadales" = "#77AADD",
  "Micrococcales" = "#44AAAA",
  "Micromonosporales" = "#99D6DD",
  "Propionibacteriales" = "#117744",
  "Pseudomonadales" = "#88CCAA",
  "Pseudonocardiales" = "#95bb72",
  "Rhizobiales" = "#fdbb6b",
  "Sphingomonadales" = "lightyellow",
  "Streptomycetales" = "#fed5a4",
  "Xanthomonadales" = "#ffc0cb",
  "Unknown" = "darkgrey",
  "Other" = "lightgrey"
)

order_order <- setdiff(
  names(colors_orders), setdiff(names(colors_orders), unique(tax_bar$Order))
)
tax_bar[!(Order %in% names(colors_orders)), Order:="Other"]
tax_bar[,Order:=factor(Order, levels = order_order)]
tax_bar <- tax_bar[order(Order)]
htmp_hiabn <- htmp_hiabn[,ASV_ID:=factor(ASV_ID, levels = tax_bar$Feature)]
htmp_hiabn <- htmp_hiabn[order(ASV_ID)]
tax_bar$Feature <- factor(tax_bar$Feature, levels = tax_bar$Feature)
p_tax <- ggplot(tax_bar, aes(x=Feature, y=1, fill=Order)) +
  geom_tile() +
  scale_fill_manual(values=colors_orders) +
  theme_void() +
  labs(fill = "Bacterial order") +
  theme(legend.position="none",
        legend.text = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        legend.key.size = unit(0.25, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "lines")
  )+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))+
  NULL
p_tax

# Heatmap ----------------------------------------------------------------------
heatmap <- ggplot(data = htmp_hiabn, 
                  mapping = aes(x = ASV_ID, y = Genotype, fill = DAA))+
  geom_tile()+
  facet_wrap2(vars(Compartment, Soil), strip = strip_nested(), ncol = 1,
              strip.position = "left")+
  # facet_wrap(~Soil + Compartment, strip.position = "left", ncol = 1)+
  scale_fill_manual(values = c("darkblue", "#902121", "white"),
                    breaks = c("Depleted", "Enriched", "NS"))+
  labs(y = "Lotus")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2,'lines'),
    strip.background = element_rect(colour = NA),
    strip.placement = "outside",
    axis.title.y = element_text(size = 8, family = "Helvetica",
                                colour = "black", face = "bold"),
    axis.text.y = element_text(size = 8, family = "Helvetica",
                               colour = "black"),
    legend.text = element_text(size = 8, family = "Helvetica"),
    legend.title = element_text(size = 8, family = "Helvetica"),
    strip.text = element_text(size = 8, family = "Helvetica",
                              face = "bold"),
    plot.margin = margin(t = 0, r = 0, b = 0.5, l = 0, unit = "lines")
  )+
  NULL; heatmap

# Annotations ------------------------------------------------------------------
## Barplot ----
bar_plot <- ggplot(data = mean_RA_res, mapping = aes(x = mean_RAs, y = Genotype))+
  geom_bar(stat = "identity")+
  facet_wrap2(vars(Compartment, Soil), strip = strip_nested(), ncol = 1,
              strip.position = "left")+
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.1, 0.2))+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.2,'lines'),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, family = "Helvetica",
                               colour = "black"),
    legend.text = element_text(size = 8, family = "Helvetica"),
    legend.title = element_text(size = 8, family = "Helvetica"),
    strip.text = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0.5, l = 0, unit = "lines")
  )+
  NULL; bar_plot

p_axis_title <- ggplot() +
  theme_void() +
  geom_text(aes(x = -10, y = 0, label = "RA in WT"),
            fontface = "bold", size = 8/.pt) +
  xlab(NULL) + ylab(NULL)

top_plot <- (p_tax | p_axis_title) + plot_layout(widths = c(0.9, 0.1))
bottom_plot <- (heatmap | bar_plot) + plot_layout(widths = c(0.9, 0.1))
main_plot <- (top_plot/bottom_plot) + plot_layout(heights = c(0.05, 0.95))
main_plot
