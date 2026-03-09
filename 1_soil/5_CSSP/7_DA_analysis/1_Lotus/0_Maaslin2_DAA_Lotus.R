### Genotype-specific differential abundance analysis using MaAsLin2 and
### structural zero detection, followed by visualisation.

# Clean up.
options(warn=-1)
rm(list=ls())

# Load packages and data -------------------------------------------------------
## Load required packages.
pkg <- c("data.table", "magrittr", "ggplot2", "vegan", "Maaslin2", "patchwork",
         "colorRamp2", "ggh4x", "ggtext", "tidyr", "dplyr")

for(pk in pkg){
  library(pk, character.only = T)
}

## Set directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load custom source file for later analysis (includes structural zero function).
source("../0_files/Structural_zeros.R")

## Read and convert input data.
ASV_table <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
colnames(ASV_table)[1] <- "ASVid"
ASV_table_full <- data.table(ASV_table)

meta_data <- fread("../../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt",
                   drop = c(5, 7, 8))

taxonomy <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv"
)
rename_tax <- function(tax_table){
  colnames(tax_table)[colnames(tax_table) == "Feature ID"] <- "ASVid"
  tax_table %>%
    separate(Taxon, into = c("Kingdom","Phylum","Class","Order",
                             "Family","Genus","Species"),
             sep = "; ", fill = "right") %>%
    mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
    replace(is.na(.), "Unknown") %>%
    select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus,
           Species, Confidence)
}
taxonomy <- rename_tax(taxonomy)
taxonomy <- as.data.table(taxonomy)

## Merge ASV table and metadata info. 
ASV_table_with_metadata <- transpose(ASV_table,
                                     make.names = "ASVid",
                                     keep.names = "SampleID")
ASV_table_with_metadata <- merge(meta_data,
                                 ASV_table_with_metadata,
                                 by = "SampleID")

## Set reference levels for the genotypes and define experimental conditions.
ASV_table_with_metadata[
  ,Genotype:=factor(Genotype, levels = c("WT", "symrk", "ccamk", 
                                         "nsp1", "nsp2"))
]

Experimental_conditions <- as.matrix(
  expand.grid(Compartment = c("Root", "Rhizosphere"),
              Soil = c("NPK", "PK", "UF"))
)

# Differential abundance analysis with MaAsLin2 --------------------------------
## Subset data for specific compartment-soil combination and prepare input for MaAsLin2.
ASV_tables_subsets <- apply(
  Experimental_conditions, 1,
  function(x){
    # Subset data for given compartment and soil condition.
    data_subset <- ASV_table_with_metadata[Compartment == x[1] & Soil == x[2]]
    
    # Extract metadata and drop unused factor levels.
    meta_data <- data.frame(data_subset[,2:5],
                            row.names = data_subset$SampleID)
    meta_data$Genotype <- droplevels(meta_data$Genotype)
    
    # Extract and transpose ASV table for MaAsLin2 input.
    ASV_table <- t(data_subset[,-(1:5)])
    colnames(ASV_table) <- data_subset$SampleID
    
    # Return ASV and metadata as a list.
    Out <- list(ASV_table = ASV_table, meta_data = meta_data)
    return(Out)
  }
)

## Run MaAsLin2 and structural zero analysis for each subset.
lapply(ASV_tables_subsets,
       function(X){
         # Ensure output folder exists.
         if(!dir.exists("Maaslin2 files")) dir.create("Maaslin2 files")
         # Run MaAsLin2.
         M <- Maaslin2(input_data = X$ASV_table,
                       input_metadata = X$meta_data,
                       output = "Maaslin2 files",
                       max_significance = 0.05,
                       plot_heatmap = F,
                       plot_scatter = F,
                       fixed_effects = "Genotype",
                       min_prevalence = 0.1,
                       reference = "WT")
         res_dt <- data.table(M$results)
         res_dt[,feature:=gsub("X", "", feature)]
         res_dt <- res_dt[order(feature)]
         
         # Run structural zeros function.
         X$meta_data$library_size <- colSums(X$ASV_table)
         
         S <- Structural_zeros3(X$ASV_table, X$meta_data,
                                group = "Genotype", ref = "WT",
                                min_reads = 20, min_present_reps = 3)
         
         # Construct a binary matrix to indicate which ASVs are DA.
         res_mat <- matrix(0, nrow = nrow(X$ASV_table), ncol = 4)
         rownames(res_mat) <- sort(rownames(X$ASV_table))
         colnames(res_mat) <- c("symrk", "ccamk", "nsp1", "nsp2")
         
         DA_symrk <- res_dt[value == "symrk" & qval < 0.05]$feature
         res_mat[DA_symrk,"symrk"] <- sign(res_dt[value == "symrk" &
                                                    qval < 0.05]$coef)
         DA_ccamk <- res_dt[value == "ccamk" & qval < 0.05]$feature
         res_mat[DA_ccamk,"ccamk"] <- sign(res_dt[value == "ccamk" &
                                                    qval < 0.05]$coef)
         DA_nsp1 <- res_dt[value == "nsp1" & qval < 0.05]$feature
         res_mat[DA_nsp1,"nsp1"] <- sign(res_dt[value == "nsp1" &
                                                  qval < 0.05]$coef)
         DA_nsp2 <- res_dt[value == "nsp2" & qval < 0.05]$feature
         res_mat[DA_nsp2,"nsp2"] <- sign(res_dt[value == "nsp2" &
                                                  qval < 0.05]$coef)
         # Integrate structural zero results.
         all_struc_zero_DA_ASVS <- unlist(S$struc_zero_DA)
         N_DA_ASV_by_genotype <- unlist(lapply(S$struc_zero_DA, length))
         names(all_struc_zero_DA_ASVS) <- rep(c("symrk", "ccamk",
                                                "nsp1", "nsp2"),
                                              N_DA_ASV_by_genotype)
                                              
         struc_zero_res_dt <- data.table(
           ASV = all_struc_zero_DA_ASVS,
           Genotype = names(all_struc_zero_DA_ASVS),
           DA_by_struc_zero = T
         )
         
         # Fill in structural zeros results from symrk.
         struc_zero_symrk <- setdiff(
           S$struc_zero_DA$symrk, res_dt[value == "symrk" & qval < 0.05]$feature
         )
         S_symrk <- S$struc_zero_table[struc_zero_symrk, c("WT", "symrk")]
         struc_zero_FC_symrk <- S_symrk[,"WT"] - S_symrk[,"symrk"]
         res_mat[names(struc_zero_FC_symrk), "symrk"] <- struc_zero_FC_symrk
         
         # Fill in structural zeros results from ccamk.
         struc_zero_ccamk <- setdiff(
           S$struc_zero_DA$ccamk, res_dt[value == "ccamk" & qval < 0.05]$feature
         )
         S_ccamk <- S$struc_zero_table[struc_zero_ccamk, c("WT", "ccamk")]
         struc_zero_FC_ccamk <- S_ccamk[,"WT"] - S_ccamk[,"ccamk"]
         res_mat[names(struc_zero_FC_ccamk), "ccamk"] <- struc_zero_FC_ccamk
         
         # Fill in structural zeros results from nsp1.
         struc_zero_nsp1 <- setdiff(
           S$struc_zero_DA$nsp1, res_dt[value == "nsp1" & qval < 0.05]$feature
         )
         S_nsp1 <- S$struc_zero_table[struc_zero_nsp1, c("WT", "nsp1")]
         struc_zero_FC_nsp1 <- S_nsp1[,"WT"] - S_nsp1[,"nsp1"]
         res_mat[names(struc_zero_FC_nsp1), "nsp1"] <- struc_zero_FC_nsp1
         
         # Fill in structural zeros results from nsp2.
         struc_zero_nsp2 <- setdiff(
           S$struc_zero_DA$nsp2, res_dt[value == "nsp2" & qval < 0.05]$feature
         )
         S_nsp2 <- S$struc_zero_table[struc_zero_nsp2, c("WT", "nsp2")]
         struc_zero_FC_nsp2 <- S_nsp2[,"WT"] - S_nsp2[,"nsp2"]
         res_mat[names(struc_zero_FC_nsp2), "nsp2"] <- struc_zero_FC_nsp2
         
         # Calculate mean RA of DA ASVs in WT samples.
         WT_samples <- rownames(X$meta_data)[X$meta_data$Genotype == "WT"]
         lib_size <- colSums(X$ASV_table[,WT_samples])
         mean_RAs <- apply(
           res_mat, 2,
           function(y){
             ASVs <- rownames(res_mat)[y != 0]
             
             ASV_subset <- X$ASV_table[ASVs,WT_samples]
             mean(rowSums(t(ASV_subset)/lib_size))
           }
         )
         
         mean_RA_res <- data.table(Compartment = X$meta_data$Compartment[1],
                                   Soil = X$meta_data$Soil[1],
                                   Genotype = names(mean_RAs),
                                   mean_RAs)
         
         # Summarise and combine results.
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
         
         # Count number of DA ASVs per condition.
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

## Export summary tables for results.
DA_ASV_amounts <- DAA_results %>% lapply(function(x) x$Amount) %>% rbindlist()
fwrite(DA_ASV_amounts, "Lotus_DA_ASVs_overview.csv")

DA_ASV_results <- DAA_results %>% lapply(function(x) x$res_dt) %>% rbindlist()
setcolorder(DA_ASV_results, c("ASV", "Soil", "Compartment", "Genotype"))
DA_ASV_results[,Genotype:=factor(Genotype, levels = c("WT", "symrk", "ccamk",
                                                      "nsp1", "nsp2"))]
DA_ASV_results <- DA_ASV_results[order(Soil, Compartment, Genotype, ASV)]
fwrite(DA_ASV_results, "Lotus_DAanalysis_results.csv")

mean_RA_res <- DAA_results %>% 
  lapply(function(x) x$mean_RA_res) %>% 
  rbindlist()

# Prepare data for heatmap visualisation----------------------------------------
## Combine DA matrices across all conditions.
res_mat_full <- t(Reduce("cbind", lapply(DAA_results, function(x) x$res_mat)))
res_mat_full <- data.table(Genotype = rownames(res_mat_full),
                           res_mat_full)
## Convert DA matrices to long format for plotting.
heatmap_data <- lapply(
  DAA_results,
  function(X){
    dt <- data.table(ASVid = rownames(X$res_mat), X$res_mat)
    dtt <- melt(dt, id.vars = 1, variable.name = "Genotype", value.name = "DAA")
    dtt[,":="(Compartment = X$res_dt$Compartment[1],
              Soil = X$res_dt$Soil[1])]
    return(dtt)
  }
)

heatmap_data <- rbindlist(heatmap_data)

## Convert numeric DA indicators to categorical variables.
heatmap_data[
  ,DAA:=fcase(DAA == 0, "NS",
              DAA == 1, "Enriched",
              DAA == -1, "Depleted")
]

## Keep only DA ASVs in heatmap.
DA_ASVs <- unique(heatmap_data[DAA != "NS"]$ASVid)
htmp_hiabn <- heatmap_data[ASVid %in% DA_ASVs]

# Taxonomy annotation for heatmap ----------------------------------------------
## Subset taxonomy table to include only DA ASVs.
tax_bar <- taxonomy[ASVid %in% unique(htmp_hiabn$ASVid)]

## Define colours for bacterial orders displayed in figure.
colors_orders <- fread("../../../../0_files/Bacterial_order_colors.csv")

# Load orders to display. These are top 20 mean RA orders in either Lotus or Hordeum
# WT across compartments and soils (as used for supplementary figures).
combined_top_orders <- readRDS("../Orders_to_display.rds")

## Define orders without assigned colours to 'Other'.
tax_bar[!(Order %in% combined_top_orders), Order:="Other"]
tax_bar[,Order:=factor(Order, levels = colors_orders$Order)]
tax_bar <- tax_bar[order(Order)]

## Match ASV ordering between taxonomy and heatmap.
htmp_hiabn <- htmp_hiabn[,ASVid:=factor(ASVid, levels = tax_bar$ASVid)]
htmp_hiabn <- htmp_hiabn[order(ASVid)]
tax_bar$ASVid <- factor(tax_bar$ASVid, levels = tax_bar$ASVid)

## Make the taxonomy barplot.
p_tax <- ggplot(tax_bar, aes(x = ASVid, y=1, fill = Order)) +
  geom_tile() +
  scale_fill_manual(values = colors_orders$Color, breaks = colors_orders$Order, drop = FALSE) +
  theme_void() +
  labs(fill = "Bacterial order") +
  theme(legend.position="none",
        legend.text = element_text(color = "black", size = 6),
        legend.title = element_text(color = "black", size = 6, face = "bold"),
        legend.key.size = unit(0.25, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "lines")
  )+
  NULL
p_tax

# Heatmap for differential abundance --------------------------------------
## Merge DA counts for y-axis annotation.
htmp_hiabn <- merge(
  htmp_hiabn, DA_ASV_amounts, by = c("Compartment", "Soil", "Genotype")
)

## Format y-axis labels to include number of DA ASVs.
htmp_hiabn[
  ,y_axis:=paste0("*", Genotype, "*", " (", DA_ASVs, ")")
]

## Define genotype ordering.
ordering <- c(
  unique(htmp_hiabn[grepl("nsp2", y_axis)]$y_axis),
  unique(htmp_hiabn[grepl("nsp1", y_axis)]$y_axis),
  unique(htmp_hiabn[grepl("ccamk", y_axis)]$y_axis),
  unique(htmp_hiabn[grepl("symrk", y_axis)]$y_axis)
)

htmp_hiabn[,y_axis:=factor(y_axis, levels = ordering)]

## Generate the heatmap.
heatmap <- ggplot(data = htmp_hiabn, 
                  mapping = aes(x = ASVid, y = y_axis, fill = DAA))+
  geom_tile()+
  facet_wrap2(vars(Compartment, Soil), strip = strip_nested(), ncol = 1,
              strip.position = "left", scales = "free_y")+
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
    axis.title.y = element_text(size = 6, family = "Helvetica",
                                colour = "black", face = "bold"),
    axis.text.y = element_markdown(size = 6, family = "Helvetica",
                                   colour = "black"),
    legend.text = element_text(size = 6, family = "Helvetica"),
    legend.title = element_text(size = 6, family = "Helvetica"),
    strip.text = element_text(size = 6, family = "Helvetica",
                              face = "bold"),
    plot.margin = margin(t = 0, r = 0, b = 0.5, l = 0, unit = "lines")
  )+
  NULL; heatmap

# Barplot of cumulative relative abundance of DA ASVs in WT --------------------
mean_RA_res[
  ,Genotype:=factor(Genotype, levels = c("nsp2", "nsp1", "ccamk", "symrk"))
]
bar_plot <- ggplot(data = mean_RA_res, 
                   mapping = aes(x = mean_RAs, y = Genotype))+
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
    axis.text.x = element_text(size = 6, family = "Helvetica",
                               colour = "black"),
    legend.text = element_text(size = 6, family = "Helvetica"),
    legend.title = element_text(size = 6, family = "Helvetica"),
    strip.text = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0.5, l = 0, unit = "lines")
  )+
  NULL; bar_plot

## Make an axis title panel.
p_axis_title <- ggplot() +
  theme_void() +
  geom_text(aes(x = -10, y = 0, label = "RA in WT"),
            fontface = "bold", size = 6/.pt) +
  xlab(NULL) + ylab(NULL)

# Save output files ------------------------------------------------------------
saveRDS(heatmap, "LotusCSSP_Askov_DA_heatmap.rds")
saveRDS(bar_plot, "LotusCSSP_Askov_DA_barplot.rds")
saveRDS(p_tax, "LotusCSSP_Askov_DA_taxonomy.rds")
saveRDS(p_axis_title, "LotusCSSP_Askov_DA_axis_title.rds")
