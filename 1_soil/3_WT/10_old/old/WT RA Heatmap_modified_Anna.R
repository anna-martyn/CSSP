# Code by Ib Thorsgaard Jensen
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Setup ----
# Loading packages
pkg <- c("data.table", "ggplot2", "magrittr", "ComplexHeatmap", "FSA", "multcompView", "colorRamp2")
for(i in pkg){
  library(i, character.only = T)
}

# Loading files
design_Lotus <- fread("Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt")
otu_table_Lotus <- fread("LotusSep_exclUFnew_ASVtable_10_4_nospike.txt")
taxonomy_Lotus <- fread("LotusSep_exclUFnew_10_4_silva138_taxonomy.txt")

design_Barley <- fread("BarleyCSSP_Askov_reseq_metadata.txt")
colnames(design_Barley)[1] <- "SampleID"
otu_table_Barley <- fread("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt")
taxonomy_Barley <- fread("Barley_Askov_Rep_10_4_taxonomy.txt")

#colors <- c(NPK = "#1b9e77", PK = "#d95f02", UF = "#666666")
colors <- c(NPK = "#6D3B00", PK = "#A06A37", UF = "#D2B48C")

# Taxonomy - Lotus
taxa <- lapply(as.list(taxonomy_Lotus$Taxon),
               function(x) x %>% strsplit(split = "; ") %>% unlist() %>% substr(start = 4, stop = 100))

taxa <- lapply(taxa, function(x) c( x, rep("Unknown", 7 - length(x)) ))
taxa <- lapply(taxa, function(x) 
  data.table(Kingdom = x[1],
             Phylum = x[2],
             Class = x[3],
             Order = x[4],
             Family = x[5],
             Genus = x[6],
             Species = x[7])
)
taxa <- rbindlist(taxa)
taxonomy_Lotus <- data.table(ASVid = taxonomy_Lotus$`ASVid`, taxa, Confidence = taxonomy_Lotus$Confidence)

# Taxonomy - Barley
taxa <- lapply(as.list(taxonomy_Barley$Taxon),
               function(x) x %>% strsplit(split = "; ") %>% unlist() %>% substr(start = 4, stop = 100))

taxa <- lapply(taxa, function(x) c( x, rep("Unknown", 7 - length(x)) ))
taxa <- lapply(taxa, function(x) 
  data.table(Kingdom = x[1],
             Phylum = x[2],
             Class = x[3],
             Order = x[4],
             Family = x[5],
             Genus = x[6],
             Species = x[7])
)
taxa <- rbindlist(taxa)
taxonomy_Barley <- data.table(ASVid = taxonomy_Barley$`ASVid`, taxa, Confidence = taxonomy_Barley$Confidence)

# Re-order columns
design_Lotus <- design_Lotus[SampleID %in% colnames(otu_table_Lotus)]
design_Lotus <- design_Lotus[Genotype == "WT" & Compartment != "nodules"]
design_Lotus[Compartment == "root", Compartment:="Root"]
design_Lotus[Compartment == "rhizosphere", Compartment:="Rhizosphere"]
otu_table_Lotus <- otu_table_Lotus[,c("ASVid", design_Lotus$SampleID), with = F]

design_Barley <- design_Barley[SampleID %in% colnames(otu_table_Barley)]
design_Barley <- design_Barley[Genotype == "WT"]
design_Barley[Compartment == "root", Compartment:="Root"]
design_Barley[Compartment == "rhizosphere", Compartment:="Rhizosphere"]
otu_table_Barley <- otu_table_Barley[,c("ASVid", design_Barley$SampleID), with = F]
N_reads <- colSums(otu_table_Barley[,-1])
sample_keep <- names(N_reads)[N_reads > 1000]
design_Barley <- design_Barley[SampleID %in% sample_keep]
otu_table_Barley <- otu_table_Barley[,c("ASVid", sample_keep), with = F]

# sanity checks
all(design_Lotus$SampleID == colnames(otu_table_Lotus)[-1])
all(otu_table_Lotus$ASVid == taxonomy_Lotus$ASV)

all(design_Barley$SampleID == colnames(otu_table_Barley)[-1])
all(otu_table_Barley$ASVid == taxonomy_Barley$ASV)

# Transformations and aggregation ----
# TSS transformation
otu_table_Lotus_RA <- data.table(otu_table_Lotus)
otu_table_Lotus_RA[,2:ncol(otu_table_Lotus):=lapply(.SD, function(x) x/sum(x)), .SDcols = 2:ncol(otu_table_Lotus)]

otu_table_Barley_RA <- data.table(otu_table_Barley)
otu_table_Barley_RA[,2:ncol(otu_table_Barley):=lapply(.SD, function(x) x/sum(x)), .SDcols = 2:ncol(otu_table_Barley)]

# Aggregating to order level
Order_table_Lotus_RA <- merge(otu_table_Lotus_RA, taxonomy_Lotus[,c("ASVid", "Order")])
setcolorder(Order_table_Lotus_RA, "Order")
Order_table_Lotus_RA[,ASVid:=NULL]
Order_table_Lotus_RA <- Order_table_Lotus_RA[,lapply(.SD, sum), Order]

Order_table_Barley_RA <- merge(otu_table_Barley_RA, taxonomy_Barley[,c("ASVid", "Order")])
setcolorder(Order_table_Barley_RA, "Order")
Order_table_Barley_RA[,ASVid:=NULL]
Order_table_Barley_RA <- Order_table_Barley_RA[,lapply(.SD, sum), Order]

# Subsetting for only orders of interest
orders_to_show <- c("Burkholderiales", "Caulobacterales", "Chloroflexales", "Flavobacteriales",
                    "Frankiales", "Gaiellales", "Gemmatimonadales", "Micrococcales",
                    "Micromonosporales", "Propionibacteriales", "Pseudomonadales",
                    "Pseudonocardiales", "Rhizobiales", "Sphingomonadales",
                    "Streptomycetales","Unknown","Xanthomonadales")

Order_table_Lotus_RA <- Order_table_Lotus_RA[Order %in% orders_to_show]
Order_table_Barley_RA <- Order_table_Barley_RA[Order %in% orders_to_show]
Order_table_RA <- merge(Order_table_Lotus_RA, Order_table_Barley_RA, by = "Order")

# Heatmap ----
df_order <- data.frame(Order_table_RA[,-1], row.names = Order_table_RA$Order)
rbind(design_Lotus[,c("SampleID", "Plant", "Compartment", "Soil")],
      design_Barley[,c("SampleID", "Plant", "Compartment", "Soil")]) -> design

anno_df <- data.frame(design[,c("Plant","Compartment", "Soil")], row.names = design$SampleID)
ha <- HeatmapAnnotation(df = anno_df,
                        col = list(Plant = c(Lotus = "#5C5C5C", Barley = "grey"),
                                   Compartment = c(Root = "#90ee90", Rhizosphere = "#6a4a3a"),
                                   Soil = colors),
                        show_annotation_name = F)
# ha <- HeatmapAnnotation(Compartment = anno_block(gp = gpar(fill = c("gold", "black", "magenta")),
#                                                  labels = c("root", "nodules", "rhizosphere"),
#                                                  labels_gp = gpar(col = "white", fontsize = 10)),
#                         Soil = anno_block(gp = gpar(fill = colors),
#                                           labels = names(colors),
#                                           labels_gp = gpar(col = "white", fontsize = 10)) )
# col_fun = colorRamp2(c(0, 0.005, 0.04, 0.15), c("white", "lightgreen", "darkgreen", "blue"))
col_fun = colorRamp2(c(0, 0.005, 0.052, 0.052001, 0.15999, 0.16, 0.34, 0.64),
                     c("darkblue", "deepskyblue", "white","#feedb4", "#F09D00", "#FFB3B2", "#ba0319", "#902121"))
# col_fun = colorRamp2(c(0, 0.04), c("white", "darkgreen"))

# [-c(1,3,4,15,16),]

empty_anno = rowAnnotation(foo = anno_empty(border = FALSE),
                           foo2 = anno_empty(border = FALSE))

anno_df$Plant <- factor(anno_df$Plant, levels = c("Lotus","Barley"))

Heatmap(matrix = as.matrix(df_order),
        border_gp = gpar(col = "black", lty = 1),
        show_heatmap_legend = F,
        col = col_fun,
        column_title = NULL,
        # name = "RA",
        # heatmap_legend_param = list(legend_direction = "horizontal"),
        column_split = anno_df,
        column_gap = unit(0, "mm"),
        # row_title = "Rhizosphere",
        # row_names_side = "left",
        cluster_columns = F,
        cluster_rows = F,
        show_column_names = F,
        top_annotation = ha,
        left_annotation = empty_anno) -> H

lgd = Legend(col_fun = col_fun, title = "RA",
              at = seq(0.16, 0.64, 0.08))

lgd2 = Legend(col_fun = col_fun,
             at = seq(0.06, 0.14, 0.02))

lgd3 = Legend(col_fun = col_fun,
             at = seq(0, 0.05, 0.01))

pdf("WT_Heatmap_2.pdf", width = 180/25.4, height = 180/25.4) 
draw(H)
draw(lgd,
     x = unit(0.09, "npc"),
     y = unit(0.87, "npc"),
     just = c("right", "top"))

draw(lgd2,
     x = unit(0.09, "npc"),
     y = unit(0.54, "npc"),
     just = c("right", "top"))

draw(lgd3, 
     x = unit(0.09, "npc"), 
     y = unit(0.33, "npc"), 
     just = c("right", "top"))
dev.off()


# df_transposed <- t(df_order)
# df_transposed <- as.data.frame(df_transposed)
# 
# 
# Heatmap(matrix = as.matrix(df_transposed),
#         border_gp = gpar(col = "black", lty = 1),
#         show_heatmap_legend = F,
#         col = col_fun,
#         row_title = NULL,
#         # name = "RA",
#         # heatmap_legend_param = list(legend_direction = "horizontal"),
#         row_split = anno_df,
#         row_gap = unit(0, "mm"),
#         # row_title = "Rhizosphere",
#         # row_names_side = "left",
#         cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         top_annotation = ha,
#         left_annotation = empty_anno) -> H2

