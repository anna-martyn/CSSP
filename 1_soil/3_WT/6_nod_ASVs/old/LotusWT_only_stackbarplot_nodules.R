options(warn=-1)

# cleanup
rm(list=ls())

# directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# files
design.file <- paste("Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", sep="")
taxonomy.file <- paste("LotusSep_exclUFnew_10_4_silva138_taxonomy.txt", sep="")
otu_table.file <- paste("LotusSep_exclUFnew_ASVtable_10_4_nospike.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

#modify taxonomy table
library(data.table)
library(magrittr)
library(stringr)
library(dplyr) 
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggh4x)
library(tidyverse)

taxa <- lapply(as.list(taxonomy$Taxon),
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
taxonomy <- data.table(ASV = taxonomy$`ASVid`, taxa, Confidence = taxonomy$Confidence)


# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- rownames(otu_table) %in% taxonomy$ASV
otu_table <- otu_table[idx, ]

idx <- match(rownames(otu_table), taxonomy$ASV)
taxonomy <- taxonomy[idx, ]

# otu_table normalization
otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))
df <- as.data.frame(otu_table_norm)

# add taxonomy info to "df"
tax_order <- taxonomy[ ,5]
tax_family <- taxonomy[ ,6]
tax_genus <- taxonomy[ ,7]
df$Order <- tax_order
df$Family <- tax_family
df$Genus <- tax_genus
ASVid <- row.names(df)
df <- cbind(ASVid, df, df, df)
row.names(df) <- NULL

# Only select nodule data and split up into NPK. PK and UF dataset.
df_nod <- df[,c(1,4,7,10,13,16,19,73,76,79,82,85,88,140,143,146,149,152,155)]
df_nod_NPK <- df_nod[,1:7]
df_nod_PK <- df_nod[,c(1,8:13)]
df_nod_UF <- df_nod[,c(1,14:19)]

# Add average value in additional column.
df_nod_NPK$Average <- rowMeans(df_nod_NPK[, c(2:7)])
df_nod_PK$Average <- rowMeans(df_nod_PK[, c(2:7)])
df_nod_UF$Average <- rowMeans(df_nod_UF[, c(2:7)])

# Check top 20 ASVs in nodules in each soil type.
NPK_sorted <- df_nod_NPK[order(-df_nod_NPK$Average), ]
NPK_top20 <- NPK_sorted[1:20, c("ASVid", "Average")]

PK_sorted <- df_nod_PK[order(-df_nod_PK$Average), ]
PK_top20 <- PK_sorted[1:20, c("ASVid", "Average")]

UF_sorted <- df_nod_UF[order(-df_nod_UF$Average), ]
UF_top20 <- UF_sorted[1:20, c("ASVid", "Average")]

combined_df <- Reduce(function(x, y) merge(x, y, by = "ASVid", all = TRUE), list(NPK_top20, PK_top20,UF_top20))

ASVs_top20all <- combined_df$ASVid

# Filter each dataframe to only include rows with ASVids in ASVs_top20all
NPK_top20all <- NPK_sorted[NPK_sorted$ASVid %in% ASVs_top20all, ]
NPK_top20all <- NPK_top20all[,c("ASVid","Average")]
colnames(NPK_top20all)[which(colnames(NPK_top20all) == "Average")] <- "Average_NPK"

PK_top20all <- PK_sorted[PK_sorted$ASVid %in% ASVs_top20all, ]
PK_top20all <- PK_top20all[,c("ASVid","Average")]
colnames(PK_top20all)[which(colnames(PK_top20all) == "Average")] <- "Average_PK"

UF_top20all <- UF_sorted[UF_sorted$ASVid %in% ASVs_top20all, ]
UF_top20all <- UF_top20all[,c("ASVid","Average")]
colnames(UF_top20all)[which(colnames(UF_top20all) == "Average")] <- "Average_UF"

# Merge the filtered dataframes by ASVid, ensuring all ASVs_top20all are retained
combined_df <- merge(
  merge(
    data.frame(ASVid = ASVs_top20all), NPK_top20all, by = "ASVid", all.x = TRUE
  ),
  merge(PK_top20all, UF_top20all, by = "ASVid", all.x = TRUE),
  by = "ASVid", all.x = TRUE
)

# Now add taxonomy information.
filtered_taxonomy <- taxonomy[taxonomy$ASV %in% ASVs_top20all, ]
colnames(filtered_taxonomy)[which(colnames(filtered_taxonomy) == "ASV")] <- "ASVid"
topASVs_nod_final <- merge(combined_df, filtered_taxonomy, by = "ASVid", all.y = TRUE)

topASVs_nod_final_1perc <- topASVs_nod_final %>%
  filter(if_any(c(Average_NPK, Average_PK, Average_UF), ~ . >= 0.01))

# stacked bar plot

# Preliminary renaming of ASVs.
prel_ASVs <- c("Unknown_1","Mesorhizobium_1","Pseudomonas_1","Mesorhizobium_2","Mesorhizobium_3","Mesorhizobium_4","Mesorhizobium_5","Mesorhizobium_6","Unknown_2")

topASVs_nod_final_1perc$Isolate <- prel_ASVs

df_long_final <- topASVs_nod_final_1perc %>%
  pivot_longer(cols = starts_with("Average_"),  # Select columns to pivot
               names_to = "Soil_type",         # Renamed to Soil_type
               values_to = "RA") %>%
  mutate(Soil_type = str_replace(Soil_type, "Average_", ""))  # Remove "Average_"


### visualization by ggplot2

levels(df_long_final$Soil_type)<- c("NPK","PK","UF") ### This is to make the letter italic for the genotype name

df_long_final$Soil_type <- factor(df_long_final$Soil_type, levels = c("NPK","PK","UF"))

isolate_levels <- c("Pseudomonas_1","Mesorhizobium_1","Mesorhizobium_2","Mesorhizobium_3","Mesorhizobium_4","Mesorhizobium_5","Mesorhizobium_6","Unknown_1","Unknown_2")
df_long_final$Isolate <- factor(df_long_final$Isolate, levels = isolate_levels)

# Define colour scale. 
## option 1: grey
# grey_palette <- gray.colors(n = 9, start = 0.9, end = 0.2)
# names(grey_palette) <- isolate_levels
# print(grey_palette)


## option 2: yellow (all)
# colors <- data.frame(group=c("Unknown_1","Mesorhizobium_1","Pseudomonas_1","Mesorhizobium_2","Mesorhizobium_3","Mesorhizobium_4","Mesorhizobium_5","Mesorhizobium_6","Unknown_2"), 
#                       colors=c("#fff9e7","#ffeeb6","#ffe385","#ffdd6c","#FFD95A","#FFD700","#FFCC00","#FFB300","#FFA500"))

## option 2: yellow, turquoise, grey
colors <- data.frame(group=c("Pseudomonas_1","Mesorhizobium_1","Mesorhizobium_2","Mesorhizobium_3","Mesorhizobium_4","Mesorhizobium_5","Mesorhizobium_6","Unknown_1","Unknown_2"), 
                     colors=c("#88CCAA","#fff9e7","#ffeeb6","#ffe385","#ffdd6c","#FFCC00","#FFB300","lightgrey","darkgrey"))

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=25),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

p2 <- ggplot(df_long_final, aes(x=Soil_type, y = RA, fill = Isolate)) +
  geom_bar(stat = "identity", width = .5) +
  #scale_fill_manual(values=grey_palette)+
  scale_fill_manual(values=as.character(colors$color)) +
  main_theme+
  ggtitle("Nodule ASVs across three soil types")+
  ylab("Relative abundance")+
  theme(legend.position = "right")+ guides(fill=guide_legend(nrow=9))+ 
  theme(axis.text.x = element_text(size=25),
        plot.title = element_text(size = 25, face = "bold"),  # Adjust the size here
        strip.text.x = element_text(size = 25, face = "bold"),
        legend.text=element_text(size=25),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=25))


p2

ggsave(paste("Lotus_WT_nodules_RA.png", sep=""), p2, width=8, height=4)
