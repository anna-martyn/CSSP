
options(warn=-1)

# cleanup
rm(list=ls())

# directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# files
design.file <- paste("BarleyCSSP_Askov_reseq_metadata.txt", sep="")
taxonomy.file <- paste("Barley_Askov_Rep_10_4_taxonomy.txt", sep="")
otu_table.file <- paste("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

#modify taxonomy table
library(data.table)
library(magrittr)

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

#filter out Lotus UF data
design <- design[1:188, ]

# re-order data matrices
idx <- design$Sample_ID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$Sample_ID, colnames(otu_table))
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
df$Order <- tax_order
ASVid <- row.names(df)
df <- cbind(ASVid, df)
row.names(df) <- NULL

### reshape data to long dataframe
df.long <- reshape(df, varying = 2:189, direction = "long", idvar = 'ASVid', timevar="sampleID", v.names=c("RA"), 
times=c(colnames(df[,2:189])), sep="")

row.names(df.long) <- NULL

df.long <- as.data.frame(df.long)

###  make a dataframe that sums RA by tax and sampleID
library(dplyr)

df.long_order <- df.long %>% select (2:4) ### select the column 2,3,4 from df.long dataframe

df.long_order <- df.long_order %>%
  group_by_(.dots = c("Order", "sampleID")) %>%   ### this is to summarize RA by both tax and sampleID
  summarise_all(funs(sum))

### Add genotype information for df.long 
df2 <- data.frame(sampleID = design$Sample_ID, Soil=design$Soil, Genotype=design$Genotype, Compartment=design$Compartment)

df.long_order <- full_join(df.long_order, df2, by="sampleID")

#subset for soil samples
idx <- df.long_order$Genotype%in% c("WT")
df.long_order_WT <- df.long_order[idx,]

### calculate mean RA for each taxa across samples and show top20
mean_RA_order_WT <- df.long_order_WT %>%
  group_by(Order) %>%
  dplyr::summarize(Mean = mean(RA, na.rm = TRUE))

mean_RA_order_WT <- as.data.frame(mean_RA_order_WT)

mean_RA_order_WT <- mean_RA_order_WT[order(-mean_RA_order_WT$Mean),]

unique(df.long_order_WT$Order)

mean_RA_order_WT %>% top_n(20)

#alternatively you can check for orders with mean RA above 1 percent
above1per <- mean_RA_order_WT[mean_RA_order_WT$Mean > 0.01, ]
unique(above1per$Order)


#these would be 25 different orders, maybe best to keep at 20 as otherwise even more colours needed to distinguish

### replace name taxa not included in top20 to "Other"

library(stringr)
library(dplyr) 

do_not_replace_WT <- c("Burkholderiales","Caulobacterales","Chloroflexales","Flavobacteriales","Frankiales","Gaiellales","Gemmatimonadales","Micrococcales","Micromonosporales","Propionibacteriales","Pseudomonadales","Pseudonocardiales","Rhizobiales","Sphingomonadales","Streptomycetales","Unknown","Xanthomonadales")

df.long_order_WT$Order <- ifelse(df.long_order_WT$Order$Order%in%do_not_replace_WT, df.long_order_WT$Order$Order, "Other")

unique(df.long_order_WT$Order)

# stacked bar plot

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggh4x)

### set order soils

df.long_order_WT$Soil <- factor(df.long_order_WT$Soil, levels = c("NPK","PK","UF"))


### visualization by ggplot2

levels(df.long_order_WT$Soil)<- c("NPK","PK","UF") ### This is to make the letter italic for the genotype name

df.long_order_WT$Soil <- factor(df.long_order_WT$Soil, levels = c("NPK","PK","UF"))

order_levels <- c("Burkholderiales","Caulobacterales","Chloroflexales","Flavobacteriales","Frankiales","Gaiellales","Gemmatimonadales","Micrococcales","Micromonosporales","Propionibacteriales","Pseudomonadales","Pseudonocardiales","Rhizobiales","Sphingomonadales","Streptomycetales","Unknown","Xanthomonadales","Other")
df.long_order_WT$Order <- factor(df.long_order_WT$Order, levels = order_levels)

levels(df.long_order_WT$Compartment)<- c("rhizosphere","root")
df.long_order_WT$Compartment <- factor(df.long_order_WT$Compartment, levels = c("rhizosphere","root"))


#chosen colours
colors <- data.frame(group=c("Burkholderiales","Caulobacterales","Chloroflexales","Flavobacteriales","Frankiales","Gaiellales","Gemmatimonadales","Micrococcales","Micromonosporales","Propionibacteriales","Pseudomonadales","Pseudonocardiales","Rhizobiales","Sphingomonadales","Streptomycetales","Unknown","Xanthomonadales","Other"), 
                     colors=c("#645394","#AA4488","#CC99BB","#ffeeef","#114477","#4477AA","#77AADD","#44AAAA","#77CCCC","#117744","#88CCAA","#CDEBC5","lightyellow","#fdbb6b","#ffd7b5","darkgrey","#ffc0cb","lightgrey")) ### color used here is from "Paired"



colors <- colors[colors$group %in%df.long_order_WT$Order, ]


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

p2 <- ggplot(df.long_order_WT, aes(x=sampleID, y = RA, fill = Order)) +
  geom_bar(stat = "identity", width = .5) +
  facet_nested(~Compartment+Soil,scales ="free_x", labeller = label_parsed)+
  scale_fill_manual(values=colors$colors)+
  main_theme+
  ggtitle("Barley")+
  ylab("Relative abundance")+
  theme(legend.position = "right")+ guides(fill=guide_legend(nrow=21))+ 
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),  # Adjust the size here
        strip.text.x = element_text(size = 25, face = "bold"),
        legend.text=element_text(size=25),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=25))


p2

ggsave(paste("Barley_WT_order_RA.png", sep=""), p2, width=12, height=8)
ggsave(paste("Barley_WT_order_RA.pdf", sep=""), p2, width=12, height=8)
