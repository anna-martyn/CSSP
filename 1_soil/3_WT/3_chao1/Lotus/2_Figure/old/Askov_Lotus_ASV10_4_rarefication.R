options(warn=-1)

# cleanup
rm(list=ls())

#
setwd("O:/Nat_MBG-PMg/Anna_Martyn/0_MARTYN_THORSGAARD/2_Soil_experiment/2_WT/Anna analysis separate/3_chao1/Lotus/rarefied/")
results.dir <- "O:/Nat_MBG-PMg/Anna_Martyn/0_MARTYN_THORSGAARD/2_Soil_experiment/2_WT/Anna analysis separate/3_chao1/Lotus/rarefied/"
figures.dir <- "O:/Nat_MBG-PMg/Anna_Martyn/0_MARTYN_THORSGAARD/2_Soil_experiment/2_WT/Anna analysis separate/3_chao1/Lotus/rarefied/"

otu_table.file <- paste(results.dir, "LotusSep_exclUFnew_ASVtable_10_4_nospike.txt", sep="")

# load data
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)

### rarefication
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")

library(BiocManager)

# BiocManager::install("phyloseq")
# install.packages("phyloseq")
library(phyloseq)

phs <- phyloseq(otu_table(otu_table, taxa_are_rows = T))
set.seed(1673967505)
RR <- rarefy_even_depth(phs)
M <- as.matrix(otu_table(RR))
ASVtable_rarefied <- as.data.frame(M)

write.csv(M, file =  "Askov_Lotus_ASVtable10_4_rarefied.csv")

