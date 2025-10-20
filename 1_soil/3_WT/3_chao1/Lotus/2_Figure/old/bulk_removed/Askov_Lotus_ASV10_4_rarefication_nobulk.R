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

#remove bulk samples as not needed for bulk part and as one sample only has 10 reads or so and will make the rarefication impossible

library(dplyr)

#check how many bulk samples there are
bulk_samples <- otu_table %>% select(contains('unplanted'))

#remove them and make new otu table
otu_table_nobulk <- otu_table %>% select(-contains('unplanted'))

### rarefication
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")

library(BiocManager)

# BiocManager::install("phyloseq")
# install.packages("phyloseq")
library(phyloseq)

phs <- phyloseq(otu_table(otu_table_nobulk, taxa_are_rows = T))
set.seed(1673967505)
RR <- rarefy_even_depth(phs)
M <- as.matrix(otu_table(RR))
ASVtable_rarefied <- as.data.frame(M)

write.csv(M, file =  "Askov_Lotus_ASVtable10_4_rarefied_nobulk.csv")

