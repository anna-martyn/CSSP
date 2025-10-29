# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the ASV table.
asv_table <- read.table("HordeumCSSP_AskovSoils_ASVtable_10_4.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char = "", skip = 1)

# Install required packages.
library(BiocManager)
BiocManager::install("phyloseq")
library(phyloseq)

# Load ASV table as phyloseq object.
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = T))

# Check minimum sequencing depth in the dataset.
min_depth <- min(sample_sums(phs))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Check second lowest, as the lowest one is only 86 reads.
depths <- sample_sums(phs)
sorted_depths <- sort(depths)
second_smallest <- sorted_depths[2]
cat("Second smallest sequencing depth:", second_smallest, "\n")

# We will rarefy removing the low-depth sample.

## Remove low-depth sample.
keep_samples <- names(depths[depths >= 1000])
cat("Removing", ncol(asv_table) - length(keep_samples), "samples with <1000 reads\n")

asv_table_filtered <- asv_table[, keep_samples]

## Rarefy the ASV table using phyloseq.
phs_filt <- phyloseq(otu_table(asv_table_filtered, taxa_are_rows = TRUE))

set.seed(1673967505)
RR_filt <- rarefy_even_depth(phs_filt)

## Convert rarefied ASV table to a dataframe.
M <- as.matrix(otu_table(RR_filt))
ASVtable_rarefied <- as.data.frame(M)

## Add ASVid column from row names.
ASVtable_rarefied$ASVid <- rownames(ASVtable_rarefied)
ASVtable_rarefied <- ASVtable_rarefied[, c("ASVid", setdiff(names(ASVtable_rarefied), "ASVid"))]

# Save the rarefied ASV table as csv and txt files.
write.csv(ASVtable_rarefied, file = "HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.csv")
write.table(ASVtable_rarefied, 
            file = "HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      