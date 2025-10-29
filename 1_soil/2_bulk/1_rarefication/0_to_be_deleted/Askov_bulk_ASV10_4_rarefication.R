# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the ASV table.
asv_table.file <- paste("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="")
asv_table <- read.table(asv_table.file, sep="\t", header=T, row.names=1, check.names=F)

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

# We will rarefy in two ways, once using the lowest depth as threshold (86), and once the second lowest (~6603).

## 1) Rarefy based on lowest value.
# Rarefy the ASV table using phyloseq.
set.seed(1673967505)
RR <- rarefy_even_depth(phs)

# Convert rarefied ASV table to a dataframe.
M <- as.matrix(otu_table(RR))
ASVtable_rarefied <- as.data.frame(M)

# Add ASVid column from row names.
ASVtable_rarefied$ASVid <- rownames(ASVtable_rarefied)
ASVtable_rarefied <- ASVtable_rarefied[, c("ASVid", setdiff(names(ASVtable_rarefied), "ASVid"))]

# Save the rarefied ASV table as csv and txt files.
write.csv(ASVtable_rarefied, file =  "Askov_ASVtable10_4_rarefied_84.csv")
write.table(ASVtable_rarefied, 
            file = "Askov_ASVtable10_4_rarefied_84.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)

## 2) Rarefy based on the second lowest value.
# Rarefy the ASV table using phyloseq.
set.seed(1673967505)
RR2 <- rarefy_even_depth(phs, sample.size = second_smallest)

# Comment: sample AM_S_3 was now removed due to low reads (86).

# Convert rarefied ASV table to a dataframe.
M2 <- as.matrix(otu_table(RR2))
ASVtable_rarefied2 <- as.data.frame(M2)

# Add ASVid column from row names.
ASVtable_rarefied2$ASVid <- rownames(ASVtable_rarefied2)
ASVtable_rarefied2 <- ASVtable_rarefied2[, c("ASVid", setdiff(names(ASVtable_rarefied2), "ASVid"))]

# Save the rarefied ASV table as csv and txt files.
write.csv(ASVtable_rarefied2, file = "Askov_ASVtable10_4_rarefied_6603.csv")
write.table(ASVtable_rarefied2, 
            file = "Askov_ASVtable10_4_rarefied_6603.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      