# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# For Hordeum rarefication was already performed in the /1_soil_exp/2_bulk/1_rarefication folder. Also chao1 was performed already.

# We will do the same for Lotus now.

# Load the ASV table.
asv_table <- read.table(
  "LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
)

# Install required packages.
library(BiocManager)
library(phyloseq)

# Load ASV table as phyloseq object.
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = T))

# Check minimum sequencing depth in the dataset.
min_depth <- min(sample_sums(phs))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Check overall depth.
depths <- sample_sums(phs)
sorted_depths <- sort(depths)
sorted_depths

# Now perform rarefication.
set.seed(1673967505)
RR <- rarefy_even_depth(phs)

# Convert rarefied ASV table to a dataframe.
M <- as.matrix(otu_table(RR))
ASVtable_rarefied <- as.data.frame(M)

# Add ASVid column from row names.
ASVtable_rarefied$ASVid <- rownames(ASVtable_rarefied)
ASVtable_rarefied <- ASVtable_rarefied[, c("ASVid", setdiff(names(ASVtable_rarefied), "ASVid"))]

# Save the rarefied ASV table as csv and txt files.
write.csv(ASVtable_rarefied, file =  "LotusCSSP_AskovSoils_ASVtable10_4_rfd.csv")
write.table(ASVtable_rarefied, 
            file = "LotusCSSP_AskovSoils_ASVtable10_4_rfd.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)
