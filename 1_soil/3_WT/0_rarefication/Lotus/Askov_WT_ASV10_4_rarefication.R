# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# For Hordeum rarefication was already performed in the /1_soil_exp/2_bulk/1_rarefication folder. Also chao1 was performed already.

# We will do the same for Lotus now.

# Load the ASV table.
asv_table <- read.table(
  "feature-table.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",  # IMPORTANT
  skip = 1            # skip the first metadata line
)

# Install required packages.
library(BiocManager)
library(phyloseq)

# Load ASV table as phyloseq object.
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = T))

# Check minimum sequencing depth in the dataset.
min_depth <- min(sample_sums(phs))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Check second and third lowest, as the lowest one is only 2 reads.
depths <- sample_sums(phs)
sorted_depths <- sort(depths)
sorted_depths

# The two unplanted samples have a very low depth (4 and 13 reads), therefore I will remove the unplanted samples (total of 3) before rarefication, as they are not used in this study anyway.
## Check which samples are unplanted samples.
cols <- grep("unplanted", colnames(asv_table), value = TRUE)
cols

## Remove samples from ASV table.
cols_to_remove <- grep("unplanted", colnames(asv_table))
asv_table_filtered <- asv_table[, -cols_to_remove, drop = FALSE]

## Re-run phyloseq with filtered ASV table.
# Load ASV table as phyloseq object.
phs_filt <- phyloseq(otu_table(asv_table_filtered, taxa_are_rows = T))

## Check minimum sequencing depth in the dataset.
min_depth_filt <- min(sample_sums(phs_filt))
cat("Lowest sequencing depth in the dataset:", min_depth_filt, "\n")

# Now perform rarefication.
set.seed(1673967505)
RR <- rarefy_even_depth(phs_filt)

# Convert rarefied ASV table to a dataframe.
M <- as.matrix(otu_table(RR))
ASVtable_rarefied <- as.data.frame(M)

# Add ASVid column from row names.
ASVtable_rarefied$ASVid <- rownames(ASVtable_rarefied)
ASVtable_rarefied <- ASVtable_rarefied[, c("ASVid", setdiff(names(ASVtable_rarefied), "ASVid"))]

# Save the rarefied ASV table as csv and txt files.
write.csv(ASVtable_rarefied, file =  "Askov_Lotus_ASVtable10_4_rarefied_unplanted_removed.csv")
write.table(ASVtable_rarefied, 
            file = "Askov_Lotus_ASVtable10_4_rarefied_unplanted_removed.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)
