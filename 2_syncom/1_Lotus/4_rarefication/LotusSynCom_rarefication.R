# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the ASV table.
asv_table <- read.table(
  "LotusSC_ASVtable_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Load required packages.
library(BiocManager)
library(phyloseq)

# Load ASV table as phyloseq object.
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = T))

# Check minimum sequencing depth in the dataset.
min_depth <- min(sample_sums(phs))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Check overview read counts in all samples.
depths <- sample_sums(phs)
sorted_depths <- sort(depths)

# We will perform the rarefication in two ways:
## 1. Including all ASVs (the ones matched to a Lotus SynCom member, and other ASVs/contaminants).
## 2. Including only ASVs matched to Lotus SynCom members.

# For this we will make a second filtered dataframe, which we also load as phyloseq file.
asv_table_filt <- asv_table[grepl("Lj", rownames(asv_table)), ]
phs_filt <- phyloseq(otu_table(asv_table_filt, taxa_are_rows = T))

min_depth <- min(sample_sums(phs_filt))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Then we will rarefy both ASV tables using phyloseq.
set.seed(1673967505)
RR <- rarefy_even_depth(phs)
RR_filt <- rarefy_even_depth(phs_filt)

# Convert rarefied ASV tables to dataframes.
M <- as.matrix(otu_table(RR))
ASVtable_rarefied <- as.data.frame(M)

M_filt <- as.matrix(otu_table(RR_filt))
ASVtable_rarefied_filt <- as.data.frame(M_filt)

# Add ASVid column from row names.
ASVtable_rarefied$ASVid <- rownames(ASVtable_rarefied)
ASVtable_rarefied <- ASVtable_rarefied[, c("ASVid", setdiff(names(ASVtable_rarefied), "ASVid"))]

ASVtable_rarefied_filt$ASVid <- rownames(ASVtable_rarefied_filt)
ASVtable_rarefied_filt <- ASVtable_rarefied_filt[, c("ASVid", setdiff(names(ASVtable_rarefied_filt), "ASVid"))]

# Save the rarefied ASV tables as csv and txt files.
write.csv(ASVtable_rarefied, file =  "LotusSC_rfd.csv")
write.table(ASVtable_rarefied, 
            file = "LotusSC_rfd.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)

write.csv(ASVtable_rarefied_filt, file =  "LotusSC_rfd_nocontaminants.csv")
write.table(ASVtable_rarefied_filt, 
            file = "LotusSC_rfd_nocontaminants.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)

