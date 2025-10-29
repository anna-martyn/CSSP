# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the ASV table.
asv_table <- read.table(
  "feature-table_BarleyCSSP_CerealSConly.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1     
)
# Install required packages.
library(BiocManager)
library(phyloseq)

# Filter ASV table to only keep samples of interest.
asv_table <- asv_table[, grepl("input|WT|symrk|ccamk|nsp1|nsp2", colnames(asv_table))]

# Load ASV table as phyloseq object.
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = T))

# Check minimum sequencing depth in the dataset.
min_depth <- min(sample_sums(phs))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Check overview read counts in all samples.
depths <- sample_sums(phs)
sorted_depths <- sort(depths)

# The input samples and one ccamk rhizosphere sample have <1000 reads. We will remove them for rarefication.
sample_depths <- colSums(asv_table)
keep_samples <- names(sample_depths[sample_depths >= 1000])
cat("Removing", ncol(asv_table) - length(keep_samples), "samples with <1000 reads\n")

asv_table_filtered <- asv_table[, keep_samples]

# Now we perform the rarefication in two ways:
## 1. Including all ASVs (the ones matched to a Lotus SynCom member, and other ASVs/contaminants).
## 2. Including only ASVs matched to Lotus SynCom members.

# For this we will make a second filtered dataframe, where contaminants are removed.
asv_table_matched <- asv_table_filtered[grepl("_", rownames(asv_table_filtered)), , drop = FALSE]

# Now we will load both asv tables as phyloseq files again.
## all ASVs
phs_all <- phyloseq(otu_table(asv_table_filtered, taxa_are_rows = TRUE))
min_depth <- min(sample_sums(phs_all))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

## filtered ASVs
phs_matched <- phyloseq(otu_table(asv_table_matched, taxa_are_rows = TRUE))
min_depth <- min(sample_sums(phs_matched))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Then we will rarefy both ASV tables using phyloseq.
set.seed(1673967505)
RR_all <- rarefy_even_depth(phs_all)

RR_matched <- rarefy_even_depth(phs_matched)

# Convert rarefied ASV tables to dataframes.
M <- as.matrix(otu_table(RR_all))
ASVtable_rarefied <- as.data.frame(M)

M_filt <- as.matrix(otu_table(RR_matched))
ASVtable_rarefied_filt <- as.data.frame(M_filt)

# Add ASVid column from row names.
ASVtable_rarefied$ASVid <- rownames(ASVtable_rarefied)
ASVtable_rarefied <- ASVtable_rarefied[, c("ASVid", setdiff(names(ASVtable_rarefied), "ASVid"))]

ASVtable_rarefied_filt$ASVid <- rownames(ASVtable_rarefied_filt)
ASVtable_rarefied_filt <- ASVtable_rarefied_filt[, c("ASVid", setdiff(names(ASVtable_rarefied_filt), "ASVid"))]

# Save the rarefied ASV tables as csv and txt files.
write.csv(ASVtable_rarefied, file =  "HordeumSynCom_rfd_less1000readsremoved.csv")
write.table(ASVtable_rarefied, 
            file = "HordeumSynCom_rfd_less1000readsremoved.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)

write.csv(ASVtable_rarefied_filt, file =  "HordeumSynCom_rfd_less1000readsremoved_nocontaminants.csv")
write.table(ASVtable_rarefied_filt, 
            file = "HordeumSynCom_rfd_less1000readsremoved_nocontaminants.txt", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE, 
            quote = FALSE)

