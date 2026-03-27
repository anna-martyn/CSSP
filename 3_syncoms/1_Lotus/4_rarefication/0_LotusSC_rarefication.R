# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
library(phyloseq)

# Loading data
asv_table <- read.table(
  "../1_data/LotusSC_ASVtable_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Creating phyloseq object
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = TRUE))

# Checking minimum sequencing depth ASV table
min_depth <- min(sample_sums(phs))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")
# All samples have a depth of >5000 reads, so no filtering required


# We will perform the rarefication in two ways:
## 1. Including all ASVs (the ones matched to a Lotus SynCom member, and other ASVs/contaminants).
## 2. Including only ASVs matched to Lotus SynCom members. (later used for final plots)

# For this we will make a second filtered dataframe, which we also load as phyloseq file.

# Filtered dataset containing only ASVs matched to SynCom
asv_table_filt <- asv_table[grepl("Lj", rownames(asv_table)), ]
phs_filt <- phyloseq(otu_table(asv_table_filt, taxa_are_rows = TRUE))

# Ensruing that minimum sequencing depth is still sufficuent high
min_depth <- min(sample_sums(phs_filt))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")
# Minimum sequencing depth is still >5000, so no filtering required

# Rarefying full and filtered ASV tables
set.seed(1673967505)
phs_rare <- rarefy_even_depth(phs)
phs_rare_filt <- rarefy_even_depth(phs_filt)

# Converting rarefied ASV tables to dataframes
asv_tabe_rare <- as.matrix(otu_table(phs_rare))
asv_tabe_rare <- as.data.frame(asv_tabe_rare)

asv_tabe_rare_filt <- as.matrix(otu_table(phs_rare_filt))
asv_tabe_rare_filt <- as.data.frame(asv_tabe_rare_filt)

# Adding 'ASVid' as a column
asv_tabe_rare$ASVid <- rownames(asv_tabe_rare)
asv_tabe_rare <- asv_tabe_rare[, c(
  "ASVid",
  setdiff(names(asv_tabe_rare), "ASVid")
)]

asv_tabe_rare_filt$ASVid <- rownames(asv_tabe_rare_filt)
asv_tabe_rare_filt <- asv_tabe_rare_filt[, c(
  "ASVid",
  setdiff(names(asv_tabe_rare_filt), "ASVid")
)]

# Saveing rarefied ASV tables
write.csv(asv_tabe_rare, file = "LotusSC_rfd.csv")
write.table(
  asv_tabe_rare,
  file = "LotusSC_rfd.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.csv(asv_tabe_rare_filt, file = "LotusSC_rfd_nocontaminants.csv")
write.table(
  asv_tabe_rare_filt,
  file = "LotusSC_rfd_nocontaminants.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

