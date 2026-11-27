# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
asv_table <- read.table(
  "../1_data/HordeumSC_ASVtable.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1     
)

# Loading packages
library(phyloseq)

# Filtering low-depth samples -------------------------------------------------
# Checking minimum sequencing depth ASV table
min_depth <- min(colSums(asv_table))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Examining sequencing depths
sample_depths <- colSums(asv_table)
sort(sample_depths)[1:5]

# Filtering out sample Hv_ccamk_HvSC_only_rhizo_1 with 917 reads
idx <- which(colnames(asv_table) == "Hv_ccamk_HvSC_only_rhizo_1")
asv_table_filtered <- asv_table[, -idx]

# Filtered ASV table including only ASVs matched to SynCom sequences
asv_table_matched <- asv_table_filtered[
  grepl("_", rownames(asv_table_filtered)),
  ,
  drop = FALSE
]

# Checking the minimum number of reads after filtering 
# all ASVs
phs_all <- phyloseq(otu_table(asv_table_filtered, taxa_are_rows = TRUE))
min_depth <- min(sample_sums(phs_all))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

## filtered ASVs
phs_matched <- phyloseq(otu_table(asv_table_matched, taxa_are_rows = TRUE))
min_depth <- min(sample_sums(phs_matched))
cat("Lowest sequencing depth in the dataset:", min_depth, "\n")

# Rarefying -------------------------------------------------------------------
# Rarefying both full ASV table and matched ASVs only
set.seed(1673967505)
phs_rare <- rarefy_even_depth(phs_all)
phs_rare_filt <- rarefy_even_depth(phs_matched)

# Converting rarefied ASV tables to dataframes
asv_table_rare <- as.matrix(otu_table(phs_rare))
asv_table_rare <- as.data.frame(asv_table_rare)

asv_table_rare_filt <- as.matrix(otu_table(phs_rare_filt))
asv_table_rare_filt <- as.data.frame(asv_table_rare_filt)

# Adding 'ASVid' as a column
asv_table_rare$ASVid <- rownames(asv_table_rare)
asv_table_rare <- asv_table_rare[, c(
  "ASVid",
  setdiff(names(asv_table_rare), "ASVid")
)]

asv_table_rare_filt$ASVid <- rownames(asv_table_rare_filt)
asv_table_rare_filt <- asv_table_rare_filt[, c(
  "ASVid",
  setdiff(names(asv_table_rare_filt), "ASVid")
)]

# Saving rarefied ASV tables
write.csv(
  asv_table_rare,
  file = "1_tables/HordeumSC_rfd.csv"
)
write.table(
  asv_table_rare,
  file = "1_tables/HordeumSC_rfd.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.csv(
  asv_table_rare_filt,
  file = "1_tables/HordeumSC_rfd_nocontaminants.csv"
)
write.table(
  asv_table_rare_filt,
  file = "1_tables/HordeumSC_rfd_nocontaminants.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

