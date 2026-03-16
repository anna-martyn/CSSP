# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
asv_table <- read.table(
  file = "../../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)

# Load packages
library(phyloseq)

# Create phyloseq object
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = TRUE))

# Filtering samples with >=1000 reads
depths <- sample_sums(phs)
keep_samples <- names(depths[depths >= 1000])
asv_table_filtered <- asv_table[, keep_samples]

# Rarefy filtered ASV table
phs_filt <- phyloseq(otu_table(asv_table_filtered, taxa_are_rows = TRUE))
set.seed(1673967505)
phs_rarefied <- rarefy_even_depth(phs_filt)

# Convert rarefied ASV table to dataframe with column for ASV IDs
asv_table_rarefied <- as.matrix(otu_table(phs_rarefied))
asv_table_rarefied <- data.frame(ASVid = rownames(asv_table_rarefied), asv_table_rarefied)

# Save rarefied ASV table
write.csv(asv_table_rarefied, file = "HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.csv")
write.table(
  asv_table_rarefied, 
  file = "HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.txt", 
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE, 
  quote = FALSE
)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      