# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
asv_table <- read.table(
  file = "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
)

# Loading packages
library(phyloseq)

# Phyloseq object
phs <- phyloseq(otu_table(asv_table, taxa_are_rows = TRUE))

# Checking minimum sequencing depth
min_depth <- min(sample_sums(phs))
min_depth > 1000 # >1000 reads in all samples, no further filtering

# Rarefying ASV table
set.seed(1673967505)
phs_rarefied <- rarefy_even_depth(phs)

# Converting rarefied ASV table to dataframe with column for ASV IDs
asv_table_rarefied <- as.matrix(otu_table(phs_rarefied))
asv_table_rarefied <- data.frame(ASVid = rownames(asv_table_rarefied), asv_table_rarefied)

# Saving rarefied ASV tables
write.csv(asv_table_rarefied, file =  "LotusCSSP_AskovSoils_ASVtable10_4_rfd.csv")
write.table(
  asv_table_rarefied,
  file = "LotusCSSP_AskovSoils_ASVtable10_4_rfd.txt", 
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE, 
  quote = FALSE
)
