# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
asv_table <- read.table(
  "LotusCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1    
)

# Removing spike-in
asv_table <- asv_table[!rownames(asv_table) %in% "85fa8bb918a926d97659d9b64ca6fedd", ]
# The ASV "85fa8bb918a926d97659d9b64ca6fedd" represents spike-in DNA, which was added during
# library preparation, but will not be used in the data analysis

# Saving ASV table without spike-in
write.table(
  asv_table,
  file = "LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
