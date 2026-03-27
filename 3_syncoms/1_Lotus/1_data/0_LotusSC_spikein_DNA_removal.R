# Clean up
options(warn = -1)
rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the ASV table.
asv_table <- read.table(
  "LotusSC_ASVtable.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)

# Removing spike-in sequences, "85fa8bb918a926d97659d9b64ca6fedd", 
# which is not used in the remaining analysis
asv_table <- asv_table[!rownames(asv_table) %in% "85fa8bb918a926d97659d9b64ca6fedd", ]

# Saving ASV table without spike-in
write.table(
  asv_table,
  file = "LotusSC_ASVtable_nospike.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
