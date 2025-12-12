# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the ASV table.
asv_table <- read.table(
  "LotusCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1    
)

# For the Lotus CSSP Askov soil samples, spike-in DNA was added during library preparation.
# As we are not using the spike-in info in our analysis, we want to remove the spike-in ASV from the ASV table for further analysis.
# The id of the spike-in ASV is "85fa8bb918a926d97659d9b64ca6fedd"
asv_table <- asv_table[!rownames(asv_table) %in% "85fa8bb918a926d97659d9b64ca6fedd", ]

# Save the ASV table where the spike-in sequence has been removed.
write.table(
  asv_table,
  file = "LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)