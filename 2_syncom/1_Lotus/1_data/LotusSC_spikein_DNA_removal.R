# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the ASV table.
asv_table <- read.table(
  "LotusSC_ASVtable.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",  # IMPORTANT
  skip = 1            # skip the first metadata line
)

# For the Lotus SynCom samples, spike-in DNA was added during library preparation.
# As we are not using the spike-in info in our analysis, we want to remove the spike-in ASV from the ASV table for further aanalysis.
# The id of the spike-in ASV is "85fa8bb918a926d97659d9b64ca6fedd"
asv_table <- asv_table[!rownames(asv_table) %in% "85fa8bb918a926d97659d9b64ca6fedd", ]

# Save the ASV table where the spike-in sequence has been removed.
write.table(
  asv_table,
  file = "LotusSC_ASVtable_nospike.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
