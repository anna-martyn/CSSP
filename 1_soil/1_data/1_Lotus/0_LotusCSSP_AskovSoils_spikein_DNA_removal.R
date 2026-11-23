# Cleaning up
options(warn = -1)
rm(list = ls())

# Install required packages.
library(Biostrings)

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

# Load representative ASV sequence file.
repseqs <- readDNAStringSet("LotusCSSP_AskovSoils_repseqs_10_4.fasta")
names(repseqs) <- sub(" .*", "", names(repseqs))

# We added spike-in DNA to the samples during library preparation, but we will not use it during data analysis.
# Therefore, we want to remove the spike ASV from the ASV table.
# Before doing so, we, however, want to check whether there are any additional spike variants present.

# Define the spike-in ASV.
spike_id  <- "85fa8bb918a926d97659d9b64ca6fedd"
spike_seq <- repseqs[[spike_id]]

# Check sequence similarity for all other sequences.
aln <- pairwiseAlignment(repseqs, spike_seq, type = "global")

spike_similarity <- data.frame(
  ASV          = names(repseqs),
  pct_identity = round(pid(aln), 2),
  edit_dist    = adist(as.character(spike_seq), as.character(repseqs))[1, ],
  length       = width(repseqs),
  total_reads  = rowSums(asv_table)[names(repseqs)]
)

# Sort output by similarity % and inspect data.
spike_similarity <- spike_similarity[order(-spike_similarity$pct_identity,
                                           spike_similarity$edit_dist), ]
rownames(spike_similarity) <- NULL

print (spike_similarity)

# Besides the expected spike-in ASV (85fa8bb918a926d97659d9b64ca6fedd), there is
# a second ASV (25f50286877b75132ef5ec0dc684eb73) with 99,74% identity to the spike-in sequence
# (1bp difference, same length of 382bp), which is also not classified as bacterial in the taxonomy file.
# All other sequences share only ~50% identity or less. Hence, this ASV is a spike variant.

# We will therefore make a new ASV table that excludes both variants.
spike_ids <- c("85fa8bb918a926d97659d9b64ca6fedd",
               "25f50286877b75132ef5ec0dc684eb73")
asv_table <- asv_table[!rownames(asv_table) %in% spike_ids, ]

# Saving ASV table without spike-in.
write.table(
  asv_table,
  file = "LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
