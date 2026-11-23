# Cleaning up
options(warn = -1)
rm(list = ls())

# Load required packages.
library(Biostrings)

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
asv_table <- read.table(
  "LotusSC_ASVtable.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)

# Load representative ASV sequence file.
repseqs <- readDNAStringSet("LotusSC_clustered_rep_seqs_99.fasta")
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

# Here we only find the expected spike-in ASV (85fa8bb918a926d97659d9b64ca6fedd), there is
# no spike-in variant to be found. All other sequences match the spike-in sequence with max. ~50% identity.

# We will therefore only remove the one spike-in ASV from the ASV table.
spike_ids <- c("85fa8bb918a926d97659d9b64ca6fedd")
asv_table <- asv_table[!rownames(asv_table) %in% spike_ids, ]

# Saving ASV table without spike-in.
write.table(
  asv_table,
  file = "LotusSC_ASVtable_nospike.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
