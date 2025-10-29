#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime tools export \
  --input-path LotusSC_new-ref-seqs-or-99.qza \
  --output-path LotusSC_new-ref-seqs-or-99

mv LotusSC_new-ref-seqs-or-99/dna-sequences.fasta LotusSC_new-ref-seqs-or-99/LotusSC_new-ref-seqs-or-99.fasta

qiime tools export \
  --input-path LotusSC_clustered_rep_seqs_99.qza \
  --output-path LotusSC_clustered_rep_seqs_99

mv LotusSC_clustered_rep_seqs_99/dna-sequences.fasta LotusSC_clustered_rep_seqs_99/LotusSC_clustered_rep_seqs_99.fasta

qiime tools export \
  --input-path LotusSC_clustered_table_99.qza \
  --output-path LotusSC_clustered_table_99

biom convert \
  -i LotusSC_clustered_table_99/feature-table.biom \
  -o LotusSC_clustered_table_99/LotusSC_ASVtable.tsv \
  --to-tsv


