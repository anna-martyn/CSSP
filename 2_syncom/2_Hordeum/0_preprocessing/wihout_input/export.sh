#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime tools export \
  --input-path HordeumSC_new-ref-seqs-or-99.qza \
  --output-path HordeumSC_new-ref-seqs-or-99

mv HordeumSC_new-ref-seqs-or-99/dna-sequences.fasta HordeumSC_new-ref-seqs-or-99/HordeumSC_new-ref-seqs-or-99.fasta

qiime tools export \
  --input-path HordeumSC_clustered_rep_seqs_99.qza \
  --output-path HordeumSC_clustered_rep_seqs_99

mv HordeumSC_clustered_rep_seqs_99/dna-sequences.fasta HordeumSC_clustered_rep_seqs_99/HordeumSC_clustered_rep_seqs_99.fasta

qiime tools export \
  --input-path HordeumSC_clustered_table_99.qza \
  --output-path HordeumSC_clustered_table_99

biom convert \
  -i HordeumSC_clustered_table_99/feature-table.biom \
  -o HordeumSC_clustered_table_99/HordeumSC_ASVtable.tsv \
  --to-tsv


