#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime vsearch cluster-features-open-reference \
  --i-table LotusSC_ASVtable.qza \
  --i-sequences LotusSC_repseqs.qza \
  --i-reference-sequences LotusSCv5v7.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table LotusSC_clustered_table_99.qza \
  --o-clustered-sequences LotusSC_clustered_rep_seqs_99.qza \
  --o-new-reference-sequences LotusSC_new-ref-seqs-or-99.qza
