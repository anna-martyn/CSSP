#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime vsearch cluster-features-open-reference \
  --i-table LotusSYM_LjSC_220200e2_ASVtable.qza \
  --i-sequences LotusSYM_LjSC_220200e2_repseqs.qza \
  --i-reference-sequences Sha_LotusSCv5v7.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table LotusSYM_LjSC_clustered_table_99.qza \
  --o-clustered-sequences LotusSYM_LjSC_clustered_rep_seqs_99.qza \
  --o-new-reference-sequences LotusSYM_LjSC_new-ref-seqs-or-99.qza
