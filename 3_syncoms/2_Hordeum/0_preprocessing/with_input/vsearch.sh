#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime vsearch cluster-features-open-reference \
  --i-table HordeumSC_ASVtable.qza \
  --i-sequences HordeumSC_repseqs.qza \
  --i-reference-sequences CerealSynCom_v5v7_noDuplicates.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table HordeumSC_clustered_table_99.qza \
  --o-clustered-sequences HordeumSC_clustered_rep_seqs_99.qza \
  --o-new-reference-sequences HordeumSC_new-ref-seqs-or-99.qza
