#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime vsearch cluster-features-open-reference \
  --i-table BarleyCSSP_SConly_ASVtable.qza \
  --i-sequences BarleyCSSP_SConly_repseqs.qza \
  --i-reference-sequences Cereal_SynCom_v5v7_exclisolates_noDup.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table BarleyCSSP_CerealSConly_clustered_table_99.qza \
  --o-clustered-sequences BarleyCSSP_CerealSConly_clustered_rep_seqs_99.qza \
  --o-new-reference-sequences BarleyCSSP_CerealSConly_new-ref-seqs-or-99.qza
