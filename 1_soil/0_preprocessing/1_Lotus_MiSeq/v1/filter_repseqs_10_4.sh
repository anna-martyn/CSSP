#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=64g
#SBATCH -c 1
#SBATCH --time=1:0:0
#SBATCH --account InRoot


qiime feature-table filter-seqs \
  --i-data ./LotusCSSP_AskovSoils_repseqs.qza \
  --i-table ./LotusCSSP_AskovSoils_ASVtable_10_4.qza \
  --o-filtered-data ./LotusCSSP_AskovSoils_repseqs_10_4.qza

