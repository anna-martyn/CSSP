#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=01:0:0
#SBATCH --account InRoot

  
qiime feature-table filter-features \
  --i-table ./LotusCSSP_AskovSoils_ASVtable.qza \
  --p-min-frequency 10 \
  --p-min-samples 4 \
  --o-filtered-table ./LotusCSSP_AskovSoils_ASVtable_10_4.qza
  
