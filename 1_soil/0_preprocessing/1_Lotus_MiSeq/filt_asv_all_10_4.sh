#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=01:0:0
#SBATCH --account InRoot

  
qiime feature-table filter-features \
  --i-table ./LotusSep_exclUFnew_ASVtable.qza \
  --p-min-frequency 10 \
  --p-min-samples 4 \
  --o-filtered-table ./LotusSep_exclUFnew_ASVtable_10_4.qza
  
