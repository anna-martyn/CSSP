#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime diversity alpha \
  --i-table Askov_ASVtable10_4_rarefied_6603.qza \
  --p-metric chao1 \
  --output-dir alpha-diversity_6603


