#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=1:0:0
#SBATCH --account InRoot

qiime diversity alpha \
  --i-table Askov_Lotus_ASVtable10_4_rarefied_unplanted_removed.qza \
  --p-metric chao1 \
  --output-dir alpha-diversity


