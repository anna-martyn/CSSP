#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime tools import \
  --input-path Askov_Lotus_ASVtable10_4_rarefied_unplanted_removed.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path Askov_Lotus_ASVtable10_4_rarefied_unplanted_removed.qza
