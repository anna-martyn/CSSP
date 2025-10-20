#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools export \
  --input-path Barley_Askov_Rep_silva138_taxonomy.qza \
  --output-path Barley_Askov_Rep_silva138_taxonomy

