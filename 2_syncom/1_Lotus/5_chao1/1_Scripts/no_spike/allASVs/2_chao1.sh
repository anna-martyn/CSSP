#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime diversity alpha \
  --i-table LotusSC_rfd.qza \
  --p-metric chao1 \
  --output-dir alpha-diversity_allASVs

qiime tools export \
  --input-path alpha-diversity_allASVs/alpha_diversity.qza \
  --output-path alpha-diversity_allASVs/exported

mv alpha-diversity_allASVs/exported/alpha-diversity.tsv \
   alpha-diversity_allASVs/LotusSC_allASVs_chao1.txt

rm -r alpha-diversity_allASVs/exported
