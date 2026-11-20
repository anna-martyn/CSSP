#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime diversity alpha \
  --i-table LotusSC_rfd_nocontaminants.qza \
  --p-metric chao1 \
  --output-dir alpha-diversity_matchedASVsonly

qiime tools export \
  --input-path alpha-diversity_matchedASVsonly/alpha_diversity.qza \
  --output-path alpha-diversity_matchedASVsonly/exported

mv alpha-diversity_matchedASVsonly/exported/alpha-diversity.tsv \
   alpha-diversity_matchedASVsonly/LotusSC_matchedASVsonly_chao1.txt

rm -r alpha-diversity_matchedASVsonly/exported
