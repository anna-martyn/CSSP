#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime diversity alpha \
  --i-table LotusCSSP_AskovSoils_ASVtable10_4_rfd.qza \
  --p-metric chao1 \
  --output-dir alpha-diversity

qiime tools export \
  --input-path alpha-diversity/alpha_diversity.qza \
  --output-path alpha-diversity/exported

mv alpha-diversity/exported/alpha-diversity.tsv \
   alpha-diversity/LotusCSSP_AskovSoils_chao1.txt

rm -r alpha-diversity/exported
