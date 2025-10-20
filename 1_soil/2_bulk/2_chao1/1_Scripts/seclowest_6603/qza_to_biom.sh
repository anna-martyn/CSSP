#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools export \
  --input-path alpha_diversity.qza \
  --output-path alpha_diversity_6603

