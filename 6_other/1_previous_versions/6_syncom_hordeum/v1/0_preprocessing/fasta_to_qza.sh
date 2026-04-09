#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime tools import \
  --input-path Cereal_SynCom_v5v7_exclisolates_noDup.fasta \
  --output-path Cereal_SynCom_v5v7_exclisolates_noDup.qza \
  --type 'FeatureData[Sequence]'
