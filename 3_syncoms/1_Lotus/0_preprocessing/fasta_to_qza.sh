#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime tools import \
  --input-path LotusSCv5v7.fasta \
  --output-path LotusSCv5v7.qza \
  --type 'FeatureData[Sequence]'
