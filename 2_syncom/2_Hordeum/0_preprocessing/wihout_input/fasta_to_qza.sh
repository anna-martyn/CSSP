#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime tools import \
  --input-path CerealSynCom_v5v7_noDuplicates.fasta \
  --output-path CerealSynCom_v5v7_noDuplicates.qza \
  --type 'FeatureData[Sequence]'
