#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime tools import \
  --input-path HordeumSynCom_rfd_less1000readsremoved_nocontaminants.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path HordeumSynCom_rfd_less1000readsremoved_nocontaminants.qza
