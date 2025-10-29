#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

biom convert -i HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.txt -o HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.biom --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path HordeumCSSP_AskovSoils_ASVtable_10_4_rfd_min1000reads.qza
