#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

biom convert -i HordeumSC_rfd_woinput_min1000reads_nocontaminants.txt -o HordeumSC_rfd_woinput_min1000reads_nocontaminants.biom --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path HordeumSC_rfd_woinput_min1000reads_nocontaminants.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path HordeumSC_rfd_woinput_min1000reads_nocontaminants.qza
