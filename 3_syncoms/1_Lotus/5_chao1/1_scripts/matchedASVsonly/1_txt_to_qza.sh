#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

biom convert -i LotusSC_rfd_nocontaminants.txt -o LotusSC_rfd_nocontaminants.biom --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path LotusSC_rfd_nocontaminants.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path LotusSC_rfd_nocontaminants.qza
