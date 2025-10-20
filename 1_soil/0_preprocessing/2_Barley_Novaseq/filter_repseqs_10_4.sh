#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=64g
#SBATCH -c 1
#SBATCH --time=1:0:0
#SBATCH --account InRoot


qiime feature-table filter-seqs \
  --i-data ./Barley_Askov_Rep_repseqs.qza \
  --i-table ./Barley_Askov_Rep_ASVtable_10_4.qza \
  --o-filtered-data ./Barley_Askov_Rep_repseqs_10_4.qza

