#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path BarleyCSSP_Askov_reseq_manifest_combined.txt \
  --output-path Barley_Askov_Rep_paired-end-demux.qza \
