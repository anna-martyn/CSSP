#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path HordeumCSSP_AskovSoils_manifest.txt \
  --output-path HordeumCSSP_AskovSoils_paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 \

qiime demux summarize \
  --i-data ./HordeumCSSP_AskovSoils_paired-end-demux.qza \
  --o-visualization ./HordeumCSSP_AskovSoils_paired-end-demux.qzv

