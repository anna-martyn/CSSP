#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path LotusCSSP_AskovSoils_manifest.txt \
  --output-path LotusCSSP_AskovSoils_paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 \

qiime demux summarize \
  --i-data ./LotusCSSP_AskovSoils_paired-end-demux.qza \
  --o-visualization ./LotusCSSP_AskovSoils_paired-end-demux.qzv

