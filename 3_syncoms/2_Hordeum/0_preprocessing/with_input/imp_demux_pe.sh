#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path HordeumSC_manifest.txt \
  --output-path HordeumSC_paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 \

qiime demux summarize \
  --i-data ./HordeumSC_paired-end-demux.qza \
  --o-visualization ./HordeumSC_paired-end-demux.qzv

