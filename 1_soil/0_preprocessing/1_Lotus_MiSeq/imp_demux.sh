#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Lotus_CSSP_AskovSoils_manifest_excl_new_bulkUF.txt \
  --output-path LotusSep_exclUFnew_paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 \

qiime demux summarize \
  --i-data ./LotusSep_exclUFnew_paired-end-demux.qza \
  --o-visualization ./LotusSep_exclUFnew_paired-end-demux.qzv

