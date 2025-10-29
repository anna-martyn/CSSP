#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot


qiime metadata tabulate \
 --m-input-file ./LotusCSSP_AskovSoils_stats.qza \
 --o-visualization  ./LotusCSSP_AskovSoils_stats.qzv \

qiime feature-table summarize \
  --i-table ./LotusCSSP_AskovSoils_ASVtable.qza \
  --m-sample-metadata-file ./LotusCSSP_AskovSoils_metadata.txt \
  --o-visualization ./LotusCSSP_AskovSoils_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./LotusCSSP_AskovSoils_repseqs.qza \
 --o-visualization ./LotusCSSP_AskovSoils_repseqs.qzv
