#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot


qiime metadata tabulate \
 --m-input-file ./HordeumCSSP_AskovSoils_stats.qza \
 --o-visualization  ./HordeumCSSP_AskovSoils_stats.qzv \

qiime feature-table summarize \
  --i-table ./HordeumCSSP_AskovSoils_ASVtable.qza \
  --m-sample-metadata-file ./HordeumCSSP_AskovSoils_metadata.txt \
  --o-visualization ./HordeumCSSP_AskovSoils_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./HordeumCSSP_AskovSoils_repseqs.qza \
 --o-visualization ./HordeumCSSP_AskovSoils_repseqs.qzv
