#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot

qiime metadata tabulate \
 --m-input-file ./LotusSYM_LjSC_220200e2_stats.qza \
 --o-visualization  ./LotusSYM_LjSC_220200e2_stats.qzv \

qiime feature-table summarize \
  --i-table ./LotusSYM_LjSC_220200e2_ASVtable.qza \
  --m-sample-metadata-file ./LotusSYM_LjSC_metadata.txt \
  --o-visualization ./LotusSYM_LjSC_220200e2_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./LotusSYM_LjSC_220200e2_repseqs.qza \
 --o-visualization ./LotusSYM_LjSC_220200e2_repseqs.qzv
