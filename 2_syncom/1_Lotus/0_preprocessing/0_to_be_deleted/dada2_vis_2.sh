#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot

qiime metadata tabulate \
 --m-input-file ./LotusSC_stats.qza \
 --o-visualization  ./LotusSC_stats.qzv \

qiime feature-table summarize \
  --i-table ./LotusSC_ASVtable.qza \
  --m-sample-metadata-file ./LotusSC_metadata.txt \
  --o-visualization ./LotusSC_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./LotusSC_repseqs.qza \
 --o-visualization ./LotusSC_repseqs.qzv
