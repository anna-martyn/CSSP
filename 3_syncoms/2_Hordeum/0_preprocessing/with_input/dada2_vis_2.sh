#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot

qiime metadata tabulate \
 --m-input-file ./HordeumSC_stats.qza  \
 --o-visualization  ./HordeumSC_stats.qzv \

qiime feature-table summarize \
  --i-table ./HordeumSC_ASVtable.qza \
  --m-sample-metadata-file ./HordeumSC_metadata_NEW.txt \
  --o-visualization ./HordeumSC_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./HordeumSC_repseqs.qza \
 --o-visualization ./HordeumSC_repseqs.qzv
