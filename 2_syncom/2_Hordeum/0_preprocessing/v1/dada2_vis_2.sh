#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot

qiime metadata tabulate \
 --m-input-file ./BarleyCSSP_SConly_stats.qza  \
 --o-visualization  ./BarleyCSSP_SConly_stats.qzv \

qiime feature-table summarize \
  --i-table ./BarleyCSSP_SConly_ASVtable.qza \
  --m-sample-metadata-file ./BarleyCSSP_SConly_metadata_NEW.txt \
  --o-visualization ./BarleyCSSP_SConly_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./BarleyCSSP_SConly_repseqs.qza \
 --o-visualization ./BarleyCSSP_SConly_repseqs.qzv
