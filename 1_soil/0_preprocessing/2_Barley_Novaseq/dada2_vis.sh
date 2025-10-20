#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot


qiime metadata tabulate \
 --m-input-file ./Barley_Askov_Rep_stats.qza \
 --o-visualization  ./Barley_Askov_Rep_stats.qzv \

qiime feature-table summarize \
  --i-table ./Barley_Askov_Rep_ASVtable.qza \
  --m-sample-metadata-file ./BarleyCSSP_Askov_reseq_metadata.txt \
  --o-visualization ./Barley_Askov_Rep_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./Barley_Askov_Rep_repseqs.qza \
 --o-visualization ./Barley_Askov_Rep_repseqs.qzv
