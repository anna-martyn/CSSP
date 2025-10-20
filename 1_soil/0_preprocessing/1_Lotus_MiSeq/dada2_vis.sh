#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH --time=02:0:0
#SBATCH --account InRoot


qiime metadata tabulate \
 --m-input-file ./LotusSep_exclUFnew_stats.qza \
 --o-visualization  ./LotusSep_exclUFnew_stats.qzv \

qiime feature-table summarize \
  --i-table ./LotusSep_exclUFnew_ASVtable.qza \
  --m-sample-metadata-file ./Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt \
  --o-visualization ./LotusSep_exclUFnew_ASVtable.qzv

qiime feature-table tabulate-seqs \
 --i-data ./LotusSep_exclUFnew_repseqs.qza \
 --o-visualization ./LotusSep_exclUFnew_repseqs.qzv
