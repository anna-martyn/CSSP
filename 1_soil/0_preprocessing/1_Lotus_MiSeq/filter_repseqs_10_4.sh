#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=64g
#SBATCH -c 1
#SBATCH --time=1:0:0
#SBATCH --account InRoot


qiime feature-table filter-seqs \
  --i-data ./LotusSep_exclUFnew_repseqs.qza \
  --i-table ./LotusSep_exclUFnew_ASVtable_10_4.qza \
  --o-filtered-data ./LotusSep_exclUFnew_repseqs_10_4.qza

