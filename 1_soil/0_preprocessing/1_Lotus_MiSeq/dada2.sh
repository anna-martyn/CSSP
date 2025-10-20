#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=15:0:0
#SBATCH --account InRoot


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs LotusSep_exclUFnew_paired-end-demux.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 240 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 30 \
  --o-table ./LotusSep_exclUFnew_ASVtable.qza \
  --o-representative-sequences ./LotusSep_exclUFnew_repseqs.qza \
  --o-denoising-stats ./LotusSep_exclUFnew_stats.qza

