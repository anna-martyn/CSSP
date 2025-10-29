#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=15:0:0
#SBATCH --account InRoot


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs LotusCSSP_AskovSoils_paired-end-demux.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 240 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 30 \
  --o-table ./LotusCSSP_AskovSoils_ASVtable.qza \
  --o-representative-sequences ./LotusCSSP_AskovSoils_repseqs.qza \
  --o-denoising-stats ./LotusCSSP_AskovSoils_stats.qza

