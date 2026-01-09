#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=15:0:0
#SBATCH --account InRoot


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./HordeumSC_paired-end-demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 30 \
  --o-table ./HordeumSC_ASVtable.qza \
  --o-representative-sequences ./HordeumSC_repseqs.qza \
  --o-denoising-stats ./HordeumSC_stats.qza
