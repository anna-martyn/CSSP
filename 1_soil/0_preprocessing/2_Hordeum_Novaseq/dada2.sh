#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=15:0:0
#SBATCH --account InRoot


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./HordeumCSSP_AskovSoils_paired-end-demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --p-max-ee-f 1 \
  --p-max-ee-r 1 \
  --p-n-threads 30 \
  --o-table ./HordeumCSSP_AskovSoils_ASVtable.qza \
  --o-representative-sequences ./HordeumCSSP_AskovSoils_repseqs.qza \
  --o-denoising-stats ./HordeumCSSP_AskovSoils_stats.qza

