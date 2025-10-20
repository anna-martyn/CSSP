#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=15:0:0
#SBATCH --account InRoot


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./LotusSYM_LjSC_paired-end-demux.qza \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 200 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 30 \
  --o-table ./LotusSYM_LjSC_220200e2_ASVtable.qza \
  --o-representative-sequences ./LotusSYM_LjSC_220200e2_repseqs.qza \
  --o-denoising-stats ./LotusSYM_LjSC_220200e2_stats.qza
