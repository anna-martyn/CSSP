#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools export \
  --input-path BarleyCSSP_CerealSConly_new-ref-seqs-or-99.qza \
  --output-path BarleyCSSP_CerealSConly_new-ref-seqs-or-99

