#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=10:0:0
#SBATCH --account InRoot

qiime tools export \
  --input-path LotusSep_exclUFnew_10_4_silva138_taxonomy.qza \
  --output-path LotusSep_exclUFnew_10_4_silva138_taxonomy

