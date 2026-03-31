#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

qiime diversity alpha \
  --i-table HordeumSynCom_rfd_less1000readsremoved_nocontaminants.qza \
  --p-metric chao1 \
  --output-dir alpha-diversity_matchedASVs


