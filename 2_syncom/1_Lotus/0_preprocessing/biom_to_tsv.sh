#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

biom convert -i feature-table.biom -o feature-table_LotusSYM_LjSC.tsv --to-tsv
