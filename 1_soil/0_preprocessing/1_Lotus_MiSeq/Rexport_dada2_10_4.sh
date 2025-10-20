#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 1
#SBATCH --time=01:0:0
#SBATCH --account InRoot

for i in *ASVtable_10_4.qza; do
qiime tools export --input-path $i --output-path .
done
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
