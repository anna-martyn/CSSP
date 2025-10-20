#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

biom convert -i Askov_ASVtable10_4_rarefied_84.txt -o Askov_ASVtable10_4_rarefied_84.biom --table-type="OTU table" --to-hdf5

