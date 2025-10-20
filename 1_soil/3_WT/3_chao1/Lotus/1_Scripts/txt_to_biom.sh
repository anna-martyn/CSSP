#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=1:0:0
#SBATCH --account InRoot

biom convert -i Askov_Lotus_ASVtable10_4_rarefied_unplanted_removed.txt -o Askov_Lotus_ASVtable10_4_rarefied_unplanted_removed.biom --table-type="OTU table" --to-hdf5

