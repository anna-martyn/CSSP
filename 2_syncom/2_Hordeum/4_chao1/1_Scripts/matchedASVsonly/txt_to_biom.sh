#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

biom convert -i HordeumSynCom_rfd_less1000readsremoved_nocontaminants.txt -o HordeumSynCom_rfd_less1000readsremoved_nocontaminants.biom --table-type="OTU table" --to-hdf5

