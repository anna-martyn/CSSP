#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

biom convert -i LotusSynCom_rfd_nounplanted_nocontaminants.txt -o LotusSynCom_rfd_nounplanted_nocontaminants.biom --table-type="OTU table" --to-hdf5

