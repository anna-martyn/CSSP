#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=12:0:00
#SBATCH --account InRoot

qiime feature-classifier classify-sklearn \
 --i-reads ./Barley_Askov_Rep_repseqs_10_4.qza \
 --i-classifier ./silva-138-ssu-nr99-v5-v7-classifier_May23.qza \
 --o-classification ./Barley_Askov_Rep_silva138_taxonomy.qza \
 
qiime metadata tabulate \
 --m-input-file ./Barley_Askov_Rep_silva138_taxonomy.qza \
 --o-visualization ./Barley_Askov_Rep_silva138_taxonomy.qzv

