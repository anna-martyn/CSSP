#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=100g
#SBATCH -c 10
#SBATCH --time=2:0:0
#SBATCH --account InRoot

# Export taxonomy
qiime tools export \
  --input-path LotusCSSP_AskovSoils_10_4_silva138_taxonomy.qza \
  --output-path LotusCSSP_AskovSoils_10_4_silva138_taxonomy

mv LotusCSSP_AskovSoils_10_4_silva138_taxonomy/taxonomy.tsv \
   LotusCSSP_AskovSoils_10_4_silva138_taxonomy/LotusCSSP_AskovSoils_taxonomy_10_4.tsv

# Export ASV tables and convert to TSV
for table_qza in *ASVtable_10_4.qza; do
    folder_name="${table_qza%.qza}_export"
    mkdir -p "$folder_name"

    qiime tools export \
        --input-path "$table_qza" \
        --output-path "$folder_name"

    # Convert exported BIOM to TSV
    biom convert \
        -i "$folder_name/feature-table.biom" \
        -o "$folder_name/LotusCSSP_AskovSoils_ASVtable_10_4.tsv" \
        --to-tsv
done  # <- closes the for loop

# Export filtered representative sequences to FASTA
qiime tools export \
  --input-path LotusCSSP_AskovSoils_repseqs_10_4.qza \
  --output-path LotusCSSP_AskovSoils_repseqs_10_4

mv LotusCSSP_AskovSoils_repseqs_10_4/dna-sequences.fasta \
   LotusCSSP_AskovSoils_repseqs_10_4/LotusCSSP_AskovSoils_repseqs_10_4.fasta
