#!/bin/bash

./reset_example.sh
cd example
mkdir -p Norfolk-1 Warrior_10
cat fastq/barcode01/*.fastq > Norfolk-1/Norfolk-1.fastq
cat fastq/barcode02/*.fastq > Warrior_10/Warrior_10.fastq
../src/reads_to_exons_concat.sh Norfolk-1/Norfolk-1.fastq Warrior_10/Warrior_10.fastq
../src/exons_concat_to_tree_imgs.sh \
    --start ../data/8_isolates_388_genes_exons.fasta.gz \
    --out_dir tree \
    --name 2_new_samples \
    */*_exons_concat.fasta
mv .metadata_264_isolates_2_new_samples.xlsx metadata_264_isolates_2_new_samples.xlsx
../src/tree_to_imgs.sh \
    --meta metadata_264_isolates_2_new_samples.xlsx \
    --out_dir tree_with_new_metadata \
    tree/RAxML_bestTree.2_new_samples.newick
