#!/bin/bash

# Bash strict mode
set -euo pipefail

threads=1
trim="yes"
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        --threads) # Number of threads to use
            threads="$2" ; shift ;;
        --trim) # Skip the trimming step
            trim="$2" ; shift ;;
        *) echo "ERROR: Unkown option: $1 " >&2
        exit 1
        ;;
    esac
    shift
done

./reset_example.sh
cd example
mkdir -p Norfolk-1 Warrior_10
cat fastq/barcode01/*.fastq > Norfolk-1/Norfolk-1.fastq
cat fastq/barcode02/*.fastq > Warrior_10/Warrior_10.fastq
../src/reads_to_exons_concat.sh --threads "$threads" --trim "$trim" Norfolk-1/Norfolk-1.fastq Warrior_10/Warrior_10.fastq
../src/exons_concat_to_tree_imgs.sh \
    --start ../data/8_isolates_388_genes_exons.fasta.gz \
    --out_dir tree \
    --name 2_new_samples \
    --threads "$threads" \
    */*_exons_concat.fasta
cp .metadata_264_isolates_2_new_samples.xlsx metadata_264_isolates_2_new_samples.xlsx
../src/tree_to_imgs.sh \
    --meta metadata_264_isolates_2_new_samples.xlsx \
    --out_dir tree_with_new_metadata \
    tree/RAxML_bestTree.2_new_samples.newick
