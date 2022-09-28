#!/bin/bash

./reset_example.sh
cd example
mkdir -p Norfolk-1 Warrior_10
cat fastq/barcode01/*.fastq > Norfolk-1/Norfolk-1.fastq
cat fastq/barcode02/*.fastq > Warrior_10/Warrior_10.fastq
../src/pipeline.sh Norfolk-1/Norfolk-1.fastq Warrior_10/Warrior_10.fastq 
