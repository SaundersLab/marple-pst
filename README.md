<img alt="MARPLE DIAGNOSTICS" src="docs/marple-green.png" width="100%"/>

Pipeline for analysis of Puccinia striiformis f.sp. tritici genomic reads

## Install

```bash
./install/install.sh
```

Before installing, please ensure the licences of miniconda and each software dependency in `env.yml` are compatible with your usage.

- Installs conda (miniconda)
- Installs dependencies from `eny.yml` into the environment `marple-pst`
- Runs tests with `test/run_tests.sh`

## Tutorial

Let's go through the pipeline using a fictional example.

In this example, 2 samples have already been sequenced (using nanpore sequencing) and basecalled using Guppy.

Our starting point is a directory for each sample containing the fastq files created by Guppy.

|sample name|reads directory        |
|-----------|-----------------------|
|Norfolk-1  |example/fastq/barcode01|
|Warrior_10 |example/fastq/barcode02|

1. Reset the example

    ```bash
    ./reset_example.sh
    ```

2. Navigate to the example directory and make a new directory for each sample.

    ```bash
    cd example
    mkdir -p Norfolk-1 Warrior_10
    ```

3. We need a single fastq for each sample, so combine all the fastq files for barcode01 into one. Do the same for barcode02.

    ```bash
    cat fastq/barcode01/*.fastq > Norfolk-1/Norfolk-1.fastq
    cat fastq/barcode02/*.fastq > Warrior_10/Warrior_10.fastq
    ```

4. Align the reads for each sample to the reference genes, extract the exons, and create a report

    ```bash
    ../src/reads_to_exons_concat.sh Norfolk-1/Norfolk-1.fastq Warrior_10/Warrior_10.fastq
    ```

5. Inspect the report by opening `example/multiqc_report.html` in your browser

6. Use the extracted and concatenated exons of the new samples to create a tree

    ```bash
    ../src/exons_concat_to_tree_imgs.sh \
        --start ../data/8_isolates_388_genes_exons.fasta.gz \
        --out_dir tree \
        --name 2_new_samples \
        */*_exons_concat.fasta
    ```

    In this example we use a starting tree input containing only 8 isolates to save time.

    When you run this with real samples, use `../data/56_isolates_388_genes_exons.fasta.gz` which
    contains 56 isolates.

7. Inspect the tree by opening `example/tree/2_new_samples_country.pdf` in your browser.

    Notice how country is shown as '?' for the new samples. We'll fix that in the next step.

8. Copy the metadata spreadsheet `data/metadata_264_isolates.xlsx` and rename it so that you have `example\metadata_264_isolates_2_new_samples.xlsx`.

    Add the new sample metadata to the spreadsheet so that the top 3 rows look like this:

    |tree_name |tree_new_name|country|nextstrain_group|region|author|other_names|genetic_group|
    |----------|-------------|-------|----------------|------|------|-----------|-------------|
    |Norfolk-1 |Norfolk-1    |UK     |                |Europe|      |           |             |
    |Warrior_10|Warrior_10   |UK     |                |Europe|      |           |             |

9. Use this metadata to make a new visualisation of the tree

    ```bash
    ../src/tree_to_imgs.sh \
        --meta metadata_264_isolates_2_new_samples.xlsx \
        --out_dir tree_with_new_metadata \
        tree/RAxML_bestTree.2_new_samples.newick
    ```

10. Inspect the new visualisation of the tree by opening `example/tree_with_new_metadata/2_new_samples_country.pdf` in your browser.

    The country of our 2 new samples is now showing correctly.

If you get stuck and want to see what the output should look like, then run this from the marple-pst directory:

```bash
./run_example.sh
```

## Pipeline

See [pipeline.md](docs/pipeline.md) for a flowchart of the pipeline.

## Test

Run all the tests (roughly 2 minutes):

```bash
./test/run_tests.sh
```

Run the tests and report code coverage:

```bash
./test/run_tests.sh coverage
```

Skip the end to end test and integration tests

```bash
./test/run_tests.sh skip_integration skip_end_to_end
```

## Uninstall

```bash
./install/uninstall.sh
```

- Deletes the `marple-pst` and its packacges
- Deletes `marple_pst_miniconda.sh` and `marple_pst_miniconda`
