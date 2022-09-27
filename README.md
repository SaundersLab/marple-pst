<img alt="MARPLE DIAGNOSTICS" src="docs/marple-green.png" width="100%"/>

<!-- Pipeline for analysis of Puccinia striiformis f.sp. tritici genomic reads -->

## Install

```bash
./install.sh
```

- Installs conda (miniconda)
- Installs dependencies from `eny.yml` into the environment `marple-pst`
- Runs tests with `test/run_tests.sh`

## Tutorial

The pipeline steps are worked through in this fictional example.

2 samples were sequenced:

|sample name|reads directory  |
|-----------|-----------------|
|Norfolk-1  |example/barcode01|
|Warrior_10 |example/barcode02|

1. Reset the example if neccessary

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

4. Run the pipeline

    ```bash
    python3 ../src/pipeline.py Norfolk-1/Norfolk-1.fastq Warrior_10/Warrior_10.fastq 
    ```

## Test

Run the tests:

```bash
./test/run_tests.sh
```

Run the tests and report code coverage:

```bash
./test/run_tests.sh COVERAGE
```

Only run the unit tests:

```bash
./test/run_tests.sh SKIP_INTEGRATION
```

## Uninstall

```bash
./uninstall.sh
```

- Deletes the `marple-pst` and its packacges
- Deletes `marple_pst_miniconda.sh` and `marple_pst_miniconda`
- Does not uninstall conda
