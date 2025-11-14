# Test Data Setup Guide

This directory contains test data for the `cellranger-snakemake` pipeline for developers as well as the tutorial in our wiki. Follow the steps below to download and set up the required test datasets.

## Overview

The test data consists of [10x Genomics PBMC 3k (no cell sorting)](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0).

## Setup Instructions

### 1. Download PBMC 3k Dataset

```bash
# Create directory for downloaded files
mkdir -p 00_TEST_DATA/
cd 00_TEST_DATA/

# Download FASTQ files (tar archive)
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_fastqs.tar

# Download library CSV file
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_library.csv

# Extract the FASTQ files
tar -xf pbmc_unsorted_3k_fastqs.tar

# Clean up
rm pbmc_unsorted_3k_fastqs.tar

cd ..
```

### 2. Verify initial download dataset Structure

After extraction, your `00_TEST_DATA` directory should contain:

```bash
tree 00_TEST_DATA/pbmc_unsorted_3k/
00_TEST_DATA/pbmc_unsorted_3k/
├── atac/
│   ├── pbmc_unsorted_3k_S3_L001_I1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L002_I1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L003_I1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L004_I1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz
│   ├── pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz
│   └── pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz
└── gex/
    ├── pbmc_unsorted_3k_S01_L003_I1_001.fastq.gz
    ├── pbmc_unsorted_3k_S01_L003_I2_001.fastq.gz
    ├── pbmc_unsorted_3k_S01_L003_R1_001.fastq.gz
    ├── pbmc_unsorted_3k_S01_L003_R2_001.fastq.gz
    ├── pbmc_unsorted_3k_S01_L004_I1_001.fastq.gz
    ├── pbmc_unsorted_3k_S01_L004_I2_001.fastq.gz
    ├── pbmc_unsorted_3k_S01_L004_R1_001.fastq.gz
    └── pbmc_unsorted_3k_S01_L004_R2_001.fastq.gz
```

### 3. Re-organize the directory structure for GEX and ATAC testing

Here we will subset and rearrange the 10x multiome ARC test dataset into GEX and ATAC test datasets. 

For this, while in `path/to/cellranger-snakemake/tests/00_TEST_DATA/pbmc_unsorted_3k` (or wherever you downloaded the dataset), run the following chunck of code:

```bash
cd path/to/cellranger-snakemake/tests/00_TEST_DATA/pbmc_unsorted_3k/

# subset ATAC test data
mkdir atac_L3
mkdir atac_L4
mv atac/pbmc_unsorted_3k_S3_L003_* atac_L3/
mv atac/pbmc_unsorted_3k_S3_L004_* atac_L4/

# subset GEX data
mkdir gex_L3
mkdir gex_L4
mv gex/pbmc_unsorted_3k_S01_L003_* gex_L3/
mv gex/pbmc_unsorted_3k_S01_L004_* gex_L4/
```

## Developer Testing Guide

This section provides comprehensive testing commands for all cellranger-snakemake workflows. You can run these individually or use the automated script at `test.sh`.

### Prerequisites

Before running tests, ensure:
1. Test data is downloaded and organized (see steps 1-3 above)
2. `cellranger-snakemake` is installed and `snakemake-run-cellranger` is in your PATH
3. For HPC testing: Update `PARTITION` and `NODELIST` variables in `test.sh`

### Basic Validation Tests

```bash
cd path/to/cellranger-snakemake/tests/

# Test CLI accessibility
snakemake-run-cellranger --help

# Test invalid config handling
snakemake-run-cellranger --workflow ARC --get-default-config 00_TEST_DATA/asdf.yaml
```

### ARC Workflow Tests

```bash
# Generate default configuration
snakemake-run-cellranger --workflow ARC --get-default-config ARC.yaml

# Dry run (validate workflow without execution)
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml --dry-run

# Generate workflow DAG visualization
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml --dag | dot -Tpdf > dag_ARC.pdf

# Local execution (single machine)
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml

# HPC execution examples
# Single-node run
clusterize "snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml" \
  --output 00_LOGS/TEST_ARC.log --nodelist $NODELIST --partition $PARTITION --mem=200G

# Multi-job cluster run (dry-run)
clusterize "snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml --dry-run \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=\$PARTITION --ntasks=2 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_ARC_dryrun.log --nodelist $NODELIST --partition $PARTITION --mem=20G

# Multi-job cluster run (full execution)
clusterize "snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=\$PARTITION --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_ARC_full.log --nodelist $NODELIST --partition $PARTITION --mem=20G
```

### GEX Workflow Tests

```bash
# Generate default configuration
snakemake-run-cellranger --workflow GEX --get-default-config GEX.yaml

# Dry run
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml --dry-run

# Generate workflow DAG
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml --dag | dot -Tpdf > dag_GEX.pdf

# Local execution
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml

# HPC execution
clusterize "snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=\$PARTITION --ntasks=2 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_GEX.log --nodelist $NODELIST --partition $PARTITION --mem=20G
```

### ATAC Workflow Tests

```bash
# Generate default configuration
snakemake-run-cellranger --workflow ATAC --get-default-config ATAC.yaml

# Dry run
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml --dry-run

# Generate workflow DAG
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml --dag | dot -Tpdf > dag_ATAC.pdf

# Local execution
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml

# HPC execution
clusterize "snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=\$PARTITION --ntasks=2 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_ATAC.log --nodelist $NODELIST --partition $PARTITION --mem=20G
```

### Automated Testing

For convenience, all tests can be run using the automated script:

```bash
# Edit HPC variables if needed
vim test.sh  # Update PARTITION and NODELIST variables

# Run all tests
bash test.sh
```
