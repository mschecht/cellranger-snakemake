# Test Data Setup Guide

This directory contains test data for the `cellranger-snakemake` pipeline for developers as well as the tutorial in our wiki. Follow the steps below to download and set up the required test datasets.

## Overview

The test data consists of [10x Genomics PBMC 3k (no cell sorting)](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0).

## Setup Instructions

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
snakemake-run-cellranger --workflow ARC --get-default-config asdf
```

### ARC Workflow Tests

```bash
# Generate default configuration
snakemake-run-cellranger --workflow ARC --get-default-config ARC.yaml

# Dry run (validate workflow without execution)
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config_filled_out.yaml --dry-run

# Generate workflow DAG visualization
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config_filled_out.yaml --dag | dot -Tpdf > dag_ARC.pdf

# Local execution (single machine)
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config_filled_out.yaml

# HPC execution
snakemake-run-cellranger --workflow ARC \
                         --config-file ARC_default_config_filled_out.yaml \
                         --snakemake-args \"--cluster 'sbatch -J {rule} --account=$ACCOUNT --partition=$PARTITION --nodelist $NODELIST --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"$
```

### GEX Workflow Tests

```bash
# Generate default configuration
snakemake-run-cellranger --workflow GEX --get-default-config GEX.yaml

# Dry run
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config_filled_out.yaml --dry-run

# Generate workflow DAG
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config_filled_out.yaml --dag | dot -Tpdf > dag_GEX.pdf

# Local execution
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config_filled_out.yaml

# HPC execution
snakemake-run-cellranger --workflow GEX \
                         --config-file GEX_default_config_filled_out.yaml \
                         --snakemake-args \"--cluster 'sbatch -J {rule} --account=$ACCOUNT --partition=$PARTITION --nodelist $NODELIST --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"$
```

### ATAC Workflow Tests

```bash
# Generate default configuration
snakemake-run-cellranger --workflow ATAC --get-default-config ATAC.yaml

# Dry run
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config_filled_out.yaml --dry-run

# Generate workflow DAG
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config_filled_out.yaml --dag | dot -Tpdf > dag_ATAC.pdf

# Local execution
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config_filled_out.yaml

# HPC execution
snakemake-run-cellranger --workflow ATAC \
                         --config-file ATAC_default_config_filled_out.yaml \
                         --snakemake-args \"--cluster 'sbatch -J {rule} --account=$ACCOUNT --partition=$PARTITION --nodelist $NODELIST --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"$
```

### Automated Testing

For convenience, all tests can be run using the automated script:

```bash
# Edit HPC variables if needed
vim test.sh  # Update PARTITION and NODELIST variables

# Run all tests
bash test.sh
```
