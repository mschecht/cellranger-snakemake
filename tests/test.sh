#!/usr/bin/env bash

# HPC variables
PARTITION=""
NODELIST=""

#############
# Basic tests
#############
snakemake-run-cellranger --help

# Invalid config test
snakemake-run-cellranger --workflow ARC --get-default-config asdf.yaml

######################
# ARC TESTS
######################
# Test making default config
snakemake-run-cellranger --workflow ARC --get-default-config ARC.yaml

# Dry run
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml --dry-run

# DAG (DOT + PDF)
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml --dag | dot -Tpdf > dag_ARC.pdf

# Local test run
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml

# HPC runs
clusterize "snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml" \
  --output 00_LOGS/TEST_ARC.log --nodelist $NODELIST --partition $PARTITION --mem=200G

clusterize "snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml --dry-run \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=$PARTITION --ntasks=2 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_ARC_dryrun.log --nodelist $NODELIST --partition $PARTITION --mem=20G

clusterize "snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=$PARTITION --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_ARC_full.log --nodelist $NODELIST --partition $PARTITION --mem=20G

######################
# GEX TESTS
######################
# Test making default config
snakemake-run-cellranger --workflow GEX --get-default-config GEX.yaml

# Dry run
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml --dry-run

# DAG (DOT + PDF)
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml --dag | dot -Tpdf > dag_GEX.pdf

# Local test run
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml

# HPC run
clusterize "snakemake-run-cellranger --workflow GEX --config-file GEX_default_config.yaml \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=$PARTITION --ntasks=2 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_GEX.log --nodelist $NODELIST --partition $PARTITION --mem=20G


######################
# ATAC TESTS
######################
# Test making default config
snakemake-run-cellranger --workflow ATAC --get-default-config ATAC.yaml

# Dry run
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml --dry-run

# DAG (DOT + PDF)
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml --dag | dot -Tpdf > dag_ATAC.pdf

# Local test run
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml

# HPC run
clusterize "snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config.yaml \
  --additional-params \"--cluster 'sbatch -J {rule} --account=pi-lbarreiro --partition=$PARTITION --ntasks=2 --cpus-per-task=12 --mem=40G' --jobs 10\"" \
  --output 00_LOGS/TEST_ATAC.log --nodelist $NODELIST --partition $PARTITION --mem=20G

