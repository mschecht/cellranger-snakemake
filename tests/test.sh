#!/usr/bin/env bash

# HPC variables
ACCOUNT=""
PARTITION=""
NODELIST=""

#############
# Basic tests
#############
snakemake-run-cellranger --help

# Invalid config test
snakemake-run-cellranger --workflow ARC --get-default-config asdf

######################
# ARC TESTS
######################
# Test making default config
snakemake-run-cellranger --workflow ARC --get-default-config ARC.yaml

# Dry run
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config_filled_out.yaml --dry-run

# DAG (DOT + PDF)
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config_filled_out.yaml --dag | dot -Tpdf > dag_ARC.pdf

# Local test run
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config_filled_out.yaml

# HPC runs
snakemake-run-cellranger --workflow ARC \
                         --config-file ARC_default_config_filled_out.yaml \
                         --snakemake-args \"--cluster 'sbatch -J {rule} --account=$ACCOUNT --partition=$PARTITION --nodelist $NODELIST --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"$

######################
# GEX TESTS
######################
# Test making default config
snakemake-run-cellranger --workflow GEX --get-default-config GEX.yaml

# Dry run
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config_filled_out.yaml --dry-run

# DAG (DOT + PDF)
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config_filled_out.yaml --dag | dot -Tpdf > dag_GEX.pdf

# Local test run
snakemake-run-cellranger --workflow GEX --config-file GEX_default_config_filled_out.yaml

# HPC run
snakemake-run-cellranger --workflow GEX \
                         --config-file GEX_default_config_filled_out.yaml \
                         --snakemake-args \"--cluster 'sbatch -J {rule} --account=$ACCOUNT --partition=$PARTITION --nodelist $NODELIST --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"$


######################
# ATAC TESTS
######################
# Test making default config
snakemake-run-cellranger --workflow ATAC --get-default-config ATAC.yaml

# Dry run
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config_filled_out.yaml --dry-run

# DAG (DOT + PDF)
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config_filled_out.yaml --dag | dot -Tpdf > dag_ATAC.pdf

# Local test run
snakemake-run-cellranger --workflow ATAC --config-file ATAC_default_config_filled_out.yaml

# HPC run
snakemake-run-cellranger --workflow ATAC \
                         --config-file ATAC_default_config_filled_out.yaml \
                         --snakemake-args \"--cluster 'sbatch -J {rule} --account=$ACCOUNT --partition=$PARTITION --nodelist $NODELIST --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10\"$