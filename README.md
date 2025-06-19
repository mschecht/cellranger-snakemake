# Description

The tool `snakemake-run-cellranger` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) wrapper for [10x Cell Ranger workflows](https://www.10xgenomics.com/software) supporting [GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest), [ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), and [ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) single-cell data processing.

# Installation instructions

## 1. Set up conda

You will need [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) to install the package. 

To check if `conda` is installed properly, run this:
```bash
$ conda --version
conda 23.7.4
```

Once you have confirmed you have conda installed, run this command to make sure you are up-to-date:
```bash
conda update conda
```

## 2. Clone the repository:
```bash
git clone https://github.com/mschecht/cellranger-snakemake.git
cd cellranger-snakemake
```

## 3. Create and activate the Conda environment:
```bash
conda env create -f environment_mini.yaml
conda activate cellranger-snakemake
```

## 4. Install `cellranger-snakemake` into the environment

> 💡 If you're developing the package, use `pip install -e .` instead for an editable install.

```bash
pip install .
```

Verify `snakemake-run-cellranger` installation:
```bash
snakemake-run-cellranger --help  # or your CLI entry point
```

...and a few other tools while we are at it:
```bash
# Check Snakemake
snakemake --version
bcftools --version
samtools --version
```

## 5. Cell Ranger Installation

This package provides a Snakemake wrapper for [10x Genomics Cell Ranger](https://www.10xgenomics.com/software), but it does not include the Cell Ranger software itself.

Check if you have [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), [Cell Ranger ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), [Cell Ranger ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) installed: 

```bash
$ cellranger --version
cellranger cellranger-7.0.1

$ cellranger-atac --version
cellranger-atac cellranger-atac-2.1.0

$ cellranger-arc --version
cellranger-arc cellranger-arc-2.0.2
```

If not, please install it :)

Once installed, link [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), [Cell Ranger ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), [Cell Ranger ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) to the conda env:
```bash
# Link Cell Ranger
SOURCE="/path/to/cellranger"
ENV_NAME="cellranger-snakemake"
TARGET="$(conda info --base)/envs/${ENV_NAME}/bin/cellranger"

ln -s "$SOURCE" "$TARGET"

# Link Cell Ranger ATAC
SOURCE="/path/to/cellranger-atac"
ENV_NAME="cellranger-snakemake"
TARGET="$(conda info --base)/envs/${ENV_NAME}/bin/cellranger-atac"

ln -s "$SOURCE" "$TARGET"

# Link Cell Ranger ARC
SOURCE="/path/to/cellranger-arc"
ENV_NAME="cellranger-snakemake"
TARGET="$(conda info --base)/envs/${ENV_NAME}/bin/cellranger-arc"

ln -s "$SOURCE" "$TARGET"
```

## 6. Verification of installation
- FIXME: add a verification step to ensure everything is set up correctly

# Quick Start

1. **Activate the environment**
```bash
conda activate cellranger-snakemake
```

2. **Create a config file:**
```bash
snakemake-run-cellranger --workflow ARC --get-default-config
```

3. **Run a workflow:**
```bash
snakemake-run-cellranger --workflow ARC --config-file ARC_default_config.yaml
```

# Developer notes

## Updating the Conda environment

If you add or update any Conda dependency, re-export the environment and commit the change:
```bash
conda env export > environment.yaml
```
