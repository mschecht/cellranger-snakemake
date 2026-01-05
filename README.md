# TLDR

The tool `cellranger-snakemake` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) wrapper for [10x Cell Ranger workflows](https://www.10xgenomics.com/software) supporting [GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest), [ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), and [ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) single-cell data processing.

# Description

Reproducibility and scalability are essential components of contemporary [FAIR](https://www.nature.com/articles/sdata201618) (Findable, Accessible, Interoperable, and Reproducible) single-cell ‘omics data analysis, yet preprocessing steps lack workflow infrastructure needed to standardize large-scale and collaborative studies. 10x Genomics’ [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) is critical software for preprocessing raw single-cell ‘omics modalities, but executing it reproducibly across hundreds or thousands of samples remains cumbersome, error-prone, and computationally ineﬃcient. We present `cellranger-snakemake`, a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow wrapper that automates, scales, and standardizes Cell Ranger preprocessing for Gene Expression ([GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest)), Chromatin accessibility ([ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest)), and multiome ([ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest)) data. The workflow supports flexible input specifications, integrated logging, and portable configuration files, making it straightforward
to deploy in high-performance computing or cloud environments. By combining Snakemake’s reproducible workflow management with Cell Ranger, `cellranger-snakemake` improves reproducibility, reduces user error, and accelerates downstream single-cell ‘omics.

# Installation instructions

## 1. Set up conda

You will need [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) to install this package. 

> 💡 Conda creates many files in your home directory by default, which can exceed file quota limits on HPC systems. If you have file number restrictions on your home directory, consider installing conda in a directory without such quotas (e.g., a scratch or project directory). These commands tell conda to use specific directories (this is mainly a note to ourselves 🙂): 

```bash
$ conda config --add pkgs_dirs /path/to/miniconda3/pkgs
$ conda config --add envs_dirs /path/to/miniconda3/envs
```

To check if `conda` is installed properly, run this:
```bash
$ conda --version
conda 25.11.1
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
# remove any old cellranger-snakemake env you might have previously installed
conda env remove --name snakemake8
conda env create -f environment.yaml
conda activate snakemake8
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

If not, please install them. Once installed, you can verify the installation using the built-in version checker:

```bash
# Check all Cell Ranger tools
snakemake-run-cellranger check-versions

# Check specific workflow requirements
snakemake-run-cellranger check-versions --workflow GEX
snakemake-run-cellranger check-versions --workflow ATAC
snakemake-run-cellranger check-versions --workflow ARC
```

**Minimum supported versions:**
- Cell Ranger (GEX): 7.0.0 or higher
- Cell Ranger ATAC: 2.0.0 or higher
- Cell Ranger ARC: 2.0.0 or higher

### Linking Cell Ranger to the conda environment

Once installed, link [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), [Cell Ranger ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), [Cell Ranger ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) to the conda env:
```bash
ENV_NAME="snakemake8"

# Link Cell Ranger
SOURCE="/path/to/cellranger"
TARGET="$(conda info --base)/envs/${ENV_NAME}/bin/cellranger"

ln -s "$SOURCE" "$TARGET"

# Link Cell Ranger ATAC
SOURCE="/path/to/cellranger-atac"
TARGET="$(conda info --base)/envs/${ENV_NAME}/bin/cellranger-atac"

ln -s "$SOURCE" "$TARGET"

# Link Cell Ranger ARC
SOURCE="/path/to/cellranger-arc"
TARGET="$(conda info --base)/envs/${ENV_NAME}/bin/cellranger-arc"

ln -s "$SOURCE" "$TARGET"
```

# Quick Start

1. **Activate the environment**
```bash
conda activate snakemake8
```

2. **Verify Cell Ranger installation (optional but recommended)**
```bash
snakemake-run-cellranger check-versions
```

3. **Create a config file:**
```bash
snakemake-run-cellranger init-config --output pipeline_config.yaml
```

4. **Run a workflow:**
```bash
snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 8
```

# Tutorial

[Check out the wiki for an in-depth tutorial with publicly available practice data!
](https://github.com/mschecht/cellranger-snakemake/wiki)

# Developer notes

## Updating the Conda environment

If you add or update any Conda dependency, re-export the environment and commit the change:
```bash
conda env export > environment.yaml
```