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

> 💡 If you need to start fresh, first remove the old environment: `conda env remove --name snakemake8`

```bash
conda env create -f environment.yaml
conda activate snakemake8
```

Verify the environment tools are working:
```bash
snakemake --version
bcftools --version
samtools --version
```

## 4. Install `cellranger-snakemake` into the environment

```bash
pip install -e .
```

> 💡 The `-e` (editable) flag means code changes take effect immediately without reinstalling. You still need to re-run `pip install -e .` after pulling changes that modify `pyproject.toml` (new dependencies or CLI commands).

Verify `snakemake-run-cellranger` installation:
```bash
snakemake-run-cellranger --help
```

## 5. Cell Ranger Installation

This package provides a Snakemake wrapper for [10x Genomics Cell Ranger](https://www.10xgenomics.com/software), but it does not include the Cell Ranger software itself.

Check if you have [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), [Cell Ranger ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), [Cell Ranger ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) installed: 

```bash
$ cellranger --version
cellranger cellranger-9.0.1

$ cellranger-atac --version
cellranger-atac cellranger-atac-2.1.0

$ cellranger-arc --version
cellranger-arc cellranger-arc-2.0.2
```
If not, please install them from the [10x Genomics download page](https://www.10xgenomics.com/support/software/cell-ranger/latest).

**Minimum supported versions:**
- Cell Ranger (GEX): 9.0.0 or higher
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

Once installed, you can verify the installation using the built-in version checker:

```bash
# Check all Cell Ranger tools
snakemake-run-cellranger check-versions

# Check specific workflow requirements
snakemake-run-cellranger check-versions --workflow GEX
snakemake-run-cellranger check-versions --workflow ATAC
snakemake-run-cellranger check-versions --workflow ARC
```

You're all set! Head to the [Quick Start](quickstart.md) to run your first workflow.