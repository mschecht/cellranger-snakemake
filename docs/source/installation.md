# Installation

## Prerequisites

- [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) or [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (recommended for faster solves)
- Git

## 1. Clone the repository

```bash
git clone https://github.com/mschecht/cellranger-snakemake.git
cd cellranger-snakemake
```

## 2. Create the conda environment

A single command installs all dependencies and the `cellranger-snakemake` CLI:

```bash
conda env create -f environment.yaml
```

> 📌 **Tip:** If you have [Mamba](https://mamba.readthedocs.io) installed, use `mamba env create -f environment.yaml` for faster dependency resolution.

> 📌 **Starting fresh?** Remove an existing environment first: `conda env remove --name snakemake8`

## 3. Verify the installation

```bash
conda activate snakemake8
```

Check that the CLI and key tools are working:

```bash
# Pipeline CLI
snakemake-run-cellranger --help

# Core tools
snakemake --version     # should print 9.x
bcftools --version      # should print 1.22+
samtools --version      # should print 1.22+
```

Check that Python packages import correctly:

```bash
python -c "import scanpy; import anndata; import muon; import celltypist; import decoupler; import scrublet; print('All imports OK')"
```

### Per-rule conda environments

Notes for developers: Rules that use **snapatac2** (ATAC analysis) and **scvi-tools** (SOLO, scANVI) run in their own isolated conda environments, created automatically by Snakemake. These environments are built once on first use and cached for subsequent runs.

No extra installation steps are needed — the CLI passes `--use-conda` to Snakemake automatically.

## 4. Cell Ranger installation

This package wraps [10x Genomics Cell Ranger](https://www.10xgenomics.com/software) but does **not** include the Cell Ranger software itself.

Check if Cell Ranger is already available:

```bash
cellranger --version         # GEX: 9.0.0+
cellranger-atac --version    # ATAC: 2.0.0+
cellranger-arc --version     # ARC: 2.0.0+
```

If not installed, download from the [10x Genomics download page](https://www.10xgenomics.com/support/software/cell-ranger/latest).

### Linking Cell Ranger to the conda environment

After installing Cell Ranger, symlink the executables into the conda environment so Snakemake can find them:

```bash
ENV_NAME="snakemake8"
CONDA_BIN="$(conda info --base)/envs/${ENV_NAME}/bin"

# Replace /path/to/ with your actual Cell Ranger installation paths
ln -s /path/to/cellranger      "$CONDA_BIN/cellranger"
ln -s /path/to/cellranger-atac "$CONDA_BIN/cellranger-atac"
ln -s /path/to/cellranger-arc  "$CONDA_BIN/cellranger-arc"
```

> **Note:** You only need to link the Cell Ranger tools for the modalities you plan to use (e.g., GEX only needs `cellranger`).

Verify with the built-in version checker:

```bash
snakemake-run-cellranger check-versions
snakemake-run-cellranger check-versions --workflow GEX   # check specific workflow
```

## Troubleshooting

### Conda is slow

Use [Mamba](https://mamba.readthedocs.io) as a drop-in replacement for faster solves:

```bash
conda install -n base -c conda-forge mamba
mamba env create -f environment.yaml
```

### HPC file quota exceeded

Conda creates many files in your home directory by default, which can exceed file quota limits. This is a reoccurrent problem for HPC users. If your home directory has file number restrictions, redirect conda's cache to a directory without quotas:

```bash
conda config --add pkgs_dirs /path/to/project_dir/conda_pkgs
conda config --add envs_dirs /path/to/project_dir/conda_envs
```

### Git pulling the development version (advanced)

After pulling new changes from the development version of the workflow, update your environment like this:

```bash
conda activate snakemake8
pip install -e .
```

If `environment.yaml` has changed, recreate the environment:

```bash
conda env remove --name snakemake8
conda env create -f environment.yaml
```

If per-rule conda environment files (`workflows/envs/*.yaml`) have changed, delete the Snakemake conda cache so they get rebuilt on next run:

```bash
rm -rf .snakemake/conda
```

---

Next: [Quick Start](quickstart.md)
