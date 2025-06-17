# Installation instructions

## Prerequisites

- Python 3.8+
- Conda (Miniconda or Anaconda)

## 1. Clone the repository:
```bash
git clone https://github.com/mschecht/cellranger-snakemake.git
cd cellranger-snakemake
```

## 2. Create and activate the Conda environment:
```bash
conda install mamba -n base -c conda-forge
mamba env create -f environment.yaml
conda activate cellranger-snakemake
```

## 3. Install `cellranger-snakemake` into the environment

**NOTE** 💡 If you're developing the package, use `pip install -e .` instead for an editable install.

```bash
pip install .
```

## 4. Verify installation
```bash
cellranger-snakemake --help  # or your CLI entry point
```

# Project Layout

```bash
$ tree .
.
├── cellranger_snakemake
│   ├── cli.py
│   ├── config_templates.py
│   ├── __init__.py
│   ├── utils
│   │   ├── __init__.py
│   │   ├── logger.py
│   └── workflows
│       ├── ARC
│       │   ├── cellrangerARC.smk
│       │   └── __init__.py
│       ├── ATAC
│       │   └── __init__.py
│       ├── GEX
│       │   └── __init__.py
│       ├── __init__.py
│       └── workflow.py
├── config.yaml
├── pyproject.toml
├── README.md
├── setup.py
```

# Developer notes

## Updating the Conda environment

If you add or update any Conda dependency, re-export the environment and commit the change:
```bash
conda env export --no-builds > environment.yaml
```
