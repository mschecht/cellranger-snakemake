# Installation instructions

## 

1. Clone the repository:
```bash
git clone https://github.com/yourusername/cellranger-snakemake.git
cd cellranger-snakemake
```
2. Create and activate the Conda environment:
```bash
conda env create -f environment.yaml
conda activate cellranger-snakemake
```

3. Verify installation
```bash
cellranger-snakemake --help  # or your CLI entry point

```

# cellranger-snakemake

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
