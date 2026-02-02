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

