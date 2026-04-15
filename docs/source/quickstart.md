# Quick Start

## 1. Activate the environment

```bash
conda activate snakemake8
```

## 2. Verify Cell Ranger installation (optional)

```bash
sc-preprocess check-versions
```

## 3. Create a config file

```bash
sc-preprocess init-config --output pipeline_config.yaml
```

## 4. Dry run (recommended)

Preview what the pipeline will do without executing anything:

```bash
sc-preprocess run --config-file pipeline_config.yaml --cores 1 --dry-run
```

Optionally, visualize the workflow as a DAG:

```bash
sc-preprocess run --config-file pipeline_config.yaml --cores 1 --dag | dot -Tpng > dag.png
```

## 5. Run the pipeline

```bash
# Local execution
sc-preprocess run --config-file pipeline_config.yaml --cores 8

# HPC execution (with SLURM profile)
sc-preprocess run --config-file pipeline_config.yaml \
                             --cores all \
                             --snakemake-args --profile HPC_profiles
```
