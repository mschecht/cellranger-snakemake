# config_templates.py

# Default configuration template for ARC pipeline
ARC_CONFIG = {
    "reference": "/path/to/reference-genome",
    "libraries": "/path/to/libraries_list.tsv",
    "HPC_mode": "",
    "mempercore": "",
    "normalize": "none",
    "directories_suffix": "none",
    "directories": {
        "LOGS_DIR": "00_LOGS",
        "CELLRANGERARC_COUNT_DIR": "01_CELLRANGERARC_COUNT",
        "CELLRANGERARC_AGGR_DIR": "02_CELLRANGERARC_AGGR"
    },
}

ARC_README_content = """# {workflow_type} Pipeline Configuration Guide

## Quick Start
1. Edit the `{config_filename}` file with your specific paths and settings
2. Run the pipeline: `snakemake-run-cellranger --workflow {workflow_type} --config {config_filename}`

## Required Configuration

### Reference Genome
```yaml
reference: /path/to/your/cellranger-arc-reference
```
- Must be a Cell Ranger ARC compatible reference
- Download from 10x Genomics or build your own
- Example: `/data/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0`

### libraries_list.tsv

This file contains metadata and paths for your cellranger-arc libraries
```yaml
libraries: /path/to/your/libraries_list.tsv
```
Create a TSV file with format:
```
$ cat libraries_list.tsv
batch	capture	CSV
3	3A	path/to/cellranger-arc-input-file.csv
```

- batch: Batch number for grouping captures
- capture: Capture identifier for 10x lane
- CSV: Path to the input CSV file for cellranger-arc. This CSV will be read by the cellranger-arc --libraries argument.

## HPC Configuration

### Scheduler Settings
```yaml
HPC_mode: 'slurm'  # Options: 'slurm', 'pbs', 'sge', or '' for local
mempercore: 8      # GB of RAM per CPU core
```

### Common HPC Configurations
- **SLURM cluster**: `HPC_mode: 'slurm'`
- **PBS/Torque**: `HPC_mode: 'pbs'`
- **SGE**: `HPC_mode: 'sge'`
- **Local machine**: `HPC_mode: ''`

## Output Directory Organization

### Option 1: Use default directory names

### Option 2: Use a Common Prefix
```yaml
directories_prefix: "results"
directories:
  LOGS_DIR: 00_LOGS
  CELLRANGERARC_COUNT_DIR: 01_CELLRANGERARC_COUNT
  CELLRANGERARC_AGGR_DIR: 02_CELLRANGERARC_AGGR
```
Results in: `00_LOGS_results`, `01_CELLRANGERARC_COUNT_results/`, etc.

### Option 2: Custom Directory Names
```yaml
directories_prefix: "none"
directories:
  LOGS_DIR: pipeline_logs
  CELLRANGERARC_COUNT_DIR: cellranger_counts
  CELLRANGERARC_AGGR_DIR: aggregated_results
```
Results in: `pipeline_logs/`, `cellranger_counts/`, `aggregated_results/`

## ⚠️ Important Notes

1. **Only modify the values** (the actual directory names)
2. **Use absolute paths** for reference and libraries when possible

## Example Configurations

### Minimal Setup (Local)
```yaml
reference: /data/references/arc-ref-GRCh38
libraries: ./my_samples.tsv
HPC_mode: ''
mempercore: 8
directories_prefix: results
```

### HPC Cluster Setup
```yaml
reference: /shared/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
libraries: /project/mylab/experiments/batch1/libraries.tsv
HPC_mode: 'slurm'
mempercore: 16
directories_prefix: batch1_analysis
```

Example of how to run the workflow on an HPC
```bash
 snakemake-run-cellranger --workflow ARC \
                          --config-file ARC_default_config_filled_out_2.yaml \
                          --additional-params '--cluster "sbatch -J {rule} --account={ACCOUNT} --partition={PARTITION} --ntasks=1 --cpus-per-task=12 --mem=50G --output={log}"'
```

## Troubleshooting

- **File not found errors**: Check that all paths are absolute and accessible
- **Memory errors**: Increase `mempercore` value
- **Permission errors**: Ensure write access to output directories
- **HPC submission failures**: Verify `HPC_mode` matches your cluster scheduler

## Getting Help

Run `snakemake-run-cellranger --help` for additional options and usage information.
"""

ATAC_CONFIG = {
    "samples": {
        "sample1": {
            "fastqs": "/path/to/fastqs/sample1",
            "libraries_csv": "/path/to/sample1_libraries.csv"
        }
    },
    "reference": "/path/to/cellranger-arc-reference",
    "output_dir": "results",
    "cellranger_arc": {
        "extra_args": "--force-cells=5000",
        "localcores": 8,
        "localmem": 64
    },
    "resources": {
        "default_threads": 8,
        "default_mem_gb": 64
    }
}

GEX_CONFIG = {
    "reference": "/path/to/cellranger-reference",
    "libraries_GEX": "/path/to/libraries_list_GEX.tsv",
    "HPC_mode": "",
    "chemistry": "auto",
    "mempercore": "",
    "normalize": "none",
    "directories_suffix": "none",
    "directories": {
        "LOGS_DIR": "00_LOGS",
        "CELLRANGERGEX_COUNT_DIR": "01_CELLRANGERGEX_COUNT",
        "CELLRANGERGEX_AGGR_DIR": "02_CELLRANGERGEX_AGGR"
    },
}

GEX_README_content = """# {workflow_type} Pipeline Configuration Guide

## Quick Start
1. Edit the `{config_filename}` file with your specific paths and settings
2. Run the pipeline: `snakemake-run-cellranger --workflow {workflow_type} --config {config_filename}`

## Required Configuration

### Reference Genome
```yaml
reference: /path/to/your/cellranger-reference
```
- Must be a Cell Ranger compatible reference
- Download from 10X Genomics or build your own
- Example: `/data/references/refdata-cellranger-GRCh38-2020-A-2.0.0`

### libraries_list.tsv

This file contains metadata and paths for your cellranger libraries
```yaml
libraries: /path/to/your/libraries_list.tsv
```
Create a TSV file with format:
```
$ cat libraries_list.tsv
batch	capture	sample	fastqs
B	1	ABC-B-1	/data/GEX/fastqs/
```

- batch: Batch number for grouping captures (no underscores allowed)
- capture: Capture identifier (no underscores allowed)
- sample: Prefix of the filenames of fastqs to select
- fastqs: Path(s) to input FASTQ data. If providing multiple paths, separate them with commas.

## HPC Configuration

### Scheduler Settings
```yaml
HPC_mode: 'slurm'  # Options: 'slurm', 'pbs', 'sge', or '' for local
mempercore: 8      # GB of RAM per CPU core
```

### Common HPC Configurations
- **SLURM cluster**: `HPC_mode: 'slurm'`
- **PBS/Torque**: `HPC_mode: 'pbs'`
- **SGE**: `HPC_mode: 'sge'`
- **Local machine**: `HPC_mode: ''`

## Output Directory Organization

### Option 1: Use default directory names

### Option 2: Use a Common Prefix
```yaml
directories_prefix: "results"
directories:
  LOGS_DIR: 00_LOGS
  CELLRANGERARC_COUNT_DIR: 01_CELLRANGERGEX_COUNT
  CELLRANGERARC_AGGR_DIR: 02_CELLRANGERGEX_AGGR
```
Results in: `00_LOGS_results`, `01_CELLRANGERGEX_COUNT_results/`, etc.

### Option 2: Custom Directory Names
```yaml
directories_prefix: "none"
directories:
  LOGS_DIR: pipeline_logs
  CELLRANGERARC_COUNT_DIR: cellranger_counts
  CELLRANGERARC_AGGR_DIR: aggregated_results
```
Results in: `pipeline_logs/`, `cellranger_counts/`, `aggregated_results/`

## ⚠️ Important Notes

1. **Only modify the values** (the actual directory names)
2. **Use absolute paths** for reference and libraries when possible

## Example Configurations

### Minimal Setup (Local)
```yaml
reference: /data/references/ref-GRCh38
libraries: ./my_samples.tsv
HPC_mode: ''
mempercore: 8
directories_prefix: results
```

### HPC Cluster Setup
```yaml
reference: /shared/references/refdata-cellranger-GRCh38-2020-A-2.0.0
libraries: /project/mylab/experiments/batch1/libraries.tsv
HPC_mode: 'slurm'
mempercore: 16
directories_prefix: batch1_analysis
```

## Troubleshooting

- **File not found errors**: Check that all paths are absolute and accessible
- **Memory errors**: Increase `mempercore` value
- **Permission errors**: Ensure write access to output directories
- **HPC submission failures**: Verify `HPC_mode` matches your cluster scheduler

## Getting Help

Run `snakemake-run-cellranger --help` for additional options and usage information.
"""