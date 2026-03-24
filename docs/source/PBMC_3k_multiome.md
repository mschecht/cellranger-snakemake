# Case study: 3K PBMC Epi Multiome ATAC + Gene Expression

The [PBMC from a Healthy Donor - No Cell Sorting (3k)](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0) dataset is a great way to explore cell cellranger snakemake 

## Learning objectives

The goal of this tutorial is to show you how to get started preprocessing your own [10X Epi Multiome ATAC + Gene Expression data](https://www.10xgenomics.com/support/epi-multiome). 

1. Learn how to set up the configuration file

2. Compose the input files

3. Run the workflow

4. Explore the resulting directory structure and where you can find processed data and metadata.

5. Import the pre-processed data into Scanpy, Suerat, or ArchR to get started with data analysis. 

## Download input data

Here we will download the [3K PBMC Epi Multiome ATAC + Gene Expression input data](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0) and the [human reference genome for Cell Ranger ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads#reference-downloads):

Download 10X Epi Multiome ATAC + Gene Expression data
```bash
mkdir tests/PBMC_3K_MULTIOME && cd tests/PBMC_3K_MULTIOME

# Input Files: 3K PBMC Epi Multiome ATAC + Gene Expression
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_fastqs.tar
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_library.csv

# Human reference genome
wget "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2024-A.tar.gz"
```

Extract files
```bash
tar -xvf pbmc_unsorted_3k_fastqs.tar
tar -xvf refdata-cellranger-arc-GRCh38-2024-A.tar.gz
```

## Set up the input files

Initialize a config file: `pipeline_config.yaml`
```bash
$ snakemake-run-cellranger init-config

Single-Cell Preprocessing Pipeline Configuration Generator

Project name: 3K_PBMC_MULTIOME_PROCESSED
Output directory (output): 3K_PBMC_MULTIOME_PROCESSED

Resource Configuration
Total memory (GB) (32):
Temporary directory ():
Directory suffix (or 'none') (none):

Cell Ranger Configuration
Enable Cell Ranger GEX? [y/n] (n): n
Enable Cell Ranger ATAC? [y/n] (n): n
Enable Cell Ranger ARC? [y/n] (n): y
  Reference genome path:
  Libraries TSV path:
  Normalization method [none/depth] (none):

Enable demultiplexing? [y/n] (n): n

Enable doublet detection? [y/n] (n): y
  Doublet detection method [scrublet/solo] (scrublet):
    Expected doublet rate (0.06):
    Min genes per cell (filter_cells_min_genes) (100):
    Min cells per gene (filter_genes_min_cells) (3):
    Min gene variability percentile (85.0):
    Number of principal components (30):

Enable cell type annotation? [y/n] (n): n

Saving configuration to 'pipeline_config.yaml'...

✓ Configuration saved to: pipeline_config.yaml
[INFO] Enabled steps: cellranger_arc, doublet_detection
```

Make a `libraries.tsv`
```bash
echo -e "batch\tcapture\tCSV" > libraries_list.tsv
echo -e "1\tA\tpbmc_unsorted_3k_library.csv" >> libraries_list.tsv
```

Update the paths in `pbmc_unsorted_3k_library.csv` to absolute paths otherwise Cell Ranger ARC will be upset:
```bash
FASTQ_DIR=$(realpath pbmc_unsorted_3k)
sed -i "s|pbmc_unsorted_3k/gex|${FASTQ_DIR}/gex|; s|pbmc_unsorted_3k/atac|${FASTQ_DIR}/atac|" pbmc_unsorted_3k_library.csv
```

Finish filling out the `pipeline_config.yaml` with paths to necessary files e.g. `libraries_list.tsv` and reference genome path:

```yaml
project_name: 3K_PBMC_MULTIOME_PROCESSED
output_dir: 3K_PBMC_MULTIOME_PROCESSED
samples: {}
resources:
  mem_gb: 32
  tmpdir: ''
directories_suffix: none
cellranger_arc:
  enabled: true
  reference: /path/to/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
  libraries: libraries_list.tsv
  normalize: none
  directories:
    LOGS_DIR: 00_LOGS
  threads: 10
  mem_gb: 64
doublet_detection:
  enabled: true
  method: scrublet
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    random_state: 0
```

## Run the tool 

### Dry-run

Before running the workflow it's best practice to run a [dry-run](https://snakemake.readthedocs.io/en/stable/executing/cli.html#useful-command-line-arguments). This is a Snakemake command that will test the workflow before running it to see which jobs will be run. This will print information for every job Snakemake plans on running. The most informative part for us is the `Job stats` section which we highlight below. This counts how many time individual Rules will be ran and is a great sanity check. For example, if you have 3 captures, then you should see Rule `cellranger_arc_count` being run three times:

```bash
# Read about this command
snakemake-run-cellranger run -h

# Dry run
$ snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 1 --dry-run
[INFO] Config validated. Enabled steps: cellranger_arc, doublet_detection
[INFO] Running Snakemake with command: snakemake --snakefile /project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/cellranger_snakemake/workflows/main.smk --configfile pipeline_config.yaml --cores 1 --use-conda --dry-run
[INFO] ============================================================
[INFO] Single-Cell Preprocessing Pipeline
[INFO] ============================================================
[INFO] Project: 3K_PBMC_MULTIOME_PROCESSED
[INFO] Output directory: 3K_PBMC_MULTIOME_PROCESSED
[INFO] Enabled steps: cellranger_arc, doublet_detection
[INFO] ============================================================
[WARNING] Row 2 in 'libraries_list.tsv': Converting relative path 'pbmc_unsorted_3k_library.csv' to absolute path '/project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/tests/PBMC_3K_MULTIOME/pbmc_unsorted_3k_library.csv'
[WARNING] Row 2 in 'libraries_list.tsv': Converting relative path 'pbmc_unsorted_3k_library.csv' to absolute path '/project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/tests/PBMC_3K_MULTIOME/pbmc_unsorted_3k_library.csv'
[INFO] libraries_list.tsv file format is valid.
[INFO] libraries_list.tsv file format is valid.
[INFO] Cell Ranger ARC: Found 1 capture(s) across 1 batch(es)
[INFO] Cell Ranger ARC: Found 1 capture(s) across 1 batch(es)
[INFO] Batch aggregation: Found 1 ARC batch(es)
[INFO] Batch aggregation: Found 1 ARC batch(es)
[INFO] Doublet Detection: Using scrublet method
[INFO] Doublet Detection: Using scrublet method
host: midway3-login4.rcc.local
Building DAG of jobs...
Job stats:
job                     count
--------------------  -------
aggregate_arc_batch         1
all                         1
cellranger_arc_aggr         1
cellranger_arc_count        1
create_arc_mudata           1
enrich_arc_metadata         1
run_scrublet                1
total                       7

...
```

### Visualize the workflow with a DAG file

Our favorite way to check the workflow before starting or for debugging is to examine the DAG file. This image show a network of jobs and dependencies for the workflow. Each node is a job and each arrow represents a dependent rule.

> 📌 **Note**: If the rules are circles then the rule has not been run yet, however, if the rules are bordered with dotted lines then it's been completed. This distinction is valuable when examining an incomplete workflow. 

```bash
snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 1 --dag | dot -Tpng > dag_arc_3k.png
```

:::{figure} _images/dag_arc_3k_incomplete.png
:alt: DAG for ARC 3k PBMC pipeline
:width: 80%

**Multiome Pipeline DAG** — DAG file showing all rules and their dependencies.
:::

Here we will break down the meaning of each rule so you can keep track of what's going on. 

**cellranger_arc_count**: Runs `cellranger-arc count` per capture, aligning GEX and ATAC reads to the reference genome and producing a joint feature-barcode matrix.

**create_arc_mudata**: Converts data from the Cell Ranger ARC output to per-capture [MuData object](https://mudata.readthedocs.io/stable/) (`.h5mu`), adding traceability metadata (`batch_id`, `capture_id`, `cell_id`).

**cellranger_arc_aggr**: Runs `cellranger-arc aggr` which aggregates all per-capture Cell Ranger ARC outputs within a batch into a single normalized count matrix.

**aggregate_arc_batch**: Merges all per-capture MuData objects into a single batch-level `.h5mu` file, verifying `cell_id` uniqueness across captures.

**run_scrublet**: Runs Scrublet doublet detection on the GEX modality of each per-capture MuData object, adding doublet scores and predictions to cell metadata.

**enrich_arc_metadata**: Joins all downstream preprocessing metadata from demultiplexing and doublet detection into the batch-level MuData object.

**all**: Final Snakemake rule that collects all expected outputs to ensure the full workflow is completed.

### Local Execution

```bash
# Local execution
snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 1
```


### Snakemake arguments

We added the parameter `--snakemake-args` to send arguments straight to `Snakemake`!

For example, a popular `Snakemake` argument is `--keep-going`, where `Snakemake` will continue running jobs even if one fails. Please note that is MUST be the last argument in the command. Here is what is looks like in practice:

```bash
snakemake-run-cellranger run --config-file pipeline_config.yaml \
                             --cores 1 \
                             --dry-run \
                             --snakemake-args --keep-going
```

Another useful argument is `--forcerun`, which forces Snakemake to re-execute a specific rule and all rules that depend on it — without re-running expensive upstream steps like Cell Ranger. This is handy when you update a script and only want to reprocess from that point forward:

```bash
snakemake-run-cellranger run --config-file pipeline_config.yaml \
                             --cores 1 \
                             --snakemake-args "--forcerun create_arc_mudata aggregate_arc_batch"
```

> 📌 **Note**: You can pass multiple rule names to `--forcerun`. Snakemake will automatically re-run all downstream rules that depend on the forced rules.

### Run jobs in parallel!

Make a Snakemake SLURM [configuration file](https://snakemake.readthedocs.io/en/v7.19.1/executing/cli.html#profiles)

> 📌 **Note**: Replace the BASH variables `SLURM_ACCOUNT` and `SLURM_PARTITION` with your SLURM appropriate setting before running the script below.

```bash
# Set your SLURM account and partition
SLURM_ACCOUNT=""   # <- replace with your account
SLURM_PARTITION=""    # <- replace with your partition

mkdir -p HPC_profiles
cat > HPC_profiles/config.yaml << EOF
executor: slurm
jobs: 10
default-resources:
- slurm_account=${SLURM_ACCOUNT}
- slurm_partition=${SLURM_PARTITION}
- runtime=720
retries: 2
latency-wait: 60
printshellcmds: true
keep-going: true
rerun-incomplete: true
EOF
```


```bash
# HPC execution - `--cores all` tell snakemake to use the `threads` assigned to each rule.
snakemake-run-cellranger run --config-file pipeline_config.yaml \
                             --cores all \
                             --snakemake-args --profile HPC_profiles --keep-going
```


## Interpreting STDOUT

After starting the program you should see an output that looks like this, let's break it down:

```console
$ snakemake-run-cellranger run --config-file pipeline_config.yaml \ 
                               --cores all \
                               --snakemake-args --profile HPC_profiles
[INFO] Config validated. Enabled steps: cellranger_arc, doublet_detection
[INFO] Running Snakemake with command: snakemake --snakefile /project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/cellranger_snakemake/workflows/main.smk --configfile pipeline_config.yaml --cores all --use-conda --profile HPC_profiles
Using profile HPC_profiles for setting default command line arguments.
[INFO] ============================================================
[INFO] Single-Cell Preprocessing Pipeline
[INFO] ============================================================
[INFO] Project: 3K_PBMC_MULTIOME_PROCESSED
[INFO] Output directory: 3K_PBMC_MULTIOME_PROCESSED
[INFO] Enabled steps: cellranger_arc, doublet_detection
[INFO] ============================================================
[WARNING] Row 2 in 'libraries_list.tsv': Converting relative path 'pbmc_unsorted_3k_library.csv' to absolute path '/project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/tests/PBMC_3K_MULTIOME/pbmc_unsorted_3k_library.csv'
[WARNING] Row 2 in 'libraries_list.tsv': Converting relative path 'pbmc_unsorted_3k_library.csv' to absolute path '/project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/tests/PBMC_3K_MULTIOME/pbmc_unsorted_3k_library.csv'
[INFO] libraries_list.tsv file format is valid.
[INFO] libraries_list.tsv file format is valid.
[INFO] Cell Ranger ARC: Found 1 capture(s) across 1 batch(es)
[INFO] Cell Ranger ARC: Found 1 capture(s) across 1 batch(es)
[INFO] Batch aggregation: Found 1 ARC batch(es)
[INFO] Batch aggregation: Found 1 ARC batch(es)
[INFO] Doublet Detection: Using scrublet method
[INFO] Doublet Detection: Using scrublet method
host: midway3-login1.rcc.local
Building DAG of jobs...
SLURM run ID: d1ae475f-13fd-4b08-a17d-256664e8ad48
MinJobAge 120s (>= 120s). 'squeue' should work reliably for status queries.
Using shell: /usr/bin/bash
Provided remote nodes: 10
Job stats:
job                     count
--------------------  -------
aggregate_arc_batch         1
all                         1
cellranger_arc_aggr         1
cellranger_arc_count        1
create_arc_mudata           1
enrich_arc_metadata         1
run_scrublet                1
total                       7

Select jobs to execute...
Execute 1 jobs...

[Tue Mar  3 10:17:20 2026]
rule cellranger_arc_count:
    input: /project/lbarreiro/SHARED/PROGRAMS/refdata-cellranger-arc-GRCh38-2020-A-2.0.0, pbmc_unsorted_3k_library.csv
    output: 3K_PBMC_MULTIOME_PROCESSED/01_CELLRANGERARC_COUNT/1_A/outs/filtered_feature_bc_matrix.h5, 3K_PBMC_MULTIOME_PROCESSED/01_CELLRANGERARC_COUNT/1_A/outs/atac_fragments.tsv.gz, 3K_PBMC_MULTIOME_PROCESSED/01_CELLRANGERARC_COUNT/1_A/outs/gex_possorted_bam.bam, 3K_PBMC_MULTIOME_PROCESSED/01_CELLRANGERARC_COUNT/1_A/outs/filtered_feature_bc_matrix/barcodes.tsv.gz, 3K_PBMC_MULTIOME_PROCESSED/01_CELLRANGERARC_COUNT/1_A/outs/web_summary.html, 3K_PBMC_MULTIOME_PROCESSED/00_LOGS/1_A_arc_count.done
    log: 3K_PBMC_MULTIOME_PROCESSED/00_LOGS/1_A_arc_count.log
    jobid: 2
    reason: Missing output files: 3K_PBMC_MULTIOME_PROCESSED/01_CELLRANGERARC_COUNT/1_A/outs/filtered_feature_bc_matrix.h5, 3K_PBMC_MULTIOME_PROCESSED/00_LOGS/1_A_arc_count.done, 3K_PBMC_MULTIOME_PROCESSED/01_CELLRANGERARC_COUNT/1_A/outs/atac_fragments.tsv.gz; Set of input files has changed since last execution
    wildcards: batch=1, capture=A
    threads: 8
    resources: mem_mb=65536, mem_mib=62500, disk_mb=1000, disk_mib=954, tmpdir=<TBD>, slurm_account=pi-lbarreiro, slurm_partition=lbarreiro, runtime=720
Shell command: None
Job 2 has been submitted with SLURM jobid 46363794 (log: /project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/tests/PBMC_3K_MULTIOME/.snakemake/slurm_logs/rule_cellranger_arc_count/1_A/46363794.log).
```

Messages from this tool will always be prefaced in brackets e.g. `[INFO]`, `[WARNING]`, `[ERROR]`.

The first `[INFO]` prints the preprocessing steps enabled in the config file. In this tutorial, we enabled Cell Ranger ARC to process the [3K PBMC Epi Multiome ATAC + Gene Expression input data](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0) and Doublet detection with Scrublet:

```
[INFO] Config validated. Enabled steps: cellranger_arc, doublet_detection
```

Next, we print the Snakemake command running under the hood for convenient debugging:

```
[INFO] Running Snakemake with command: snakemake --snakefile /project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/cellranger_snakemake/workflows/main.smk --configfile pipeline_config.yaml --cores all --use-conda --profile HPC_profiles
```

After that, we print some more `[INFO]` about the run:

```
[INFO] ============================================================
[INFO] Single-Cell Preprocessing Pipeline
[INFO] ============================================================
[INFO] Project: 3K_PBMC_MULTIOME_PROCESSED
[INFO] Output directory: 3K_PBMC_MULTIOME_PROCESSED
[INFO] Enabled steps: cellranger_arc, doublet_detection
[INFO] ============================================================
```

Additionally, we always let you know with a `[WARNING]` if we change anything from your input files on the fly: 

```
[WARNING] Row 2 in 'libraries_list.tsv': Converting relative path 'pbmc_unsorted_3k_library.csv' to absolute path '/project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/tests/PBMC_3K_MULTIOME/pbmc_unsorted_3k_library.csv'
[WARNING] Row 2 in 'libraries_list.tsv': Converting relative path 'pbmc_unsorted_3k_library.csv' to absolute path '/project/lbarreiro/USERS/mschechter/github/cellranger-snakemake/tests/PBMC_3K_MULTIOME/pbmc_unsorted_3k_library.csv'
```

Finally, we print `[INFO]` from every entry job so you can fact check your workflow i.e. are these number of `batches` and `samples` you were expecting to preprocess?

```
[INFO] libraries_list.tsv file format is valid.
[INFO] libraries_list.tsv file format is valid.
[INFO] Cell Ranger ARC: Found 1 capture(s) across 1 batch(es)
[INFO] Cell Ranger ARC: Found 1 capture(s) across 1 batch(es)
[INFO] Batch aggregation: Found 1 ARC batch(es)
[INFO] Batch aggregation: Found 1 ARC batch(es)
[INFO] Doublet Detection: Using scrublet method
[INFO] Doublet Detection: Using scrublet method
```

Everything else are messages directly from Snakemake runing the workflow you configured! If you are new to Snakemake please take some time to orient yourself: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html



Log file paths for every Rule will be printing in the Snakemake stdout like this: 

```
log: 3K_PBMC_MULTIOME_PROCESSED/00_LOGS/1_A_arc_count.log
```

For example, you could explore the log for that Cell Ranger ARC job by printing the log file like this: 

```bash
$ cat 3K_PBMC_MULTIOME_PROCESSED/00_LOGS/1_A_arc_count.log
Martian Runtime - v4.0.5
Serving UI at http://midway3-0323.rcc.local:44865?auth=nTZ0qkEGdaSUvGtvxqySJq1uPiYiQHQQ5nIFFKYaxgM

Running preflight checks (please wait)...
Checking FASTQ folder...
Checking reference...
Checking reference_path (/project/lbarreiro/SHARED/PROGRAMS/refdata-cellranger-arc-GRCh38-2020-A-2.0.0) on midway3-0323.rcc.local...
Checking optional arguments...

...

2026-03-03 10:30:40 [runtime] (chunks_complete) ID.1_A.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._GEX_MATRIX_COMPUTER.MAKE_SHARD
2026-03-03 10:30:40 [runtime] (run:local)       ID.1_A.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._GEX_MATRIX_COMPUTER.MAKE_SHARD.fork0.join
2026-03-03 10:30:43 [runtime] (update)          ID.1_A.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._ATAC_MATRIX_COMPUTER.ALIGN_ATAC_READS.fork0 chunks running (0/9 completed)
```

## Examine the output directory structure

After successfully completing the workflow, you should see this resulting directory structure. Let's break it down: 

```bash
$ tree -L 2 3K_PBMC_MULTIOME_PROCESSED/
3K_PBMC_MULTIOME_PROCESSED/
├── 00_LOGS
│   ├── 1_A_arc_count.done
│   ├── 1_A_arc_count.log
│   ├── 1_A_arc_mudata.done
│   ├── 1_A_arc_mudata.log
│   ├── 1_arc_aggr.done
│   ├── 1_arc_aggr.log
│   ├── 1_arc_batch_aggregation.done
│   ├── 1_arc_batch_aggregation.log
│   ├── 1_arc_enrichment.done
│   ├── 1_arc_enrichment.log
│   ├── 1_A_scrublet.done
│   └── 1_A_scrublet.log
├── 01_CELLRANGERARC_COUNT
│   └── 1_A
├── 02_CELLRANGERARC_AGGR
│   └── 1_aggregation.csv
├── 03_ANNDATA
│   └── 1_A.h5mu
├── 04_BATCH_OBJECTS
│   └── 1_arc.h5mu
├── 06_DOUBLET_DETECTION
│   └── 1_A_scrublet.tsv.gz
└── 08_FINAL
    ├── 1_arc.h5mu
    ├── 1_arc_obs_summary.tsv.gz
    └── 1_arc_obs.tsv.gz
```

`00_LOGS`

This directory contains all the `.log` and `.done` files created throughout the workflow and are organized by `Batch_Capture_modality_rule`. The `.log` files will contain any STDOUT printed from every step of the workflow. This is allows you to dive in and interrogate any step your single-cell preprocessing. 

A quick way to find errors if you are debugging the workflow is to run:

```bash
grep -R "error" 3K_PBMC_MULTIOME_PROCESSED/00_LOGS
```

The `.done` files are an internal checklist to keep track of a subset of rules that finished (don't worry about it unless you are a developer and want to contribute to the code base).

`01_CELLRANGERARC_COUNT`

Here you will find all of the `Cell Ranger count` outputs for each individual capture.

`02_CELLRANGERARC_AGGR`

This will be the aggregated count matrices across batches. In this tutorial there is only one capture so you won't find any processed data here.

`03_ANNDATA`

Here you will find `MuData` objects for every capture. In this case it will be Muon because multiome.

`04_BATCH_OBJECTS`

`MuData` object from the aggregated object

`06_DOUBLET_DETECTION`

Doublet detection outputs from `Scrublet`

`08_FINAL`

`1_arc.h5mu`
`1_arc_obs_summary.tsv.gz`
`1_arc_obs.tsv.gz`

## Load the output for downstream analysis

### Examine barcode metadata

```python
# Preview cell metadata
mdata.obs[["batch_id", "capture_id", "cell_id", "scrublet_score", "predicted_doublet"]].head()
```

### Muon

The final MuData object in `08_FINAL/` contains both GEX and ATAC modalities with all preprocessing metadata joined in. Load it with:

```python 
import muon as mu

mdata = mu.read("3K_PBMC_MULTIOME_PROCESSED/08_FINAL/1_arc.h5mu")

# Verify traceability metadata is present and unique
assert "batch_id" in mdata.obs.columns
assert "capture_id" in mdata.obs.columns
assert "cell_id" in mdata.obs.columns
assert mdata.obs["cell_id"].is_unique, "cell_id must be unique!"

# Inspect the data
print(mdata)
print(f"\nTotal cells: {mdata.n_obs}")
print(f"\nModalities: {list(mdata.mod.keys())}")
print(f"\nGEX shape: {mdata['gex'].shape}")
print(f"\nATAC shape: {mdata['atac'].shape}")
print(f"\nObs columns: {list(mdata.obs.columns)}")
```

### Scanpy 

```python
# Access GEX modality for Scanpy downstream analysis
adata_gex = mdata["gex"]
```

### SnapATAC2

### Seurat

### ArchR

