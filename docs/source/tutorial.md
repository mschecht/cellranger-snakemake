# Tutorial and test case for development for GEX, ATAC, and ARC workflow

In this section, we will run `cellranger-snakemake` using a test dataset from Cell Ranger (derived from fasta files are used by the internal testing tool for cellranger called [cellranger testrun](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in#testrun)) to get new users ready to go as well as developers who need test cases for each of the workflow modes. Follow the steps to set up the test dataset, and run basic commands. 

> 📌 **Note**: This section works the same for either GEX or ATAC and has a few modifications for ARC which we note below.

## 1. Launch the conda environment

If you haven't installed `cellranger-snakemake` already, please refer to the [installation instructions](https://github.com/mschecht/cellranger-snakemake?tab=readme-ov-file#installation-instructions) on the main page or see our [installation instructions](installation.md).

Activate the `cellranger-snakemake` conda environment with the following command:

```bash
conda activate snakemake8
```

## 2. Explore the help menus

Here is how you can check out the help menu for all positional arugments:

```bash
# Read about positional arguments
snakemake-run-cellranger --help
```

To learn more about a specific positional argument, include the argument and `--help` like this:

```bash
# Help menu for run positional argument
snakemake-run-cellranger run --help
```

## 3. Generate input files for test

The command `snakemake-run-cellranger generate-test-data` conveniently creates a directory containing all the input files necessary you need to run the test dataset. 

```bash
# Read about test data set
snakemake-run-cellranger generate-test-data -h

snakemake-run-cellranger generate-test-data GEX --output-dir tests/00_TEST_DATA_GEX
snakemake-run-cellranger generate-test-data ATAC --output-dir tests/00_TEST_DATA_ATAC
snakemake-run-cellranger generate-test-data ARC --output-dir tests/00_TEST_DATA_ARC
```

This should have produced the following file structure:

```bash
$ tree tests/00_TEST_DATA_GEX
tests/00_TEST_DATA_GEX
├── HPC_profiles
│   └── config.yaml
├── libraries_list_gex.tsv
├── reference_gex.txt
└── test_config_gex.yaml
```


## 4. Input files

Let's walk through the input files necessary to run the workflow!

### `config.yaml`

This YAML file contains all the bells and whistles needed to run the underlying snakemake workflow!

To generate the config files for a workflow customized to your data run this command below. It will interactively ask you which steps you plan on running and automatically produce a config file.

```bash
snakemake-run-cellranger init-config
```

You can also run this command to generate a default config yaml with every configuration available:

```bash
snakemake-run-cellranger init-config --get-default-config
```

For this tutorial, here is the test config yaml file: 

```bash
$ cat tests/00_TEST_DATA_GEX/test_config_gex.yaml
project_name: test_gex
output_dir: test_output_gex
resources:
  mem_gb: 64
  tmpdir: ''
directories_suffix: none
cellranger_gex:
  enabled: true
  reference: /path/to/cellranger-9.0.1/external/cellranger_tiny_ref
  libraries: tests/00_TEST_DATA_GEX/libraries_list_gex.tsv
  chemistry: auto
  normalize: none
  create-bam: false
  threads: 10
  mem_gb: 64
demultiplexing:
  enabled: true
  method: vireo
  vireo:
    donors: 2
    cellsnp:
      vcf: /path/to/vcf/file.vcf.gz
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
doublet_detection:
  enabled: false
  method: scrublet
  scrublet:
    expected_doublet_rate: 0.06
    min_counts: 2
    min_cells: 3
celltype_annotation:
  enabled: false
  method: celltypist
  celltypist:
    model: Immune_All_Low.pkl
    majority_voting: false
```

### `libraries_list.tsv`

This input file is a TSV file that contains the metadata and paths for your cellranger libraries. Here is the format:

| batch | capture | sample  | fastqs                   |
|-------|---------|---------|--------------------------|
| A     | 1       | ABC-A-1 | path/to/data/GEX/fastqs/ |
| A     | 2       | IJK-A-2 | path/to/data/GEX/fastqs/ |
| B     | 1       | XYZ-A-1 | path/to/data/GEX/fastqs/ |

Column descriptions:
- `batch`: batch ID for grouping captures
- `capture`: capture identifier or lane on the 10X chip
- `sample`: prefix of the filenames of FASTQs to select
- `fastqs`: full path(s) to where the input FASTQ files are located - if providing multiple paths, separate them with commas.


**Note**: For the ARC workflow, the input file is a little bit different. You will need to create a tab-separated file that contains the metadata and paths for to cellranger ARC [library csv files](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/single-library-analysis#create-a-libraries-csv-file) (files that contain paths the ATAC and GEX FASTQ files). This file, which we will call `libraries_list_ARC.tsv` during this tutorial, follows the following format:

| batch | capture | CSV                                |
|-------|---------|------------------------------------|
| A     | 1       | path/to/data/ATAC/ARC_library.csv/ |
| A     | 2       | path/to/data/ATAC/ARC_library.csv/ |
| A     | 3       | path/to/data/ATAC/ARC_library.csv/ |

Column descriptions:
- `batch`: batch ID for grouping captures
- `capture`: capture identifier or lane on the 10X chip
- `CSV`: path to ARC library CSV (contains paths to fastas for both GEX and ATAC)


### `HPC_profiles/`

The `HPC_profiles/` directory contains another `config.yaml` that configures the cloud computing and HPC infrastructure settings to help `snakemake` launch parallel jobs. This config would be the argument for `snakemake --profile HPC_profiles`. You can read more about it [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). See [section 7](#7-launching-on-hpc) for detailed usage.

For this test dataset, we made the default HPC profile config to be compatible with [SLURM](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html). However, you can [install another executor](https://snakemake.github.io/snakemake-plugin-catalog/index.html) to match you local HPC/cloud computing infrastructure. 

```bash
$ cat tests/00_TEST_DATA_GEX/HPC_profiles/config.yaml
executor: slurm
jobs: 10
default-resources:
- slurm_account=pi-lbarreiro
- slurm_partition=lbarreiro-hm
- runtime=720
retries: 2
latency-wait: 60
printshellcmds: true
keep-going: true
rerun-incomplete: true
```

## 5. Run a dry run

Before you run the workflow it's a good idea to see how many jobs will be run to make sure your input files contain all the paths.

```bash
# Read about this command
snakemake-run-cellranger run -h

# Dry run
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dry-run
```

You can also visualize this with a [dag file](https://en.wikipedia.org/wiki/Directed_acyclic_graph):

```bash
# Generate workflow DAG
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dag | dot -Tpng > dag.png
```

## 6. Run the tool!

```bash
# Remove previous test runs
rm -rf 1_L00*
rm -r test_output_gex

# Local execution
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1
```

The flag `--snakemake-args` passes and arguments after it directly to `snakemake`. Please note that this flag has to be the very last flag in the command:

```bash
# Local execution - add more arguments to snakemake
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --snakemake-args --jobs 2
```

(7-launching-on-hpc)=
## 7. Launching on HPC

To launch on the HPC, we will use the `--snakemake-args` command to pass additional arguments to snakemake to let it know we are going to use an HPC. The `--snakemake-args` must be the LAST argument and anything after it will be snakemake arguments passed directly to `snakemake`.

> **Note**: If the directory gets locked, you can unlock it by running:
> `snakemake-run-cellranger run --config-file <your_config.yaml> --cores 1 --snakemake-args --unlock`

The argument we will be passing straight to `snakemake` will be `--profile`. The provides `snakemake` with a path to a configurgation file that contains parameters fro runnign the is workflow on an HPC or cloud computing environment. Run `snakemake -h` to read more detail.

The command `snakemake-run-cellranger generate-test-data` you ran above already produced a boiler plate config yaml file filled out for SLURM here:

```bash
$ cat tests/00_TEST_DATA_GEX/profiles/config.yaml
executor: 'slurm'
jobs: 1
default-resources:
- slurm_account=
- slurm_partition=
- mem_mb=
- runtime=
retries: 2
latency-wait: 60
```

You read about HPC executor functionality [here](https://snakemake.github.io/snakemake-plugin-catalog/). Fill out this config with HPC/cloud computing info that works for you! We made autogenerated an example for `Slurm`.

What is the difference between `--cores` and `--jobs`? The `--cores` command assigns the number of CPUs per jobs while the `--jobs` argument controls how many parallel jobs can be run at the same time.

```bash
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --snakemake-args --unlock

# HPC execution - `--cores all` tell snakemake to use the `threads` assigned to each rule.
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml \
                             --cores all \
                             --snakemake-args --profile tests/00_TEST_DATA_GEX/HPC_profiles
```