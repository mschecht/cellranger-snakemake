import os
import sys
import json
import pandas as pd

from collections import defaultdict
from cellranger_snakemake.utils import utils
from cellranger_snakemake.utils.custom_logger import custom_logger

# Snakemake workflow for Cell Ranger Arc: https://www.10xgenomics.com/support/software/cell-ranger-arc/latest
# 
# Here is an example of the command:
# snakemake --snakefile cellrangerarc.smk --configfile config.json --cores 1
#
# You will need a config file that contains paths to required files
# The config.json should look like this:
# {
#  "reference": "path/to/reference_genome_directory",
#  "libraries": "path/to/libraries_list.tsv",
#  "HPC_mode": "local",  # optional: local, slurm, etc.
#  "mempercore": 8,      # optional: memory per core in GB
#  "normalize": "depth"  # optional: normalization method
# }
# The libraries_list.tsv should have the following columns:
#  - batch: Batch number which groups samples together
#  - sample_ID: The ID of the sample
#  - CSV: The path to the CSV file containing the library information

# Get directories
dirs_dict = utils.get_directories(config)

# Create directories if they don't exist
for dir_path in dirs_dict.values():
    os.makedirs(dir_path, exist_ok=True)

reference_genome = config["reference"]
libraries_file = config["libraries"]
jobmode = f"--jobmode={config['HPC_mode']}" if config.get("HPC_mode") else ""
mempercore = f"--mempercore={config['mempercore']}" if config.get("mempercore") else ""
normalize = f"--normalize={config['normalize']}" if config.get("normalize") else ""

if mempercore and not jobmode:
    raise ValueError("You need to set the jobmode in the config file if you want to use the memory per core option.")

def sanity_check_libraries_list_tsv(filepath, log_file=None):
    """
    Check the format of the libraries list TSV file.
    Args:
        filepath (str): Path to the TSV file.
        log_file (str, optional): Path to the log file. Defaults to None.
    Returns:
        pd.DataFrame: Valid dataframe or exits on error
    """

    # Read the file
    try:
        df = pd.read_csv(filepath, sep="\t")
        df.columns = df.columns.str.strip()
    except Exception as e:
        custom_logger.error(f"Pandas could not read your libraries filepath here: '{filepath}'. "
                     f"This was the error: {e}")
        sys.exit(1)

    # Validate headers
    expected_columns = {"batch", "sample_ID", "CSV"}
    actual_columns = set(df.columns)
    if expected_columns != actual_columns:
        custom_logger.error(f"Expected columns {expected_columns}, found {actual_columns}")
        sys.exit(1)

    valid = True
    for idx, row in df.iterrows():
        csv_path = row["CSV"]
        error_location = f"Row {idx + 2} in '{filepath}'"
        if not isinstance(csv_path, str):
            custom_logger.error(f"{error_location}: CSV path is not a string: {csv_path}")
            valid = False
        elif not os.path.isabs(csv_path):
            custom_logger.error(f"{error_location}: Path is not absolute: {csv_path}")
            valid = False
        elif not csv_path.endswith(".csv"):
            custom_logger.error(f"{error_location}: The path '{csv_path}' does not end with '.csv'. Please ensure all CSV paths are correctly specified and have the proper file extension.")
            valid = False
        elif not os.path.exists(csv_path):
            custom_logger.error(f"{error_location}: CSV file does not exist: {csv_path}")
            valid = False

    if valid:
        custom_logger.info("libraries_list.tsv file format is valid.")
        return df
    else:
        custom_logger.error("Some errors were found in the file. Please check the log.")
        sys.exit(1)

df = sanity_check_libraries_list_tsv(libraries_file)

# Collect a summary and validate library CSV files
summary_dict = {}
batch_to_samples = defaultdict(list)

for idx, row in df.iterrows():
    batch = row["batch"]
    # FIXME: change to "capture_ID" to be consistent with 10x Genomics terminology
    sample_id = row["sample_ID"] 
    csv_path = row["CSV"]
    
    # Read and validate the library CSV
    try:
        lib_df = pd.read_csv(csv_path)
        required_lib_columns = {"fastqs", "sample", "library_type"}
        if not required_lib_columns.issubset(lib_df.columns):
            raise ValueError(f"Library CSV missing required columns: {required_lib_columns - set(lib_df.columns)}")
        
        # Extract paths for ATAC and GEX
        atac_rows = lib_df[lib_df["library_type"] == "Chromatin Accessibility"]
        gex_rows = lib_df[lib_df["library_type"] == "Gene Expression"]
        
        if len(atac_rows) == 0:
            raise ValueError("No 'Chromatin Accessibility' library found")
        if len(gex_rows) == 0:
            raise ValueError("No 'Gene Expression' library found")
            
        atac_path = atac_rows["fastqs"].iloc[0]
        gex_path = gex_rows["fastqs"].iloc[0]
        
    except Exception as e:
        custom_logger.error(f"Error reading library CSV for sample {sample_id}: {e}")
        sys.exit(1)
    
    out_dir = os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], sample_id, "outs")

    summary_dict[sample_id] = {
        "batch": batch,
        "Library_path": csv_path,
        "ATAC_path": atac_path,
        "GEX_path": gex_path,
        "Output_dir": out_dir
    }

    batch_to_samples[batch].append(sample_id)

# Print summary to stderr to avoid interfering with DAG output

custom_logger.info(f"Found {len(summary_dict)} samples across {len(batch_to_samples)} batches:")
for batch, samples in batch_to_samples.items():
    custom_logger.info(f"Batch {batch}: {len(samples)} samples")

rule all:
    input:
        expand(
            os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.done"),
            batch=batch_to_samples.keys()
        )

rule cellranger_arc_count:
    """
    [cellranger-arc count](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/command-line-arguments#count)

    This program counts ATAC and gene expression reads from a single library against a reference genome.
    """
    input:
        reference = reference_genome,
        library_csv = lambda wc: summary_dict[wc.sample_id]["Library_path"]
    output:
        done_flag = touch(os.path.join(dirs_dict["LOGS_DIR"], "{sample_id}_cellranger_arc_count.done")),
        atac_fragments = os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], "{sample_id}", "outs", "atac_fragments.tsv.gz"),
        per_barcode_metrics = os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], "{sample_id}", "outs", "per_barcode_metrics.csv"),
        gex_molecule_info = os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], "{sample_id}", "outs", "gex_molecule_info.h5")
    params:
        count_dir = dirs_dict["CELLRANGERARC_COUNT_DIR"],
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{sample_id}_cellranger_arc_count.log")
    threads: 8
    resources:
        mem_gb = 64
    run:
        Library_path = summary_dict[wildcards.sample_id]["Library_path"]
        os.makedirs(params.count_dir, exist_ok=True)

        shell(f"""
            cellranger-arc count --id={wildcards.sample_id} \
                                 --reference={reference_genome} \
                                 --libraries={Library_path} \
                                 {jobmode} \
                                 {mempercore} \
                                 >> {log} 2>&1
            mkdir -p {params.count_dir}
            mv {wildcards.sample_id} {params.count_dir}
            """)

rule cellranger_arc_aggr_csv:
    """
    Generate a CSV file for use with `cellranger-arc aggr`, combining outputs from multiple samples.

    This rule collects the following files from each `cellranger-arc count` sample directory:
        - `outs/atac_fragments.tsv.gz`
        - `outs/per_barcode_metrics.csv`
        - `outs/gex_molecule_info.h5`

    It assembles a CSV with one row per sample, suitable as input to:
        `cellranger-arc aggr --csv`.

    Output CSV columns:
        - library_id
        - atac_fragments
        - per_barcode_metrics
        - gex_molecule_info
    """
    log: os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr_csv.log")
    input:
        # Require counts for all samples in batch before aggregation
        done_flags = lambda wildcards: [
            os.path.join(dirs_dict["LOGS_DIR"], f"{sample}_cellranger_arc_count.done")
            for sample in batch_to_samples[int(wildcards.batch)]
        ],
        # Also require the actual output files
        atac_fragments = lambda wildcards: [
            os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], sample, "outs", "atac_fragments.tsv.gz")
            for sample in batch_to_samples[int(wildcards.batch)]
        ],
        per_barcode_metrics = lambda wildcards: [
            os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], sample, "outs", "per_barcode_metrics.csv")
            for sample in batch_to_samples[int(wildcards.batch)]
        ],
        gex_molecule_info = lambda wildcards: [
            os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], sample, "outs", "gex_molecule_info.h5")
            for sample in batch_to_samples[int(wildcards.batch)]
        ]
    output:
        aggr_csv = os.path.join(dirs_dict["CELLRANGERARC_AGGR_DIR"], "{batch}", "{batch}_aggr.csv")
    run:
        batch = int(wildcards.batch)
        aggr_rows = []

        for sample in batch_to_samples[batch]:
            row = {
                "library_id": sample,
                "atac_fragments": os.path.abspath(os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], sample, "outs", "atac_fragments.tsv.gz")),
                "per_barcode_metrics": os.path.abspath(os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], sample, "outs", "per_barcode_metrics.csv")),
                "gex_molecule_info": os.path.abspath(os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], sample, "outs", "gex_molecule_info.h5"))
            }
            aggr_rows.append(row)

        # Create output directory
        os.makedirs(os.path.dirname(output.aggr_csv), exist_ok=True)
        
        aggr_csv_df = pd.DataFrame(aggr_rows)
        aggr_csv_df.to_csv(output.aggr_csv, index=False)
        
        print(f"Created aggregation CSV for batch {batch} with {len(aggr_rows)} samples", file=sys.stderr)

rule cellranger_arc_aggr:
    """
    Run [cellranger-arc aggr](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/command-line-arguments) 
    to aggregate multiple single-cell multiome libraries.

    This rule only runs if there is more than one sample in a batch. If only one sample is present,
    the aggregation step is skipped.
    """
    log: os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.log")
    input:
        aggr_csv = rules.cellranger_arc_aggr_csv.output.aggr_csv
    output:
        done_flag = touch(os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.done")),
        # Define expected outputs from aggregation
        aggr_output = directory(os.path.join(dirs_dict["CELLRANGERARC_AGGR_DIR"], "{batch}", "outs"))
    threads: 8
    resources:
        mem_gb = 64
    run:
        batch_samples = batch_to_samples[int(wildcards.batch)]
        
        if len(batch_samples) > 1:
            # Change to aggregation directory
            ID = os.path.join(dirs_dict['CELLRANGERARC_COUNT_DIR'], wildcards.batch) 
            shell(f"""
                cellranger-arc aggr --id={ID} \
                                    --reference={reference_genome} \
                                    --csv={{input.aggr_csv}} \
                                    {jobmode} \
                                    {normalize} \
                                    >> {{log}} 2>&1
            """)
            print(f"Aggregated {len(batch_samples)} samples for batch {wildcards.batch}", file=sys.stderr)
        else:
            print(f"Batch {wildcards.batch} has only one sample ({batch_samples[0]}). Skipping cellranger-arc aggr step.", file=sys.stderr)
            # Create empty output directory to satisfy the rule
            os.makedirs(output.aggr_output, exist_ok=True)
            with open(os.path.join(output.aggr_output, "single_sample_batch.txt"), "w") as f:
                f.write(f"This batch contained only one sample: {batch_samples[0]}\n")
                f.write("Aggregation was skipped.\n")