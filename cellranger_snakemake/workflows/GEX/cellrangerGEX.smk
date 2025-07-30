import os
import sys
import json
import pandas as pd

from collections import defaultdict
from cellranger_snakemake.utils import utils
from cellranger_snakemake.utils.custom_logger import custom_logger

# Snakemake workflow for Cell Ranger: https://www.10xgenomics.com/support/software/
#
# The libraries_list_GEX.tsv should have the following columns:
#  - batch: Batch number which groups samples together
#  - ID: The ID of the sample
#  - Sample: Prefix of the subset of the fastq files
#  - Fastqs: Path to all fastq files

# Get directories
dirs_dict = utils.get_directories(config)

# Create directories if they don't exist
for dir_path in dirs_dict.values():
    os.makedirs(dir_path, exist_ok=True)

reference_genome = config["reference"]
libraries_file = config["libraries_GEX"]
jobmode = f"--jobmode={config['HPC_mode']}" if config.get("HPC_mode") else ""
mempercore = f"--mempercore={config['mempercore']}" if config.get("mempercore") else ""
normalize = f"--normalize={config['normalize']}" if config.get("normalize") else ""

if mempercore and not jobmode:
    raise ValueError("You need to set the jobmode in the config file if you want to use the memory per core option.")

def sanity_check_libraries_list_GEX_tsv(filepath, log_file=None):
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
    expected_columns = {"ID", "batch", "sample", "Fastqs"}
    actual_columns = set(df.columns)
    if expected_columns != actual_columns:
        custom_logger.error(f"Expected columns {expected_columns}, found {actual_columns}")
        sys.exit(1)

    valid = True
    for idx, row in df.iterrows():
        fastq_path = row["Fastqs"]
        error_location = f"Row {idx + 2} in '{filepath}'"
        if not isinstance(fastq_path, str):
            custom_logger.error(f"{error_location}: Path is not a string: {fastq_path}")
            valid = False

        # Try to identify if Fastqs column has more than one path
        DELIMITERS = [",", ";", ":", " "]
        # First, check for disallowed delimiters
        for delim in DELIMITERS[1:]:  # skip comma (allowed)
            if delim in fastq_path:
                custom_logger.error(f"{error_location}: Path contains invalid delimiter '{delim}'. Only comma-separated paths are allowed.")
                sys.exit(1)

        # Split the paths by comma
        suspected_paths = [path.strip() for path in fastq_path.split(",") if path.strip()]
        # Validate each individual path
        for path in suspected_paths:
            if not os.path.isabs(path):
                custom_logger.error(f"{error_location}: Path is not absolute: {path}")
                valid = False
            elif not os.path.exists(path):
                custom_logger.error(f"{error_location}: Path does not exist: {path}")
                valid = False

    if valid:
        custom_logger.info("libraries_list.tsv file format is valid.")
        return df
    else:
        custom_logger.error("Some errors were found in the file. Please check the log.")
        sys.exit(1)

df = sanity_check_libraries_list_GEX_tsv(libraries_file)

# Collect a summary
summary_dict = {}
batch_to_samples = defaultdict(list)
for idx, row in df.iterrows():
    ID = row["ID"]
    batch = str(row["batch"])  # convert to string for consistency
    sample = row["sample"]
    fastqs = row["Fastqs"]
    out_dir = os.path.join(dirs_dict["CELLRANGERGEX_COUNT_DIR"], ID)

    # Initialize nested dicts if needed
    if ID not in summary_dict:
        summary_dict[ID] = {}

    if batch not in summary_dict[ID]:
        summary_dict[ID][batch] = {
            "sample": sample,
            "fastqs": set(),
            "output_dir": out_dir
        }

    # Add fastqs (assuming multiple fastqs per batch)
    summary_dict[ID][batch]["fastqs"].add(fastqs)

    # batch_to_samples remains a dict of sets
    batch_to_samples[ID].append(batch)

# Print summary to stderr to avoid interfering with DAG output
custom_logger.info(f"Found {sum(len(batches) for batches in batch_to_samples.values())} batch(es) across {len(summary_dict)} sample(s):")
for batch, samples in batch_to_samples.items():
    custom_logger.info(f"Sample {batch}: {len(samples)} batch(es)")
batch_to_samples_str = {sample: {str(b) for b in batches} for sample, batches in batch_to_samples.items()} # make sure batches are strings

# Set of files to be expected once all rules are finished
done_files = [
    os.path.join(dirs_dict["LOGS_DIR"], f"{ID}_aggr.done")
    for ID in summary_dict
]

rule all:
    input:
        done_files

rule cellranger_gex_count:
    """
    [cellranger count](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)

    This program counts gene expression (targeted or whole-transcriptome) and/or feature barcode reads from a single sample and GEM well.
    """
    input:
        reference = reference_genome
    output:
        done = os.path.join(dirs_dict["LOGS_DIR"], "{ID}_{batch}_cellranger_gex_count.done"),
        molecule_info = os.path.join(dirs_dict["CELLRANGERGEX_COUNT_DIR"], "{ID}/{ID}_{batch}/outs/molecule_info.h5")
    params:
        sample = lambda wildcards: summary_dict[wildcards.ID][wildcards.batch]["sample"],
        fastqs = lambda wildcards: list(summary_dict[wildcards.ID][wildcards.batch]["fastqs"])[0],
        outdir = lambda wildcards: summary_dict[wildcards.ID][wildcards.batch]["output_dir"]
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{ID}_{batch}_cellranger_gex_count.log")
    threads: 8
    resources:
        mem_gb = 64
    run:
        shell(f"""
        cellranger count \
            --id={wildcards.ID}_{wildcards.batch} \
            --sample={params.sample} \
            --fastqs={params.fastqs} \
            --transcriptome={input.reference} \
            {jobmode} \
            >> {log} 2>&1

        mkdir -p {params.outdir}
        mv {wildcards.ID}_{wildcards.batch} {params.outdir}
        touch {output.done}
        """)

rule cellranger_gex_aggr_csv:
    """
    Generate a CSV file for use with `cellranger aggr`, combining outputs from multiple batches per sample ID.
    """
    input:
        # Only fetch inputs for the current ID
        done_flags = lambda wildcards: [
            os.path.join(dirs_dict["LOGS_DIR"], f"{wildcards.ID}_{batch}_cellranger_gex_count.done")
            for batch in summary_dict[wildcards.ID]
        ],

        molecule_info = lambda wildcards: [
            os.path.join(
                dirs_dict["CELLRANGERGEX_COUNT_DIR"],
                f"{wildcards.ID}/{wildcards.ID}_{batch}/outs/molecule_info.h5"
            )
            for batch in summary_dict[wildcards.ID]
        ]
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{ID}_cellranger_aggr_csv.log")
    output:
        aggr_csv = os.path.join(dirs_dict["CELLRANGERGEX_AGGR_DIR"], "{ID}/{ID}_aggr.csv"),
        done_flag = os.path.join(dirs_dict["LOGS_DIR"], "{ID}_aggr.done")
    run:
        sample = wildcards.ID
        aggr_rows = []

        for batch in batch_to_samples[sample]:
            row = {
                "sample_id": batch,
                "molecule_h5": os.path.abspath(os.path.join(summary_dict[sample][batch]["output_dir"], f"{sample}_{batch}/outs/molecule_info.h5"))
            }
            aggr_rows.append(row)

        # Create output directory
        os.makedirs(os.path.dirname(output.aggr_csv), exist_ok=True)
        
        # Write CSV
        aggr_csv_df = pd.DataFrame(aggr_rows)
        aggr_csv_df.to_csv(output.aggr_csv, index=False)

        shell(f"touch {output.done_flag}")
        shell(f'echo "Created aggregation CSV for ID {wildcards.ID} with {len(batch_to_samples[wildcards.ID])} batches" | tee -a {log}')
