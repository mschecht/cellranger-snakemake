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
    out_dir = os.path.join(dirs_dict["CELLRANGERGEX_COUNT_DIR"], ID, batch)

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
batch_to_samples_str = {sample: {str(b) for b in batches} for sample, batches in batch_to_samples.items()}

# Create output directories
for sample, batches_dict in summary_dict.items(): 
    for batch, batch_info in batches_dict.items():
        os.makedirs(batch_info["output_dir"], exist_ok=True)

# Set of files to be expected once all rules are finished
done_files = [
    os.path.join(dirs_dict["LOGS_DIR"], f"{ID}_{batch}_cellranger_gex_count.done")
    for ID in summary_dict
    for batch in summary_dict[ID]
]

rule all:
    input:
        done_files

rule cellranger_gex_count:
    input:
        reference = reference_genome
    output:
        done = os.path.join(dirs_dict["LOGS_DIR"], "{ID}_{batch}_cellranger_gex_count.done")
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
            --id={wildcards.ID} \
            --sample={params.sample} \
            --fastqs={params.fastqs} \
            --transcriptome={input.reference} \
            {jobmode} \
            >> {log} 2>&1
        mv {wildcards.ID} {params.outdir}
        touch {output.done}
        """)