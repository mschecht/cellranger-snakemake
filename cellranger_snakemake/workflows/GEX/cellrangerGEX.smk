import os
import sys
import shutil
import json
import pandas as pd

from collections import defaultdict
from cellranger_snakemake.utils import utils
from cellranger_snakemake.utils.custom_logger import custom_logger

# Snakemake workflow for Cell Ranger: https://www.10xgenomics.com/support/software/
#
# The libraries_list_GEX.tsv should have the following columns:
#  - batch: Batch number which groups captures together
#  - capture_ID: The capture_ID of the sample
#  - Sample: Prefix of the subset of the fastq files
#  - fastqs: Path to all fastq files

# Get directories
dirs_dict = utils.get_directories(config)

# Create directories if they don't exist
for dir_path in dirs_dict.values():
    os.makedirs(dir_path, exist_ok=True)

reference_genome = config["reference"]
libraries_file = config["libraries_GEX"]
jobmode = f"--jobmode={config['HPC_mode']}" if config.get("HPC_mode") else ""
chemistry = f"--chemistry={config['chemistry']}" if config.get("chemistry") else ""
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
    expected_columns = {"batch", "capture", "sample", "fastqs"}
    actual_columns = set(df.columns)
    if expected_columns != actual_columns:
        custom_logger.error(f"Expected columns {expected_columns}, found {actual_columns}")
        sys.exit(1)

    valid = True
    for idx, row in df.iterrows():
        if utils.has_underscore(str(row["batch"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Batch name '{row['batch']}' cannot contain underscores.")
            valid = False

        if utils.has_underscore(str(row["capture"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Capture name '{row['capture']}' cannot contain underscores.")
            valid = False
            
        fastq_path = row["fastqs"]
        error_location = f"Row {idx + 2} in '{filepath}'"
        if not isinstance(fastq_path, str):
            custom_logger.error(f"{error_location}: Path is not a string: {fastq_path}")
            valid = False

        # Try to identify if fastqs column has more than one path
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
        custom_logger.info(f"{libraries_file} file format is valid.")
        return df
    else:
        custom_logger.error(f"Some errors were found in {libraries_file}.")
        sys.exit(1)

df = sanity_check_libraries_list_GEX_tsv(libraries_file)

# Collect a summary
summary_dict = {}
capture_to_batch = defaultdict(list)
for idx, row in df.iterrows():
    batch = str(row["batch"])  # convert to string for consistency
    capture = row["capture"]
    sample = row["sample"]
    fastqs = row["fastqs"]
    out_dir = os.path.join(dirs_dict["CELLRANGERGEX_COUNT_DIR"], batch)

    # Initialize nested dicts if needed
    if batch not in summary_dict:
        summary_dict[batch] = {}

    if capture not in summary_dict[batch]:
        summary_dict[batch][capture] = {
            "sample": sample,
            "fastqs": set(),
            "output_dir": out_dir
        }

    # Add fastqs (assuming multiple fastqs per batch)
    summary_dict[batch][capture]["fastqs"].add(fastqs)

    # batch_to_samples remains a dict of sets
    capture_to_batch[batch].append(capture)

# Print summary to stderr to avoid interfering with DAG output
custom_logger.info(f"Found {sum(len(captures) for captures in capture_to_batch.values())} capture(s) across {len(summary_dict)} batch(es):")
for batch, captures in capture_to_batch.items():
    custom_logger.info(f"Batch {batch}: {len(captures)} capture(s)")
capture_to_batch_str = {capture: {str(b) for b in batch} for capture, batch in capture_to_batch.items()} # make sure batches are strings

# Set of files to be expected once all rules are finished
done_files = [
    os.path.join(dirs_dict["LOGS_DIR"], f"{batch}_organize.done")
    for batch in summary_dict
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
        done = os.path.join(dirs_dict["LOGS_DIR"], "{batch}_{capture}_cellranger_gex_count.done")
    params:
        sample = lambda wildcards: summary_dict[wildcards.batch][wildcards.capture]["sample"],
        fastqs = lambda wildcards: ",".join(list(summary_dict[wildcards.batch][wildcards.capture]["fastqs"]))
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{batch}_{capture}_cellranger_gex_count.log")
    threads: 8
    resources:
        mem_gb = 64
    run:
        shell(f"""
        cellranger count \
            --id={wildcards.batch}_{wildcards.capture} \
            --sample={params.sample} \
            --fastqs={params.fastqs} \
            --transcriptome={input.reference} \
            --chemistry={chemistry} \
            {jobmode} \
            >> {log} 2>&1

        touch {output.done}
        """)

rule cellranger_gex_aggr_csv:
    """
    Generate a CSV file for use with `cellranger aggr`, combining outputs from multiple batches per sample ID.
    """
    input:
        # Only fetch inputs for the current ID
        done_flags = lambda wildcards: [os.path.join(dirs_dict["LOGS_DIR"], f"{wildcards.ID}_{batch}_cellranger_gex_count.done") 
        for batch in summary_dict[wildcards.ID]]
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{ID}_cellranger_aggr_csv.log")
    output:
        aggr_csv = "{ID}_aggr.csv",
        done_flag = os.path.join(dirs_dict["LOGS_DIR"], "{ID}_aggr.done")
    run:
        sample = wildcards.ID
        aggr_rows = []

        for batch in batch_to_samples[sample]:
            row = {
                "sample_id": batch,
                "molecule_h5": os.path.abspath(f"{sample}_{batch}/outs/molecule_info.h5")}
            aggr_rows.append(row)

        # Write CSV
        aggr_csv_df = pd.DataFrame(aggr_rows)
        aggr_csv_df.to_csv(output.aggr_csv, index=False)

        shell(f"touch {output.done_flag}")
        shell(f'echo "Created aggregation CSV for ID {wildcards.ID} with {len(batch_to_samples[wildcards.ID])} batches" | tee -a {log}')

rule cellranger_aggr:
    """
    Run [cellranger aggr](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-aggr) 
    to aggregate multiple single-cell libraries.

    This rule only runs if there is more than one batch in a sample. If only one sample is present, the aggregation step is skipped.
    """
    log: 
        os.path.join(dirs_dict["LOGS_DIR"], "{ID}_cellranger_aggr.log")
    input:
        aggr_csv = rules.cellranger_gex_aggr_csv.output.aggr_csv
    output:
        done_flag = os.path.join(dirs_dict["LOGS_DIR"], "{ID}_cellranger_aggr.done")
    threads: 8
    resources:
        mem_gb = 64
    run:
        batch_samples = batch_to_samples[(wildcards.ID)]
        
        if len(batch_samples) > 1:
            shell(f"""
                cellranger aggr --id={{wildcards.ID}} \
                                --csv={{input.aggr_csv}} \
                                {jobmode} \
                                {normalize} \
                                >> {{log}} 2>&1

                touch {output.done_flag}     
                echo "Aggregated {len(batch_samples)} batches for sample {wildcards.ID}." | tee -a {log}
            """)

        else:
            # Create empty output directory to satisfy the rule
            with open("single_batch_sample.txt"), "w" as f:
                f.write(f"This sample contained only one batch: {batch_samples[wildcards.ID][0]}\n")
                f.write("Aggregation was skipped.\n")
            shell(f"touch {output.done_flag}")
            shell(f'echo "Sample {wildcards.ID} has only one batch ({batch_samples[wildcards.ID][0]}). Skipping cellranger aggr step." | tee -a {log}')

rule cellranger_gex_organize:
    """
    Organize all folders generated by cellranger.
    """
    input:
        done_flag = rules.cellranger_aggr.output.done_flag
    params:
        aggr_csv = rules.cellranger_gex_aggr_csv.output.aggr_csv,
        aggr_output_dir = dirs_dict["CELLRANGERGEX_AGGR_DIR"],
        count_outputs = lambda wildcards: [f"{wildcards.ID}_{batch}" for batch in summary_dict[wildcards.ID]],
        count_dir = lambda wildcards: [summary_dict[wildcards.ID][batch]["output_dir"] for batch in summary_dict[wildcards.ID]]
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{ID}_cellranger_organize.log")
    output:
        done_flag = os.path.join(dirs_dict["LOGS_DIR"], "{ID}_organize.done")
    run:
        # Create cellranger output subdirectory for each sample ID
        for d in params.count_dir:
            os.makedirs(d, exist_ok=True)

        # Move cellranger count directories into their respective final output directory
        for src, dest in zip(params.count_outputs, params.count_dir):
            shutil.move(src, dest)

        # Move cellranger aggr directories into final output directory
        shutil.move(params.aggr_csv, wildcards.ID)
        shutil.move(wildcards.ID, params.aggr_output_dir)

        shell(f"touch {output.done_flag}")