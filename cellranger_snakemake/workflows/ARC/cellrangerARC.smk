import os
import json
import logging
import pandas as pd

from collections import defaultdict

# Snakmake workflow to Cell Ranger Arc: https://www.10xgenomics.com/support/software/cell-ranger-arc/latest
# 
# Here is an example of the command:
# snakemake --snakefile cellrangerarc.smk --configfile config.json --cores 1
#
# You will need a config file that contains paths to anvi'o files such as external-genomes.txt and metagenomes.txt
# The config.json should look like this:
# {
#  "reference": "path/to/reference_genome_directory",
#  "libraries": "path/to/libraries_list.tsv",
# }
# The libraries_list.tsv should have the following columns:
#  - batch: Batch number which groups samples together
#  - sample_ID: The ID of the sample
#  - CSV: The path to the CSV file containing the library information

dirs_dict = {"LOGS_DIR": "00_LOGS",
            "CELLRANGERARC_COUNT_DIR": "01_CELLRANGERARC_COUNT",
            "CELLRANGERARC_AGGR_DIR": "02_CELLRANGERARC_AGGR"}  

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
    """

    # Set up logging
    logger = logging.getLogger("LibrariesListChecker")
    logger.setLevel(logging.DEBUG)
    logger.handlers = []  # avoid duplicate logs in Jupyter or re-runs

    # Log to console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    logger.addHandler(console_handler)

    # Optional: log to file
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    # Read the file
    try:
        df = pd.read_csv(filepath, sep="\t")
        df.columns = df.columns.str.strip()
    except Exception as e:
        logger.error(f"Pandas could not read your libraries filepath here: '{filepath}'",
                     f"This was the error: {e}")
        return False

    # Validate headers
    expected_columns = {"batch", "sample_ID", "CSV"}
    actual_columns = set(df.columns)
    if expected_columns != actual_columns:
        logger.error(f"Expected columns {expected_columns}, found {actual_columns}")
        return False

    valid = True
    for idx, row in df.iterrows():
        csv_path = row["CSV"]
        if not isinstance(csv_path, str):
            logger.error(f"Row {idx + 2}: CSV path is not a string: {csv_path}")
            valid = False
        elif not os.path.isabs(csv_path):
            logger.error(f"Row {idx + 2}: Path is not absolute: {csv_path}")
            valid = False
        elif not csv_path.endswith(".csv"):
            logger.error(f"Row {idx + 2}: Path does not end with .csv: {csv_path}")
            valid = False

    if valid:
        logger.info("File format is valid.")
        return df
    else:
        logger.error(f"Some errors were found in the file. Please check the log.")
        sys.exit()

df = sanity_check_libraries_list_tsv(libraries_file)

# Collect a summary
summary_dict = {}
batch_to_samples = defaultdict(list)

for idx, row in df.iterrows():
    batch = row["batch"]
    sample_id = row["sample_ID"]
    csv_path = row["CSV"]
    df = pd.read_csv(csv_path)
    atac_path = df.loc[df["library_type"] == "Chromatin Accessibility", "fastqs"].values[0]
    gex_path = df.loc[df["library_type"] == "Gene Expression", "fastqs"].values[0]
    out_dir = os.path.join(sample_id, "outs") 

    summary_dict[sample_id] = {
        "batch": batch,
        "Library_path": csv_path,
        "ATAC_path": atac_path,
        "GEX_path": gex_path,
        "Output_dir": out_dir
    }

    batch_to_samples[batch].append(sample_id)

rule all:
    input:
        expand(
            os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.done"),
            batch=batch_to_samples.keys()
        )


# rule cellranger_arc_mkfastq:
#     """
#     The UChicago sequencing core facility does this step.
#     """
#     log: os.path.join(dirs_dict["LOGS_DIR"], "{sample_id}_cellranger_arc_mkfastq.log")
#     input:
#         atac_path = lambda wildcards: summary_dict[wildcards.sample_id]["ATAC_path"],
#         gex_path = lambda wildcards: summary_dict[wildcards.sample_id]["GEX_path"]
#     output:
#         touch(os.path.join(dirs_dict["LOGS_DIR"], "{sample_id}_cellranger_arc_mkfastq.done"))
#     run:
#         atac_path = input.atac_path
#         gex_path = input.gex_path

#         shell("cellranger-arc mkfastq --id={sample_id} \
#                                         --run={atac_path} \
#                                         --output-dir={gex_path} \
#                                         --localcores=12 \
#                                         --localmem=400")

rule cellranger_arc_count:
    """
    [cellranger-arc count](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/command-line-arguments#count)

    This program is used to count ATAC and gene expression reads from a single library against a reference genome.
    """
    log: os.path.join(dirs_dict["LOGS_DIR"], "{sample_id}_cellranger_arc.log")
    input:  
    output: touch(os.path.join(dirs_dict["LOGS_DIR"], "{sample_id}_cellranger_arc_count.done"))
    run:
        Library_path = summary_dict[wildcards.sample_id]["Library_path"]
        
        shell("""
                cellranger-arc count --id={wildcards.sample_id} \
                                     --reference={reference_genome} \
                                     --libraries={Library_path} \
                                     {jobmode} \
                                     {mempercore} \
                                     >> {log} 2>&1
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

    Wildcards:
        - {batch}: Batch identifier, used to select samples and name the output

    Assumes the following are available globally:
        - dirs_dict["LOGS_DIR"]: Path to `.done` flags confirming count completion
        - batch_to_samples: Dictionary mapping batch IDs to sample IDs
        - summary_dict: Dictionary of sample metadata (not used directly here)
    """
    log: os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr_csv.log")
    input:
        done_flags=lambda wildcards: [
            os.path.join(dirs_dict["LOGS_DIR"], "{sample}_cellranger_arc_count.done")
            for sample in batch_to_samples[wildcards.batch]
        ]
    output:
        aggr_csv = os.path.join("{batch}", "{batch}_aggr.csv")
    run:
        batch = int(wildcards.batch)
        aggr_rows = []

        for sample in batch_to_samples[batch]:
            row = {
                "library_id": sample,
                "atac_fragments": os.path.abspath(os.path.join(sample, "outs", "atac_fragments.tsv.gz")),
                "per_barcode_metrics": os.path.abspath(os.path.join(sample, "outs", "per_barcode_metrics.csv")),
                "gex_molecule_info": os.path.abspath(os.path.join(sample, "outs", "gex_molecule_info.h5"))
            }
            aggr_rows.append(row)

        aggr_csv_df = pd.DataFrame(aggr_rows)
        aggr_csv_df.to_csv(output.aggr_csv, index=False)


rule cellranger_arc_aggr:
    """
    Run [cellranger-arc agg](https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/latest/analysis/running-pipelines/command-line-arguments) to aggregate multiple single-cell multiome libraries.

    This rule only runs if there is more than one sample in a batch. If only one sample is present,
    the aggregation step is skipped.

    Inputs:
        - aggr_csv: A CSV file created by `cellranger_arc_aggr_csv` containing metadata for each library
          including paths to ATAC fragments, per-barcode metrics, and GEX molecule info.

    Outputs:
        - A `.done` flag file indicating that the aggregation step has been completed for this batch.

    Parameters:
        - {batch}: The batch identifier used for naming output files and selecting samples to aggregate.
        - {reference_genome}: Path to the reference genome required by `cellranger-arc`.
        - {jobmode}, {normalize}: Optional additional parameters passed to the `cellranger-arc aggr` command.
        - dirs_dict["LOGS_DIR"]: Directory where log and `.done` files are written.
    """
    log: os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.log")
    input:
        aggr_csv = rules.cellranger_arc_aggr_csv.output.aggr_csv
    output:
        touch(os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.done"))
    run:
        if len(batch_to_samples[int(wildcards.batch)]) > 1:
            shell("""
                cellranger-arc aggr --id={wildcards.batch} \
                                    --reference={reference_genome} \
                                    --csv={input.aggr_csv} \
                                    {jobmode} \
                                    {normalize} \
                                    >> {log} 2>&1
            """)
        else:
            print(f"Batch {wildcards.batch} has only one sample. Skipping cellranger-arc aggr step.")
            pass
