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
#  - batch: Batch number which groups captures together
#  - capture: The ID of the capture
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

df = utils.sanity_check_libraries_list_tsv(
    libraries_file, 
    expected_columns={"batch", "capture", "CSV"}, 
    path_column="CSV", 
    file_extension=".csv"
)

# Collect a summary and validate library CSV files
summary_dict = {}
batch_to_captures = defaultdict(list)

for idx, row in df.iterrows():
    batch = row["batch"]
    capture = row["capture"] 
    csv_path = row["CSV"]
    
    # Read and validate the library CSV
    try:
        lib_df = pd.read_csv(csv_path)
        required_lib_columns = {"fastqs", "sample", "library_type"}
        if not required_lib_columns.issubset(lib_df.columns):
            raise ValueError(f"Library CSV missing required columns: {required_lib_columns - set(lib_df.columns)}")
        
        # Validate and convert fastqs paths to absolute paths
        for lib_idx, lib_row in lib_df.iterrows():
            fastq_path = lib_row["fastqs"]
            library_type = lib_row["library_type"]
            
            if not isinstance(fastq_path, str):
                raise ValueError(f"Library CSV row {lib_idx + 1}: fastqs path is not a string: {fastq_path}")
            
            # Convert relative path to absolute
            if not os.path.isabs(fastq_path):
                abs_fastq_path = os.path.abspath(fastq_path)
                custom_logger.warning(f"CSV '{csv_path}', row {lib_idx + 1} ({library_type}): Converting relative path '{fastq_path}' to absolute path '{abs_fastq_path}'")
                lib_df.at[lib_idx, "fastqs"] = abs_fastq_path
                fastq_path = abs_fastq_path
            
            # Validate path exists
            if not os.path.exists(fastq_path):
                raise ValueError(f"Library CSV row {lib_idx + 1} ({library_type}): fastqs path does not exist: {fastq_path}")
        
        # Save the updated CSV with absolute paths
        lib_df.to_csv(csv_path, index=False)
        custom_logger.info(f"Updated CSV file '{csv_path}' with absolute paths")
        
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
        custom_logger.error(f"Error reading library CSV for capture {capture}: {e}")
        sys.exit(1)
    
    out_dir = os.path.join(dirs_dict["CELLRANGERARC_COUNT_DIR"], capture, "outs")

    summary_dict[capture] = {
        "batch": batch,
        "Library_path": csv_path,
        "ATAC_path": atac_path,
        "GEX_path": gex_path,
        "Output_dir": out_dir
    }

    batch_to_captures[batch].append(capture)

# Print summary to stderr to avoid interfering with DAG output

custom_logger.info(f"Found {len(summary_dict)} captures across {len(batch_to_captures)} batches:")
for batch, captures in batch_to_captures.items():
    custom_logger.info(f"Batch {batch}: {len(captures)} captures")

rule all:
    input:
        expand(
            os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.done"),
            batch=batch_to_captures.keys()
        )

rule cellranger_arc_count:
    """
    [cellranger-arc count](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/command-line-arguments#count)

    This program counts ATAC and gene expression reads from a single library against a reference genome.
    """
    input:
        reference = reference_genome,
        library_csv = lambda wc: summary_dict[wc.capture]["Library_path"]
    output:
        done_flag = touch(os.path.join(dirs_dict["LOGS_DIR"], "{capture}_cellranger_arc_count.done")),
    params:
        count_dir = dirs_dict["CELLRANGERARC_COUNT_DIR"],
    log:
        os.path.join(dirs_dict["LOGS_DIR"], "{capture}_cellranger_arc_count.log")
    threads: 8
    resources:
        mem_gb = 64
    run:
        Library_path = summary_dict[wildcards.capture]["Library_path"]
        os.makedirs(params.count_dir, exist_ok=True)

        shell(f"""
            cellranger-arc count --id={wildcards.capture} \
                                 --reference={reference_genome} \
                                 --libraries={Library_path} \
                                 {jobmode} \
                                 {mempercore} \
                                 >> {log} 2>&1
            # mkdir -p {params.count_dir}
            # mv {wildcards.capture} {params.count_dir}
            """)

rule cellranger_arc_aggr_csv:
    """
    Generate a CSV file for use with `cellranger-arc aggr`, combining outputs from multiple captures.

    This rule collects the following files from each `cellranger-arc count` capture directory:
        - `outs/atac_fragments.tsv.gz`
        - `outs/per_barcode_metrics.csv`
        - `outs/gex_molecule_info.h5`

    It assembles a CSV with one row per capture, suitable as input to:
        `cellranger-arc aggr --csv`.

    Output CSV columns:
        - library_id
        - atac_fragments
        - per_barcode_metrics
        - gex_molecule_info
    """
    log: os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr_csv.log")
    input:
        # Require counts for all captures in batch before aggregation
        done_flags = lambda wildcards: [
            os.path.join(dirs_dict["LOGS_DIR"], f"{capture}_cellranger_arc_count.done")
            for capture in batch_to_captures[int(wildcards.batch)]
        ]
    output:
        aggr_csv = os.path.join(dirs_dict["CELLRANGERARC_AGGR_DIR"], "{batch}", "{batch}_aggr.csv")
    run:
        batch = int(wildcards.batch)
        aggr_rows = []

        for capture in batch_to_captures[batch]:
            row = {
                "library_id": capture,
                "atac_fragments": os.path.abspath(os.path.join(capture, "outs", "atac_fragments.tsv.gz")),
                "per_barcode_metrics": os.path.abspath(os.path.join(capture, "outs", "per_barcode_metrics.csv")),
                "gex_molecule_info": os.path.abspath(os.path.join(capture, "outs", "gex_molecule_info.h5"))
            }
            aggr_rows.append(row)

        # Create output directory
        os.makedirs(os.path.dirname(output.aggr_csv), exist_ok=True)
        
        aggr_csv_df = pd.DataFrame(aggr_rows)
        aggr_csv_df.to_csv(output.aggr_csv, index=False)
        
        print(f"Created aggregation CSV for batch {batch} with {len(aggr_rows)} captures", file=sys.stderr)

rule cellranger_arc_aggr:
    """
    Run [cellranger-arc aggr](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/command-line-arguments) 
    to aggregate multiple single-cell multiome libraries.

    This rule only runs if there is more than one capture in a batch. If only one capture is present,
    the aggregation step is skipped.
    """
    log: os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.log")
    input:
        aggr_csv = rules.cellranger_arc_aggr_csv.output.aggr_csv
    output:
        done_flag = touch(os.path.join(dirs_dict["LOGS_DIR"], "{batch}_cellranger_arc_aggr.done")),
        aggr_output = directory(os.path.join(dirs_dict["CELLRANGERARC_AGGR_DIR"], "{batch}", "outs"))
    threads: 8
    resources:
        mem_gb = 64
    run:
        batch_captures = batch_to_captures[int(wildcards.batch)]
        
        if len(batch_captures) > 1:
            shell(f"""
                cellranger-arc aggr --id={wildcards.batch} \
                                    --reference={reference_genome} \
                                    --csv={{input.aggr_csv}} \
                                    {jobmode} \
                                    {normalize} \
                                    >> {{log}} 2>&1
            """)
            print(f"Aggregated {len(batch_captures)} captures for batch {wildcards.batch}", file=sys.stderr)
        else:
            print(f"Batch {wildcards.batch} has only one capture ({batch_captures[0]}). Skipping cellranger-arc aggr step.", file=sys.stderr)
            # Create empty output directory to satisfy the rule
            os.makedirs(output.aggr_output, exist_ok=True)
            with open(os.path.join(output.aggr_output, "single_capture_batch.txt"), "w") as f:
                f.write(f"This batch contained only one capture: {batch_captures[0]}\n")
                f.write("Aggregation was skipped.\n")