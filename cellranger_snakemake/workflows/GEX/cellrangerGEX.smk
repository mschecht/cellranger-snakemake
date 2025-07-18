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

    custom_logger.info(f"{libraries_file} file format is valid.")

    return df

df = sanity_check_libraries_list_GEX_tsv(libraries_file)

# Collect a summary and validate library CSV files
summary_dict = {}
batch_to_samples = defaultdict(list)

for idx, row in df.iterrows():
    batch = row["batch"]
    ID = row["ID"] 
    sample = row["sample"]
    Fastqs = row["Fastqs"]

    summary_dict[batch] = {
        "ID": ID,
        "sample": sample,
        "Fastqs": Fastqs
    }

print(summary_dict)

