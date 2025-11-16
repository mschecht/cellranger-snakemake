import os
import sys

import pandas as pd

from cellranger_snakemake.utils.custom_logger import custom_logger

# Function to apply suffix if specified
def get_directories(config):
    """Get directory paths with optional suffix applied"""
    dirs = {}
    base_dirs = config.get("directories")
    suffix = config.get("directories_suffix", "none")
    
    for key, value in base_dirs.items():
        if suffix != "none" and suffix != "":
            # Apply suffix to directory path
            dirs[key] = f"{value}_{suffix}"
        else:
            # Use directory path as-is
            dirs[key] = value
    
    return dirs

def has_underscore(s: str) -> bool:
    """
    Return True if the string contains at least one underscore, else False.
    """
    return "_" in s


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
    expected_columns = {"batch", "capture", "sample", "fastqs"}
    actual_columns = set(df.columns)
    if expected_columns != actual_columns:
        custom_logger.error(f"Expected columns {expected_columns}, found {actual_columns}")
        sys.exit(1)

    valid = True
    for idx, row in df.iterrows():
        if has_underscore(str(row["batch"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Batch name '{row['batch']}' cannot contain underscores.")
            valid = False

        if has_underscore(str(row["capture"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Capture name '{row['capture']}' cannot contain underscores.")
            valid = False
            
        fastq_path = row["fastqs"]
        error_location = f"Row {idx + 2} in '{filepath}'"
        if not isinstance(fastq_path, str):
            custom_logger.error(f"{error_location}: Path is not a string: {fastq_path}")
            valid = False
            continue

        # Try to identify if fastqs column has more than one path
        DELIMITERS = [",", ";", ":", " "]
        # First, check for disallowed delimiters
        for delim in DELIMITERS[1:]:  # skip comma (allowed)
            if delim in fastq_path:
                custom_logger.error(f"{error_location}: Path contains invalid delimiter '{delim}'. Only comma-separated paths are allowed.")
                sys.exit(1)

        # Split the paths by comma
        suspected_paths = [path.strip() for path in fastq_path.split(",") if path.strip()]
        # Convert relative paths to absolute and validate existence
        absolute_paths = []
        for path in suspected_paths:
            if not os.path.isabs(path):
                # Convert relative path to absolute
                abs_path = os.path.abspath(path)
                custom_logger.warning(f"{error_location}: Converting relative path '{path}' to absolute path '{abs_path}'")
                path = abs_path
            
            if not os.path.exists(path):
                custom_logger.error(f"{error_location}: Path does not exist: {path}")
                valid = False
            else:
                absolute_paths.append(path)
        
        # Update the dataframe with absolute paths
        if absolute_paths:
            df.at[idx, "fastqs"] = ",".join(absolute_paths)

    if valid:
        custom_logger.info(f"{filepath} file format is valid.")
        return df
    else:
        custom_logger.error(f"Some errors were found in {filepath}.")
        sys.exit(1)