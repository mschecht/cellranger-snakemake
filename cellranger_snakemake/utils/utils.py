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


def sanity_check_libraries_list_tsv(filepath, log_file=None, expected_columns=None, path_column="fastqs", file_extension=None):
    """
    Check the format of the libraries list TSV file.
    Args:
        filepath (str): Path to the TSV file.
        log_file (str, optional): Path to the log file. Defaults to None.
        expected_columns (set, optional): Expected column names. Defaults to {"batch", "capture", "sample", "fastqs"}.
        path_column (str, optional): Name of the column containing file paths. Defaults to "fastqs".
        file_extension (str, optional): Required file extension (e.g., ".csv"). If None, checks for directory existence.
    Returns:
        pd.DataFrame: Valid dataframe or exits on error
    """
    # Set default expected columns if not provided
    if expected_columns is None:
        expected_columns = {"batch", "capture", "sample", "fastqs"}

    # Read the file
    try:
        df = pd.read_csv(filepath, sep="\t")
        df.columns = df.columns.str.strip()

    except Exception as e:
        custom_logger.error(f"Pandas could not read your libraries filepath here: '{filepath}'. "
                            f"This was the error: {e}")
        sys.exit(1)

    # Validate headers
    actual_columns = set(df.columns)
    if expected_columns != actual_columns:
        custom_logger.error(f"Expected columns {expected_columns}, found {actual_columns}")
        sys.exit(1)

    valid = True
    for idx, row in df.iterrows():
        if has_underscore(str(row["batch"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Batch name '{row['batch']}' cannot contain underscores.")
            valid = False

        if "capture" in expected_columns and has_underscore(str(row["capture"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Capture name '{row['capture']}' cannot contain underscores.")
            valid = False
            
        file_path = row[path_column]
        error_location = f"Row {idx + 2} in '{filepath}'"
        if not isinstance(file_path, str):
            custom_logger.error(f"{error_location}: Path is not a string: {file_path}")
            valid = False
            continue

        # Handle different path formats (fastqs vs CSV)
        if path_column == "fastqs":
            # For GEX/ATAC: Handle comma-separated fastq directories
            DELIMITERS = [",", ";", ":", " "]
            # First, check for disallowed delimiters
            for delim in DELIMITERS[1:]:  # skip comma (allowed)
                if delim in file_path:
                    custom_logger.error(f"{error_location}: Path contains invalid delimiter '{delim}'. Only comma-separated paths are allowed.")
                    sys.exit(1)

            # Split the paths by comma
            suspected_paths = [path.strip() for path in file_path.split(",") if path.strip()]
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
                df.at[idx, path_column] = ",".join(absolute_paths)
        else:
            # For ARC CSV or other single file formats
            # Check file extension if specified
            if file_extension and not file_path.endswith(file_extension):
                custom_logger.error(f"{error_location}: The path '{file_path}' does not end with '{file_extension}'. Please ensure all paths have the proper file extension.")
                valid = False
                continue
            
            # Convert relative path to absolute
            if not os.path.isabs(file_path):
                abs_path = os.path.abspath(file_path)
                custom_logger.warning(f"{error_location}: Converting relative path '{file_path}' to absolute path '{abs_path}'")
                file_path = abs_path
                df.at[idx, path_column] = abs_path
            
            if not os.path.exists(file_path):
                custom_logger.error(f"{error_location}: File does not exist: {file_path}")
                valid = False

    if valid:
        custom_logger.info(f"{filepath} file format is valid.")
        return df
    else:
        custom_logger.error(f"Some errors were found in {filepath}.")
        sys.exit(1)