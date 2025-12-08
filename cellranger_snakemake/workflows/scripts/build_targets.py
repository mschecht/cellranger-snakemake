"""Build target files for rule all based on pipeline configuration."""

import pandas as pd
from pathlib import Path


def build_all_targets(config, enabled_steps):
    """
    Build list of final output files based on enabled steps.
    
    Args:
        config: Snakemake config dictionary
        enabled_steps: List of enabled pipeline step names
        
    Returns:
        list: Paths to final output files
    """
    targets = []
    
    # Cell Ranger outputs
    if "cellranger_gex" in enabled_steps:
        targets.extend(get_cellranger_gex_outputs(config))
    if "cellranger_atac" in enabled_steps:
        targets.extend(get_cellranger_atac_outputs(config))
    if "cellranger_arc" in enabled_steps:
        targets.extend(get_cellranger_arc_outputs(config))
    
    # Demux outputs
    if "demultiplexing" in enabled_steps:
        targets.extend(get_demux_outputs(config))
    
    # Doublet detection outputs
    if "doublet_detection" in enabled_steps:
        targets.extend(get_doublet_outputs(config))
    
    # Annotation outputs
    if "celltype_annotation" in enabled_steps:
        targets.extend(get_annotation_outputs(config))
    
    return targets


def parse_libraries_file(libraries_path):
    """
    Parse libraries TSV file to get sample information.
    
    Args:
        libraries_path: Path to libraries TSV file
        
    Returns:
        list: Sample/capture names
    """
    df = pd.read_csv(libraries_path, sep="\t")
    
    # For GEX/ATAC, samples are in 'capture' or 'sample' column
    if 'capture' in df.columns:
        return df['capture'].unique().tolist()
    elif 'sample' in df.columns:
        return df['sample'].unique().tolist()
    else:
        raise ValueError(f"Libraries file must have 'capture' or 'sample' column: {libraries_path}")


def get_cellranger_gex_outputs(config):
    """
    Get GEX output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to GEX done files
    """
    import os
    gex_config = config["cellranger_gex"]
    output_dir = config.get("output_dir", "./output")
    dirs = gex_config.get("directories", {})
    logs_dir = os.path.join(output_dir, dirs.get("LOGS_DIR", "00_LOGS"))
    
    # Parse libraries to get batches
    df = pd.read_csv(gex_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(f"{logs_dir}/{batch}_gex_aggr.done")
    
    return outputs


def get_cellranger_atac_outputs(config):
    """
    Get ATAC output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to ATAC done files
    """
    import os
    atac_config = config["cellranger_atac"]
    output_dir = config.get("output_dir", "./output")
    dirs = atac_config.get("directories", {})
    logs_dir = os.path.join(output_dir, dirs.get("LOGS_DIR", "00_LOGS"))
    
    # Parse libraries to get batches
    df = pd.read_csv(atac_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(f"{logs_dir}/{batch}_atac_aggr.done")
    
    return outputs

def get_cellranger_arc_outputs(config):
    """
    Get ARC output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to ARC done files
    """
    import os
    arc_config = config["cellranger_arc"]
    output_dir = config.get("output_dir", "./output")
    dirs = arc_config.get("directories", {})
    logs_dir = os.path.join(output_dir, dirs.get("LOGS_DIR", "00_LOGS"))
    
    # Parse libraries to get batches
    df = pd.read_csv(arc_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(f"{logs_dir}/{batch}_arc_aggr.done")
    
    return outputs


def get_demux_outputs(config):
    """
    Get demultiplexing output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to demux output files
    """
    demux_config = config["demultiplexing"]
    method = demux_config["method"]
    output_dir = config.get("output_dir", "./output")
    
    # Get sample IDs
    sample_ids = list(config.get("samples", {}).keys())
    
    outputs = []
    for sample in sample_ids:
        outputs.append(f"{output_dir}/demux/{method}/{sample}_demux_results.csv")
    
    return outputs


def get_doublet_outputs(config):
    """
    Get doublet detection output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to doublet detection output files
    """
    doublet_config = config["doublet_detection"]
    method = doublet_config["method"]
    output_dir = config.get("output_dir", "./output")
    
    # Get sample IDs
    sample_ids = list(config.get("samples", {}).keys())
    
    outputs = []
    for sample in sample_ids:
        outputs.append(f"{output_dir}/doublets/{method}/{sample}_doublet_scores.csv")
        outputs.append(f"{output_dir}/doublets/{method}/{sample}_doublet_predictions.csv")
    
    return outputs


def get_annotation_outputs(config):
    """
    Get cell type annotation output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to annotation output files
    """
    annot_config = config["celltype_annotation"]
    method = annot_config["method"]
    output_dir = config.get("output_dir", "./output")
    
    # Get sample IDs
    sample_ids = list(config.get("samples", {}).keys())
    
    outputs = []
    for sample in sample_ids:
        outputs.append(f"{output_dir}/annotation/{method}/{sample}_cell_types.csv")
        outputs.append(f"{output_dir}/annotation/{method}/{sample}_annotated.h5ad")
    
    return outputs
