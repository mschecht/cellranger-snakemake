"""Build target files for rule all based on pipeline configuration."""

import os
import pandas as pd

from pathlib import Path
from cellranger_snakemake.config_validator import parse_output_directories

def build_all_targets(config, enabled_steps):
    """
    Build list of final output files based on enabled steps.
    
    Args:
        config: Snakemake config dictionary
        enabled_steps: List of enabled pipeline step names
        
    Returns:
        list: Paths to final output files
    """
    # Defensive check - ensure enabled_steps is a list
    if not enabled_steps:
        return []
    
    if not isinstance(enabled_steps, list):
        enabled_steps = list(enabled_steps)
    
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
    if not config.get("cellranger_gex"):
        return []
    
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    gex_config = config["cellranger_gex"]
    
    # Parse libraries to get batches
    df = pd.read_csv(gex_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(os.path.join(logs_dir, f"{batch}_gex_aggr.done"))
    
    return outputs


def get_cellranger_atac_outputs(config):
    """
    Get ATAC output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to ATAC done files
    """
    if not config.get("cellranger_atac"):
        return []
    
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    atac_config = config["cellranger_atac"]
    
    # Parse libraries to get batches
    df = pd.read_csv(atac_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(os.path.join(logs_dir, f"{batch}_atac_aggr.done"))
    
    return outputs

def get_cellranger_arc_outputs(config):
    """
    Get ARC output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to ARC done files
    """
    if not config.get("cellranger_arc"):
        return []
    
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    arc_config = config["cellranger_arc"]
    
    # Parse libraries to get batches
    df = pd.read_csv(arc_config["libraries"], sep="\t")
    batches = df['batch'].unique().tolist()
    
    # Return done files for each batch
    outputs = []
    for batch in batches:
        outputs.append(os.path.join(logs_dir, f"{batch}_arc_aggr.done"))
    
    return outputs


def get_demux_outputs(config):
    """
    Get demultiplexing output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to demux output files
    """
    if not config.get("demultiplexing"):
        return []
    
    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    demux_config = config["demultiplexing"]
    method = demux_config["method"]
    
    outputs = []
    
    # Try to get batches from cellranger GEX if available
    if config.get("cellranger_gex"):
        gex_config = config["cellranger_gex"]
        libraries_path = gex_config["libraries"]
        df = pd.read_csv(libraries_path, sep="\t")
        batches = df['batch'].unique().tolist()
        captures = df['capture'].unique().tolist()
        
        # For vireo: add cellsnp-lite and vireo outputs per batch-capture
        if method == "vireo":
            for batch in batches:
                for capture in captures:
                        outputs.append(os.path.join(logs_dir, f"cellsnp_output_{batch}_{capture}.done"))
                        outputs.append(os.path.join(logs_dir, f"vireo_output_{batch}_{capture}.done"))
    
    return outputs


def get_doublet_outputs(config):
    """
    Get doublet detection output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to doublet detection output files
    """
    if not config.get("doublet_detection"):
        return []
    
    output_dirs = parse_output_directories(config)
    doublet_config = config["doublet_detection"]
    method = doublet_config["method"]
    doublet_dir = output_dirs["doublet_detection_dir"]
    
    outputs = []
    
    # Get sample IDs from cellranger GEX if available
    if config.get("cellranger_gex"):
        gex_config = config["cellranger_gex"]
        df = pd.read_csv(gex_config["libraries"], sep="\t")
        captures = df['capture'].unique().tolist()
        
        for sample in captures:
            outputs.append(os.path.join(doublet_dir, f"{sample}_doublet_results.csv"))
    
    return outputs


def get_annotation_outputs(config):
    """
    Get cell type annotation output file paths.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Paths to annotation output files
    """
    if not config.get("celltype_annotation"):
        return []
    
    output_dirs = parse_output_directories(config)
    annot_config = config["celltype_annotation"]
    method = annot_config["method"]
    annotation_dir = output_dirs["celltype_annotation_dir"]
    
    outputs = []
    
    # Get sample IDs from cellranger GEX if available
    if config.get("cellranger_gex"):
        gex_config = config["cellranger_gex"]
        df = pd.read_csv(gex_config["libraries"], sep="\t")
        captures = df['capture'].unique().tolist()
        
        for sample in captures:
            outputs.append(os.path.join(annotation_dir, f"{sample}_annotations.csv"))
    
    return outputs
