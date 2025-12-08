"""Configuration parsing utilities for Snakemake workflows."""

import sys
from pathlib import Path


def get_enabled_steps(config):
    """
    Extract enabled pipeline steps from config.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Names of enabled pipeline steps
    """
    steps = []
    
    # Check each step
    for step in ["cellranger_gex", "cellranger_atac", "cellranger_arc", 
                 "demultiplexing", "doublet_detection", "celltype_annotation"]:
        if config.get(step) and config[step].get("enabled", True):
            steps.append(step)
    
    return steps


def get_step_method(config, step):
    """
    Get the method selected for a pipeline step.
    
    Args:
        config: Snakemake config dictionary
        step: Name of the pipeline step
        
    Returns:
        str: Method name or None
    """
    step_config = config.get(step, {})
    return step_config.get("method")


def get_step_config(config, step):
    """
    Get configuration for a specific pipeline step.
    
    Args:
        config: Snakemake config dictionary
        step: Name of the pipeline step
        
    Returns:
        dict: Step configuration
    """
    return config.get(step, {})


def get_samples(config):
    """
    Get sample information from config.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        dict: Sample metadata dictionary
    """
    return config.get("samples", {})


def get_sample_ids(config):
    """
    Get list of sample IDs.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Sample IDs
    """
    return list(config.get("samples", {}).keys())


def get_output_dir(config):
    """
    Get base output directory.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        str: Output directory path
    """
    return config.get("output_dir", "./output")


def get_resources(config):
    """
    Get resource configuration.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        dict: Resource configuration
    """
    return config.get("resources", {})


def get_hpc_config(config):
    """
    Get HPC configuration.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        dict: HPC configuration
    """
    return config.get("hpc", {})
