#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Snakemake Cell Ranger ARC Workflow Runner
A command-line tool to run Cell Ranger ARC workflows with Snakemake

[Cell Ranger software documention](https://www.10xgenomics.com/software#data-processing-pipelines)
"""

import os
import sys
import yaml
import logging
import argparse
import subprocess

from pathlib import Path
from config_templates import ARC_CONFIG, ATAC_CONFIG, RNAseq_CONFIG

__version__ = "1.0.0"
__description__ = "Snakemake wrapper for Cell Ranger ARC workflows"

class ColorFormatter(logging.Formatter):
    # ANSI escape codes for colors
    COLORS = {
        'DEBUG': "\033[37m",    # White
        'INFO': "\033[36m",     # Cyan
        'WARNING': "\033[33m",  # Yellow
        'ERROR': "\033[31m",    # Red
        'CRITICAL': "\033[41m", # Red background
    }
    RESET = "\033[0m"

    def format(self, record):
        color = self.COLORS.get(record.levelname, self.RESET)
        message = super().format(record)
        return f"{color}{message}{self.RESET}"

# Initialize logger
logger = logging.getLogger("snakemake_runner")
logger.setLevel(logging.DEBUG)  # Adjust level as needed

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

formatter = ColorFormatter("[%(levelname)s] %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)


# Default configuration template
# ARC_CONFIG = {
#     "reference": "/path/to/reference-genome",
#     "libraries": "/path/to/libraries_list",
#     "HPC_mode": "",
#     "mempercore": 8,
#     "normalize": "none",
# }

# ATAC_CONFIG = {
#     "samples": {
#         "sample1": {
#             "fastqs": "/path/to/fastqs/sample1",
#             "libraries_csv": "/path/to/sample1_libraries.csv"
#         }
#     },
#     "reference": "/path/to/cellranger-arc-reference",
#     "output_dir": "results",
#     "cellranger_arc": {
#         "extra_args": "--force-cells=5000",
#         "localcores": 8,
#         "localmem": 64
#     },
#     "resources": {
#         "default_threads": 8,
#         "default_mem_gb": 64
#     }
# }

# RNA_CONFIG = {
#     "samples": {
#         "sample1": {
#             "fastqs": "/path/to/fastqs/sample1",
#             "libraries_csv": "/path/to/sample1_libraries.csv"
#         }
#     },
#     "reference": "/path/to/cellranger-arc-reference",
#     "output_dir": "results",
#     "cellranger_arc": {
#         "extra_args": "--force-cells=5000",
#         "localcores": 8,
#         "localmem": 64
#     },
#     "resources": {
#         "default_threads": 8,
#         "default_mem_gb": 64
#     }
# }

# Default Snakefile template
DEFAULT_SNAKEFILE = '''
import os
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Get all samples
SAMPLES = list(config["samples"].keys())

rule all:
    input:
        expand("{output_dir}/{sample}/outs/web_summary.html", 
               output_dir=config["output_dir"], 
               sample=SAMPLES)

rule cellranger_arc_count:
    input:
        fastqs = lambda wildcards: config["samples"][wildcards.sample]["fastqs"],
        libraries = lambda wildcards: config["samples"][wildcards.sample]["libraries_csv"],
        reference = config["reference"]
    output:
        directory("{output_dir}/{sample}"),
        web_summary = "{output_dir}/{sample}/outs/web_summary.html",
        filtered_feature_bc_matrix = directory("{output_dir}/{sample}/outs/filtered_feature_bc_matrix"),
        raw_feature_bc_matrix = directory("{output_dir}/{sample}/outs/raw_feature_bc_matrix")
    params:
        sample_id = "{sample}",
        extra = config["cellranger_arc"].get("extra_args", ""),
        output_prefix = lambda wildcards, output: os.path.dirname(output[0])
    threads: 
        config["cellranger_arc"].get("localcores", 8)
    resources:
        mem_gb = config["cellranger_arc"].get("localmem", 64)
    shell:
        """
        cd {params.output_prefix}
        cellranger-arc count \\
            --id={params.sample_id} \\
            --reference={input.reference} \\
            --libraries={input.libraries} \\
            --fastqs={input.fastqs} \\
            --localcores={threads} \\
            --localmem={resources.mem_gb} \\
            {params.extra}
        """
'''

def sanity_check(args):
    """Check for invalid argument combinations."""
    if args.get_default_config not in (False, True):
        filename = args.get_default_config
        if not isinstance(filename, str) or not filename.endswith('.yaml'):
            logger.error("`--get-default-config` must be a .yaml filename")
            sys.exit(1)


def write_default_config(workflow, filename=None):
    """Print or save the default configuration YAML based on the workflow and filename."""
    
    if workflow == "arc":
        config_data = ARC_CONFIG
        default_filename = 'ARC_default_config.yaml'
    elif workflow == "atac":
        config_data = ATAC_CONFIG
        default_filename = 'ATAC_default_config.yaml'
    elif workflow == "rna":    
        config_data = RNAseq_CONFIG
        default_filename = 'RNA_default_config.yaml'
    else:
        logger.error(f"Unknown workflow type: {workflow}")
        sys.exit(1)

    if filename is None:
        # Print YAML to stdout
        print(yaml.dump(config_data, indent=2))
        logger.info(f"Writing '{workflow}' default config YAML to '{default_filename}'")
    else:
        # Write YAML to filename
        logger.info(f"Writing '{workflow}' default config YAML to '{filename}'")
        with open(filename, 'w') as f:
            yaml.dump(config_data, f, indent=2)


def create_default_files(output_dir="."):
    """Create default Snakefile and config.yaml"""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Write Snakefile
    snakefile_path = output_path / "Snakefile"
    with open(snakefile_path, 'w') as f:
        f.write(DEFAULT_SNAKEFILE.strip())
    
    # Write config.yaml
    config_path = output_path / "config.yaml"
    with open(config_path, 'w') as f:
        f.write("# Configuration for Cell Ranger ARC Snakemake workflow\n")
        f.write("# Modify the paths and parameters below\n\n")
        yaml.dump(DEFAULT_CONFIG, f, default_flow_style=False, indent=2)
    
    print(f"Created default files in {output_path}:")
    print(f"  - {snakefile_path}")
    print(f"  - {config_path}")
    print("\nEdit config.yaml with your specific paths and parameters, then run:")
    print(f"  snakemake-run-cellranger.py --run")

def run_snakemake(config_file="config.yaml", snakefile="Snakefile", dry_run=False, cores=None, additional_args=""):
    """Run the Snakemake workflow"""
    
    # Check if files exist
    if not os.path.exists(config_file):
        print(f"Error: Configuration file '{config_file}' not found.")
        print("Run with --get-default-config to create a template.")
        sys.exit(1)
    
    if not os.path.exists(snakefile):
        print(f"Error: Snakefile '{snakefile}' not found.")
        print("Run with --create-default-files to create template files.")
        sys.exit(1)
    
    # Build snakemake command
    cmd = ["snakemake"]
    
    if snakefile != "Snakefile":
        cmd.extend(["-s", snakefile])
    
    if config_file != "config.yaml":
        cmd.extend(["--configfile", config_file])
    
    if dry_run:
        cmd.append("-n")
    
    if cores:
        cmd.extend(["-j", str(cores)])
    else:
        cmd.extend(["-j", "1"])  # Default to 1 core
    
    # Add additional arguments
    if additional_args:
        cmd.extend(additional_args.split())
    
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        if not dry_run:
            print("Workflow completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error running Snakemake: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: Snakemake not found. Please install Snakemake first.")
        print("  conda install -c bioconda snakemake")
        sys.exit(1)