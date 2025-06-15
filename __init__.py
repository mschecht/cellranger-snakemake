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
from config_templates import ARC_CONFIG, ATAC_CONFIG, GEX

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

def sanity_check(args):
    """Check for invalid argument combinations."""
    if args.get_default_config not in (False, True):
        filename = args.get_default_config
        if not isinstance(filename, str) or not filename.endswith('.yaml'):
            logger.error("`--get-default-config` must be a .yaml filename")
            sys.exit(1)


def write_default_config(workflow, filename=None):
    """Print or save the default configuration YAML based on the workflow and filename."""
    
    if workflow == "ARC":
        config_data = ARC_CONFIG
        default_filename = 'ARC_default_config.yaml'
    elif workflow == "ATAC":
        config_data = ATAC_CONFIG
        default_filename = 'ATAC_default_config.yaml'
    elif workflow == "GEX":    
        config_data = GEX
        default_filename = 'GEX_default_config.yaml'
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