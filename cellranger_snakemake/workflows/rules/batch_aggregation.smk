"""Batch aggregation rules - merge per-capture objects into batch-level objects."""

import os
import sys
from pathlib import Path
import pandas as pd

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger
from cellranger_snakemake.config_validator import parse_output_directories

# Get output directories
OUTPUT_DIRS = parse_output_directories(config)
ANNDATA_DIR = OUTPUT_DIRS["anndata_dir"]
BATCH_DIR = os.path.join(config.get("output_dir", "output"), "04_BATCH_OBJECTS")
LOGS_DIR = OUTPUT_DIRS["logs_dir"]

# Helper function to get all captures for a batch
def get_captures_for_batch(batch_id, modality):
    """Get all captures for a given batch."""
    if modality == "gex" and config.get("cellranger_gex"):
        libraries_file = config["cellranger_gex"]["libraries"]
    elif modality == "atac" and config.get("cellranger_atac"):
        libraries_file = config["cellranger_atac"]["libraries"]
    elif modality == "arc" and config.get("cellranger_arc"):
        libraries_file = config["cellranger_arc"]["libraries"]
    else:
        return []

    df = pd.read_csv(libraries_file, sep="\t")
    captures = df[df['batch'] == int(batch_id)]['capture'].tolist()
    return captures


# ============================================================================
# GEX Batch Aggregation
# ============================================================================

if config.get("cellranger_gex"):

    # Get unique batches
    libraries_file = config["cellranger_gex"]["libraries"]
    gex_df = pd.read_csv(libraries_file, sep="\t")
    GEX_BATCHES = gex_df['batch'].unique().tolist()

    custom_logger.info(f"Batch aggregation: Found {len(GEX_BATCHES)} GEX batch(es)")

    def get_gex_percapture_objects(wildcards):
        """Get all per-capture .h5ad files for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "gex")
        return [os.path.join(ANNDATA_DIR, f"{wildcards.batch}_{cap}.h5ad") for cap in captures]

    def get_gex_percapture_done(wildcards):
        """Get all per-capture done flags for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "gex")
        return [os.path.join(LOGS_DIR, f"{wildcards.batch}_{cap}_gex_anndata.done") for cap in captures]

    rule aggregate_gex_batch:
        """Aggregate per-capture GEX objects into batch-level object."""
        input:
            h5ads = get_gex_percapture_objects,
            done_flags = get_gex_percapture_done
        output:
            h5ad = os.path.join(BATCH_DIR, "{batch}_gex.h5ad"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_gex_batch_aggregation.done"))
        params:
            batch = "{batch}",
            modality = "gex"
        log:
            os.path.join(LOGS_DIR, "{batch}_gex_batch_aggregation.log")
        script:
            "../scripts/aggregate_batch.py"


# ============================================================================
# ATAC Batch Aggregation
# ============================================================================

if config.get("cellranger_atac"):

    # Get unique batches
    libraries_file = config["cellranger_atac"]["libraries"]
    atac_df = pd.read_csv(libraries_file, sep="\t")
    ATAC_BATCHES = atac_df['batch'].unique().tolist()

    custom_logger.info(f"Batch aggregation: Found {len(ATAC_BATCHES)} ATAC batch(es)")

    def get_atac_percapture_objects(wildcards):
        """Get all per-capture .h5ad files for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "atac")
        return [os.path.join(ANNDATA_DIR, f"{wildcards.batch}_{cap}.h5ad") for cap in captures]

    def get_atac_percapture_done(wildcards):
        """Get all per-capture done flags for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "atac")
        return [os.path.join(LOGS_DIR, f"{wildcards.batch}_{cap}_atac_anndata.done") for cap in captures]

    rule aggregate_atac_batch:
        """Aggregate per-capture ATAC objects into batch-level object."""
        input:
            h5ads = get_atac_percapture_objects,
            done_flags = get_atac_percapture_done
        output:
            h5ad = os.path.join(BATCH_DIR, "{batch}_atac.h5ad"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_atac_batch_aggregation.done"))
        params:
            batch = "{batch}",
            modality = "atac"
        log:
            os.path.join(LOGS_DIR, "{batch}_atac_batch_aggregation.log")
        script:
            "../scripts/aggregate_batch.py"


# ============================================================================
# ARC Batch Aggregation
# ============================================================================

if config.get("cellranger_arc"):

    # Get unique batches
    libraries_file = config["cellranger_arc"]["libraries"]
    arc_df = pd.read_csv(libraries_file, sep="\t")
    ARC_BATCHES = arc_df['batch'].unique().tolist()

    custom_logger.info(f"Batch aggregation: Found {len(ARC_BATCHES)} ARC batch(es)")

    def get_arc_percapture_objects(wildcards):
        """Get all per-capture .h5mu files for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "arc")
        return [os.path.join(ANNDATA_DIR, f"{wildcards.batch}_{cap}.h5mu") for cap in captures]

    def get_arc_percapture_done(wildcards):
        """Get all per-capture done flags for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "arc")
        return [os.path.join(LOGS_DIR, f"{wildcards.batch}_{cap}_arc_anndata.done") for cap in captures]

    rule aggregate_arc_batch:
        """Aggregate per-capture ARC objects into batch-level object."""
        input:
            h5mus = get_arc_percapture_objects,
            done_flags = get_arc_percapture_done
        output:
            h5mu = os.path.join(BATCH_DIR, "{batch}_arc.h5mu"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_arc_batch_aggregation.done"))
        params:
            batch = "{batch}",
            modality = "arc"
        log:
            os.path.join(LOGS_DIR, "{batch}_arc_batch_aggregation.log")
        script:
            "../scripts/aggregate_batch.py"
