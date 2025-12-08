"""Doublet detection workflow rules for identifying doublets in scRNA-seq data."""

import os
import sys
from pathlib import Path

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger


# Get doublet detection config
if config.get("doublet_detection"):
    DOUBLET_CONFIG = config["doublet_detection"]
    DOUBLET_METHOD = DOUBLET_CONFIG["method"]
    DOUBLET_PARAMS = DOUBLET_CONFIG.get("parameters", {})
    DOUBLET_SAMPLES = config.get("samples", [])
    
    custom_logger.info(f"Doublet Detection: Using {DOUBLET_METHOD} method")


# ============================================================================
# SCRUBLET
# ============================================================================

if config.get("doublet_detection") and DOUBLET_METHOD == "scrublet":
    
    rule scrublet:
        """Run Scrublet for doublet detection."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            scores = "{sample}/scrublet/doublet_scores.txt",
            calls = "{sample}/scrublet/doublet_calls.txt",
            done = touch("{sample}/scrublet/{sample}_scrublet.done")
        params:
            expected_rate = DOUBLET_PARAMS.get("expected_doublet_rate", 0.06),
            sim_doublet_ratio = DOUBLET_PARAMS.get("sim_doublet_ratio", 2.0),
            n_neighbors = DOUBLET_PARAMS.get("n_neighbors", None),
            outdir = "{sample}/scrublet"
        log:
            "{sample}/scrublet/{sample}_scrublet.log"
        script:
            "../scripts/run_scrublet.py"


# ============================================================================
# DOUBLETFINDER (R)
# ============================================================================

if config.get("doublet_detection") and DOUBLET_METHOD == "doubletfinder":
    
    rule doubletfinder:
        """Run DoubletFinder for doublet detection."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            seurat = "{sample}/doubletfinder/seurat_object.rds",
            classifications = "{sample}/doubletfinder/classifications.txt",
            done = touch("{sample}/doubletfinder/{sample}_doubletfinder.done")
        params:
            expected_rate = DOUBLET_PARAMS.get("expected_doublet_rate", 0.06),
            pn = DOUBLET_PARAMS.get("pN", 0.25),
            pk = DOUBLET_PARAMS.get("pK", None),
            outdir = "{sample}/doubletfinder"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/doubletfinder/{sample}_doubletfinder.log"
        script:
            "../scripts/run_doubletfinder.R"


# ============================================================================
# SCDS (R - Bioconductor)
# ============================================================================

if config.get("doublet_detection") and DOUBLET_METHOD == "scds":
    
    rule scds:
        """Run scds for doublet detection."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            scores = "{sample}/scds/doublet_scores.txt",
            calls = "{sample}/scds/doublet_calls.txt",
            done = touch("{sample}/scds/{sample}_scds.done")
        params:
            outdir = "{sample}/scds"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/scds/{sample}_scds.log"
        script:
            "../scripts/run_scds.R"


# ============================================================================
# SCDBLFINDER (R - Bioconductor)
# ============================================================================

if config.get("doublet_detection") and DOUBLET_METHOD == "scdblfinder":
    
    rule scdblfinder:
        """Run scDblFinder for doublet detection."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            sce = "{sample}/scdblfinder/sce_object.rds",
            scores = "{sample}/scdblfinder/doublet_scores.txt",
            calls = "{sample}/scdblfinder/doublet_calls.txt",
            done = touch("{sample}/scdblfinder/{sample}_scdblfinder.done")
        params:
            clusters = DOUBLET_PARAMS.get("clusters", None),
            dbr = DOUBLET_PARAMS.get("dbr", None),
            outdir = "{sample}/scdblfinder"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/scdblfinder/{sample}_scdblfinder.log"
        script:
            "../scripts/run_scdblfinder.R"
