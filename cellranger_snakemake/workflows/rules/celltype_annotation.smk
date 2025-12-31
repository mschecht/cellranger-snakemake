"""Cell type annotation workflow rules for automated cell type identification."""

import os
import sys
from pathlib import Path

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger


# Get annotation config
if config.get("celltype_annotation"):
    ANNOT_CONFIG = config["celltype_annotation"]
    ANNOT_METHOD = ANNOT_CONFIG["method"]
    ANNOT_PARAMS = ANNOT_CONFIG.get("parameters", {})
    ANNOT_SAMPLES = config.get("samples", [])
    
    custom_logger.info(f"Cell Type Annotation: Using {ANNOT_METHOD} method")


# ============================================================================
# CELLTYPIST
# ============================================================================

if config.get("celltype_annotation") and ANNOT_METHOD == "celltypist":
    
    rule celltypist:
        """Run Celltypist for cell type annotation."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            predictions = "{sample}/celltypist/predicted_labels.csv",
            probabilities = "{sample}/celltypist/probability_matrix.csv",
            done = touch("{sample}/celltypist/{sample}_celltypist.done")
        params:
            model = ANNOT_PARAMS.get("model", "Immune_All_Low.pkl"),
            majority_voting = ANNOT_PARAMS.get("majority_voting", False),
            over_clustering = ANNOT_PARAMS.get("over_clustering", None),
            outdir = "{sample}/celltypist"
        log:
            "{sample}/celltypist/{sample}_celltypist.log"
        script:
            "../scripts/run_celltypist.py"


# ============================================================================
# AZIMUTH (R - Seurat)
# ============================================================================

if config.get("celltype_annotation") and ANNOT_METHOD == "azimuth":
    
    rule azimuth:
        """Run Azimuth for reference-based annotation."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            seurat = "{sample}/azimuth/seurat_annotated.rds",
            predictions = "{sample}/azimuth/predictions.csv",
            done = touch("{sample}/azimuth/{sample}_azimuth.done")
        params:
            reference = ANNOT_PARAMS.get("reference", "pbmc"),
            outdir = "{sample}/azimuth"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/azimuth/{sample}_azimuth.log"
        script:
            "../scripts/run_azimuth.R"


# ============================================================================
# SINGLER (R - Bioconductor)
# ============================================================================

if config.get("celltype_annotation") and ANNOT_METHOD == "singler":
    
    rule singler:
        """Run SingleR for reference-based annotation."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            results = "{sample}/singler/singler_results.rds",
            predictions = "{sample}/singler/predictions.csv",
            done = touch("{sample}/singler/{sample}_singler.done")
        params:
            reference = ANNOT_PARAMS.get("reference_name"),
            reference_data = ANNOT_PARAMS.get("reference_data", None),
            label_type = ANNOT_PARAMS.get("label_type", "label.main"),
            outdir = "{sample}/singler"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/singler/{sample}_singler.log"
        script:
            "../scripts/run_singler.R"


# ============================================================================
# SCTYPE (R)
# ============================================================================

if config.get("celltype_annotation") and ANNOT_METHOD == "sctype":
    
    rule sctype:
        """Run ScType for marker-based annotation."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            seurat = "{sample}/sctype/seurat_annotated.rds",
            predictions = "{sample}/sctype/predictions.csv",
            done = touch("{sample}/sctype/{sample}_sctype.done")
        params:
            tissue = ANNOT_PARAMS.get("tissue_type"),
            db_url = ANNOT_PARAMS.get("marker_db_url", None),
            custom_markers = ANNOT_PARAMS.get("custom_marker_file", None),
            outdir = "{sample}/sctype"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/sctype/{sample}_sctype.log"
        script:
            "../scripts/run_sctype.R"
