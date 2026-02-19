"""Cell type annotation workflow rules for automated cell type identification."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger
from cellranger_snakemake.config_validator import parse_output_directories

# Get centralized output directories
OUTPUT_DIRS = parse_output_directories(config)
ANNDATA_DIR = OUTPUT_DIRS["anndata_dir"]
ANNOTATION_DIR = OUTPUT_DIRS["celltype_annotation_dir"]
LOGS_DIR = OUTPUT_DIRS["logs_dir"]


# ============================================================================
# CELLTYPIST
# ============================================================================

if config.get("celltype_annotation") and config["celltype_annotation"]["method"] == "celltypist":
    CELLTYPIST_PARAMS = config["celltype_annotation"]["celltypist"]

    custom_logger.info("Cell Type Annotation: Using celltypist method")

    rule run_celltypist:
        """Run Celltypist cell type annotation on per-capture AnnData object."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.done")
        output:
            tsv = os.path.join(ANNOTATION_DIR, "{batch}_{capture}_celltypist.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_celltypist.done"))
        params:
            model = CELLTYPIST_PARAMS.get("model"),
            majority_voting = CELLTYPIST_PARAMS.get("majority_voting", False),
            over_clustering = CELLTYPIST_PARAMS.get("over_clustering", None),
            min_prop = CELLTYPIST_PARAMS.get("min_prop", 0.5)
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_celltypist.log")
        script:
            "../scripts/run_celltypist.py"


# ============================================================================
# SCANVI (scvi-tools)
# ============================================================================

if config.get("celltype_annotation") and config["celltype_annotation"]["method"] == "scanvi":
    SCANVI_PARAMS = config["celltype_annotation"]["scanvi"]

    custom_logger.info("Cell Type Annotation: Using scANVI (scvi-tools) method")

    rule run_scanvi:
        """Run scANVI reference-based annotation on per-capture AnnData object."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            reference = SCANVI_PARAMS.get("reference_path"),
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.done")
        output:
            tsv = os.path.join(ANNOTATION_DIR, "{batch}_{capture}_scanvi.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_scanvi.done"))
        params:
            label_key = SCANVI_PARAMS.get("label_key"),
            n_hidden = SCANVI_PARAMS.get("n_hidden", 128),
            n_latent = SCANVI_PARAMS.get("n_latent", 30),
            n_layers = SCANVI_PARAMS.get("n_layers", 2),
            max_epochs = SCANVI_PARAMS.get("max_epochs", 400)
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_scanvi.log")
        conda:
            "../envs/scvi-tools.yaml"
        script:
            "../scripts/run_scanvi.py"


# ============================================================================
# DECOUPLER MARKERS
# ============================================================================

if config.get("celltype_annotation") and config["celltype_annotation"]["method"] == "decoupler_markers":
    DECOUPLER_PARAMS = config["celltype_annotation"]["decoupler_markers"]

    custom_logger.info("Cell Type Annotation: Using decoupler marker-based method")

    rule run_decoupler_markers:
        """Run decoupler marker-based annotation on per-capture AnnData object."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.done")
        output:
            tsv = os.path.join(ANNOTATION_DIR, "{batch}_{capture}_decoupler.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_decoupler.done"))
        params:
            marker_database = DECOUPLER_PARAMS.get("marker_database"),
            method = DECOUPLER_PARAMS.get("method", "ulm"),
            min_score = DECOUPLER_PARAMS.get("min_score", 0.5)
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_decoupler.log")
        script:
            "../scripts/run_decoupler_markers.py"


# ============================================================================
# CELLTYPIST CUSTOM (train custom model)
# ============================================================================

if config.get("celltype_annotation") and config["celltype_annotation"]["method"] == "celltypist_custom":
    CELLTYPIST_CUSTOM_PARAMS = config["celltype_annotation"]["celltypist_custom"]

    custom_logger.info("Cell Type Annotation: Using celltypist with custom trained model")

    # Note: Custom training would be a separate step that produces a model file
    # Then this rule uses that model just like regular celltypist
    # For now, this is a placeholder - actual implementation would need a training rule first

    rule run_celltypist_custom:
        """Run Celltypist with custom-trained model on per-capture AnnData object."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            model = CELLTYPIST_CUSTOM_PARAMS.get("training_data"),  # Pre-trained model path
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.done")
        output:
            tsv = os.path.join(ANNOTATION_DIR, "{batch}_{capture}_celltypist_custom.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_celltypist_custom.done"))
        params:
            model = CELLTYPIST_CUSTOM_PARAMS.get("training_data"),
            majority_voting = CELLTYPIST_CUSTOM_PARAMS.get("majority_voting", False),
            over_clustering = CELLTYPIST_CUSTOM_PARAMS.get("over_clustering", None),
            min_prop = CELLTYPIST_CUSTOM_PARAMS.get("min_prop", 0.5)
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_celltypist_custom.log")
        script:
            "../scripts/run_celltypist.py"  # Same script, different model
