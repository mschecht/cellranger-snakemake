"""Doublet detection workflow rules for identifying doublets in scRNA-seq data."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger
from cellranger_snakemake.config_validator import parse_output_directories

# Get centralized output directories
OUTPUT_DIRS = parse_output_directories(config)
ANNDATA_DIR = OUTPUT_DIRS["anndata_dir"]
DOUBLET_DIR = OUTPUT_DIRS["doublet_detection_dir"]
LOGS_DIR = OUTPUT_DIRS["logs_dir"]


# ============================================================================
# SCRUBLET
# ============================================================================

if config.get("doublet_detection") and config["doublet_detection"]["method"] == "scrublet":
    SCRUBLET_PARAMS = config["doublet_detection"]["scrublet"]

    custom_logger.info("Doublet Detection: Using scrublet method")

    rule run_scrublet:
        """Run Scrublet doublet detection on per-capture AnnData object."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.done")
        output:
            tsv = os.path.join(DOUBLET_DIR, "{batch}_{capture}_scrublet.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_scrublet.done"))
        params:
            expected_doublet_rate = SCRUBLET_PARAMS.get("expected_doublet_rate", 0.06),
            min_counts = SCRUBLET_PARAMS.get("min_counts", 2),
            min_cells = SCRUBLET_PARAMS.get("min_cells", 3),
            min_gene_variability_pctl = SCRUBLET_PARAMS.get("min_gene_variability_pctl", 85.0),
            n_prin_comps = SCRUBLET_PARAMS.get("n_prin_comps", 30)
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_scrublet.log")
        script:
            "../scripts/run_scrublet.py"


# ============================================================================
# SOLO (scvi-tools)
# ============================================================================

if config.get("doublet_detection") and config["doublet_detection"]["method"] == "solo":
    SOLO_PARAMS = config["doublet_detection"]["solo"]

    custom_logger.info("Doublet Detection: Using SOLO (scvi-tools) method")

    rule run_solo:
        """Run SOLO doublet detection on per-capture AnnData object."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.done")
        output:
            tsv = os.path.join(DOUBLET_DIR, "{batch}_{capture}_solo.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_solo.done"))
        params:
            n_hidden = SOLO_PARAMS.get("n_hidden", 128),
            n_latent = SOLO_PARAMS.get("n_latent", 64),
            n_layers = SOLO_PARAMS.get("n_layers", 1),
            learning_rate = SOLO_PARAMS.get("learning_rate", 1e-3),
            max_epochs = SOLO_PARAMS.get("max_epochs", 400)
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_solo.log")
        conda:
            "../envs/scvi-tools.yaml"
        script:
            "../scripts/run_solo.py"
