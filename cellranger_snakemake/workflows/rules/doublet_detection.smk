"""Doublet detection workflow rules for identifying doublets in scRNA-seq data."""

import os
import sys
from pathlib import Path
from tempfile import gettempdir

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

    # Determine which modality's anndata done flag and file extension to use
    if config.get("cellranger_gex"):
        ANNDATA_DONE_SUFFIX = "gex_anndata"
        ANNDATA_EXT = "h5ad"
    elif config.get("cellranger_atac"):
        ANNDATA_DONE_SUFFIX = "atac_anndata"
        ANNDATA_EXT = "h5ad"
    elif config.get("cellranger_arc"):
        ANNDATA_DONE_SUFFIX = "arc_mudata"
        ANNDATA_EXT = "h5mu"
    else:
        ANNDATA_DONE_SUFFIX = "gex_anndata"
        ANNDATA_EXT = "h5ad"

    rule run_scrublet:
        """Run Scrublet doublet detection on per-capture AnnData object.
        For ARC (MuData), scrublet is run on the GEX (rna) modality."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}." + ANNDATA_EXT),
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_" + ANNDATA_DONE_SUFFIX + ".done")
        output:
            tsv = os.path.join(DOUBLET_DIR, "{batch}_{capture}_scrublet.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_scrublet.done"))
        params:
            filter_cells_min_genes = SCRUBLET_PARAMS.get("filter_cells_min_genes", 1),
            filter_genes_min_cells = SCRUBLET_PARAMS.get("filter_genes_min_cells", 3),
            expected_doublet_rate = SCRUBLET_PARAMS.get("expected_doublet_rate", 0.06),
            min_gene_variability_pctl = SCRUBLET_PARAMS.get("min_gene_variability_pctl", 85.0),
            n_prin_comps = SCRUBLET_PARAMS.get("n_prin_comps", 30),
            sim_doublet_ratio = SCRUBLET_PARAMS.get("sim_doublet_ratio", 2.0),
            threshold = SCRUBLET_PARAMS.get("threshold", None),
            n_neighbors = SCRUBLET_PARAMS.get("n_neighbors", None),
            random_state = SCRUBLET_PARAMS.get("random_state", 0)
        threads: config["doublet_detection"].get("threads", 10)
        resources:
            mem_mb = config["doublet_detection"].get("mem_gb", 64) * 1024, # NOTE: SLURM resources uses M  B
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
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

    # Determine which modality's anndata done flag and file extension to use
    if config.get("cellranger_gex"):
        ANNDATA_DONE_SUFFIX = "gex_anndata"
        ANNDATA_EXT = "h5ad"
    elif config.get("cellranger_atac"):
        ANNDATA_DONE_SUFFIX = "atac_anndata"
        ANNDATA_EXT = "h5ad"
    elif config.get("cellranger_arc"):
        ANNDATA_DONE_SUFFIX = "arc_mudata"
        ANNDATA_EXT = "h5mu"
    else:
        ANNDATA_DONE_SUFFIX = "gex_anndata"
        ANNDATA_EXT = "h5ad"

    rule run_solo:
        """Run SOLO doublet detection on per-capture AnnData object.
        For ARC (MuData), SOLO is run on the GEX (rna) modality."""
        input:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}." + ANNDATA_EXT),
            anndata_done = os.path.join(LOGS_DIR, "{batch}_{capture}_" + ANNDATA_DONE_SUFFIX + ".done")
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
