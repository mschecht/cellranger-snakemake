"""Cell Ranger workflow rules for all modalities (GEX, ATAC, ARC)."""

# Standard library imports
import os
import sys
import shutil

import pandas as pd

from pathlib import Path
from tempfile import gettempdir
from collections import defaultdict
from cellranger_snakemake import utils as u

def parse_cellranger_config(config, modality_key, has_chemistry=True):
    """
    Parse Cell Ranger configuration for a specific modality.
    
    Args:
        config: Snakemake config dictionary
        modality_key: Key for the modality (e.g., "cellranger_gex")
        has_chemistry: Whether this modality supports chemistry parameter
        
    Returns:
        dict: Parsed configuration with standardized keys
    """
    mod_config = config[modality_key]
    output_dir = config.get("output_dir", "output")
    
    # Get directories with defaults
    dirs = mod_config.get("directories", {})
    
    # Determine directory names based on modality
    if "gex" in modality_key:
        prefix = "CELLRANGERGEX"
    elif "atac" in modality_key:
        prefix = "CELLRANGERATAC"
    elif "arc" in modality_key:
        prefix = "CELLRANGERARC"
    
    parsed = {
        "reference": mod_config["reference"],
        "libraries": mod_config["libraries"],
        "normalize": mod_config.get("normalize", "none"),
        "create_bam": mod_config.get("create-bam", False),
        "logs_dir": os.path.join(output_dir, dirs.get("LOGS_DIR", "00_LOGS")),
        "count_dir": os.path.join(output_dir, dirs.get(f"{prefix}_COUNT_DIR", f"01_{prefix}_COUNT")),
        "aggr_dir": os.path.join(output_dir, dirs.get(f"{prefix}_AGGR_DIR", f"02_{prefix}_AGGR")),
    }
    
    if has_chemistry:
        parsed["chemistry"] = mod_config.get("chemistry", "auto")
    
    return parsed


# ============================================================================
# CELL RANGER GEX
# ============================================================================

if config.get("cellranger_gex"):
    gex_cfg = parse_cellranger_config(config, "cellranger_gex", has_chemistry=True)
    
    GEX_REFERENCE = gex_cfg["reference"]
    GEX_LIBRARIES = gex_cfg["libraries"]
    GEX_CHEMISTRY = gex_cfg["chemistry"]
    GEX_CREATE_BAM = gex_cfg["create_bam"]
    GEX_NORMALIZE = gex_cfg["normalize"]
    GEX_LOGS_DIR = gex_cfg["logs_dir"]
    GEX_COUNT_DIR = gex_cfg["count_dir"]
    GEX_AGGR_DIR = gex_cfg["aggr_dir"]
    
    # Parse libraries file
    gex_df = u.sanity_check_libraries_list_tsv(
        GEX_LIBRARIES,
        expected_columns={"batch", "capture", "sample", "fastqs"},
        path_column="fastqs",
        file_extension=None
    )
    
    # Build sample mappings - convert batch to string for consistency
    gex_samples = gex_df["capture"].unique().tolist()
    gex_batch_to_samples = {str(k): v for k, v in gex_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    u.custom_logger.info(f"Cell Ranger GEX: Found {len(gex_samples)} sample(s) across {len(gex_batch_to_samples)} batch(es)")
    
    
    rule cellranger_gex_count:
        """Run Cell Ranger count for gene expression data."""
        input:
            reference = GEX_REFERENCE,
            fastqs = lambda wc: gex_df[gex_df["capture"] == wc.capture]["fastqs"].iloc[0]
        output:
            h5 = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
            summary = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "web_summary.html"),
            done = touch(os.path.join(GEX_LOGS_DIR, "{batch}_{capture}_gex_count.done"))
        params:
            outdir = GEX_COUNT_DIR,
            chemistry = GEX_CHEMISTRY,
            create_bam = GEX_CREATE_BAM,
            sample_name = lambda wc: gex_df[gex_df["capture"] == wc.capture]["sample"].iloc[0]
        threads: config["cellranger_gex"].get("threads", 10)
        resources:
            mem_mb = config["cellranger_gex"].get("mem_gb", 64) * 1024,
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(GEX_LOGS_DIR, "{batch}_{capture}_gex_count.log")
        run:
            output_id = f"{wildcards.batch}_{wildcards.capture}"
            create_bam_str = str(params.create_bam).lower() # Convert Python boolean to lowercase string for Cell Ranger

            shell(
                f"""
                cellranger count \\
                    --id={output_id} \\
                    --transcriptome={input.reference} \\
                    --fastqs={input.fastqs} \\
                    --sample={params.sample_name} \\
                    --create-bam={create_bam_str} \\
                    --chemistry={params.chemistry} \\
                    2>&1 > {log}
                """
            )
            # Move output to final location
            if os.path.exists(output_id):
                final_path = os.path.join(params.outdir, output_id)
                if os.path.exists(final_path):
                    shutil.rmtree(final_path)
                shutil.move(output_id, final_path)
    
    
    rule cellranger_gex_aggr:
        """Aggregate multiple GEX samples."""
        input:
            samples = lambda wc: expand(
                os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
                batch=wc.batch,
                capture=gex_batch_to_samples[wc.batch]
            )
        output:
            h5 = os.path.join(GEX_AGGR_DIR, "{batch}", "outs", "count", "filtered_feature_bc_matrix.h5"),
            done = touch(os.path.join(GEX_LOGS_DIR, "{batch}_gex_aggr.done"))
        params:
            outdir = GEX_AGGR_DIR,
            normalize = GEX_NORMALIZE,
            csv = os.path.join(GEX_AGGR_DIR, "{batch}_aggregation.csv")
        threads: config["cellranger_gex"].get("threads", 10)
        resources:
            mem_mb = config["cellranger_gex"].get("mem_gb", 64) * 1024, # NOTE: SLURM resources uses M  B
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(GEX_LOGS_DIR, "{batch}_gex_aggr.log")
        run:
            # Create aggregation CSV
            captures = gex_batch_to_samples[wildcards.batch]
            aggr_data = []
            for capture in captures:
                molecule_h5 = os.path.join(
                    GEX_COUNT_DIR,
                    f"{wildcards.batch}_{capture}",
                    "outs",
                    "molecule_info.h5"
                )
                aggr_data.append({"sample_id": f"{wildcards.batch}_{capture}", "molecule_h5": molecule_h5})
            
            pd.DataFrame(aggr_data).to_csv(params.csv, index=False)
            
            # Only run aggregation if more than one capture
            if len(captures) > 1:
                shell(
                    f"""
                    cellranger aggr \\
                        --id={wildcards.batch} \\
                        --csv={params.csv} \\
                        --normalize={params.normalize} \\
                        2>&1 > {log}
                    """
                )
                # Move output to final location
                if os.path.exists(wildcards.batch):
                    final_path = os.path.join(params.outdir, wildcards.batch)
                    if os.path.exists(final_path):
                        shutil.rmtree(final_path)
                    shutil.move(wildcards.batch, final_path)
            else:
                u.custom_logger.info(f"Batch {wildcards.batch} has only one capture ({captures[0]}). Skipping cellranger aggr step.")


# ============================================================================
# CELL RANGER ATAC
# ============================================================================

if config.get("cellranger_atac"):
    atac_cfg = parse_cellranger_config(config, "cellranger_atac", has_chemistry=True)
    
    ATAC_REFERENCE = atac_cfg["reference"]
    ATAC_LIBRARIES = atac_cfg["libraries"]
    ATAC_CHEMISTRY = atac_cfg["chemistry"]
    ATAC_NORMALIZE = atac_cfg["normalize"]
    ATAC_LOGS_DIR = atac_cfg["logs_dir"]
    ATAC_COUNT_DIR = atac_cfg["count_dir"]
    ATAC_AGGR_DIR = atac_cfg["aggr_dir"]
    
    # Parse libraries file
    atac_df = u.sanity_check_libraries_list_tsv(
        ATAC_LIBRARIES,
        expected_columns={"batch", "capture", "sample", "fastqs"},
        path_column="fastqs",
        file_extension=None
    )
    
    atac_samples = atac_df["capture"].unique().tolist()
    atac_batch_to_samples = {str(k): v for k, v in atac_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    u.custom_logger.info(f"Cell Ranger ATAC: Found {len(atac_samples)} sample(s) across {len(atac_batch_to_samples)} batch(es)")
    
    
    rule cellranger_atac_count:
        """Run Cell Ranger ATAC count."""
        input:
            reference = ATAC_REFERENCE,
            fastqs = lambda wc: atac_df[atac_df["capture"] == wc.capture]["fastqs"].iloc[0]
        output:
            h5 = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_peak_bc_matrix.h5"),
            summary = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "web_summary.html"),
            done = touch(os.path.join(ATAC_LOGS_DIR, "{batch}_{capture}_atac_count.done"))
        params:
            outdir = ATAC_COUNT_DIR,
            sample_name = lambda wc: atac_df[atac_df["capture"] == wc.capture]["sample"].iloc[0]
        threads: 8
        resources:
            mem_gb = RESOURCES.get("mem_gb", 64),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ATAC_LOGS_DIR, "{batch}_{capture}_atac_count.log")
        run:
            output_id = f"{wildcards.batch}_{wildcards.capture}"
            shell(
                f"""
                cellranger-atac count \\
                    --id={output_id} \\
                    --reference={input.reference} \\
                    --fastqs={input.fastqs} \\
                    --sample={params.sample_name} \\
                    2>&1 > {log}
                """
            )
            # Move output to final location
            if os.path.exists(output_id):
                final_path = os.path.join(params.outdir, output_id)
                if os.path.exists(final_path):
                    shutil.rmtree(final_path)
                shutil.move(output_id, final_path)
    
    
    rule cellranger_atac_aggr:
        """Aggregate multiple ATAC samples."""
        input:
            reference = ATAC_REFERENCE,
            samples = lambda wc: expand(
                os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_peak_bc_matrix.h5"),
                batch=wc.batch,
                capture=atac_batch_to_samples[wc.batch]
            )
        output:
            h5 = os.path.join(ATAC_AGGR_DIR, "{batch}", "outs", "filtered_peak_bc_matrix.h5"),
            done = touch(os.path.join(ATAC_LOGS_DIR, "{batch}_atac_aggr.done"))
        params:
            outdir = ATAC_AGGR_DIR,
            normalize = ATAC_NORMALIZE,
            csv = os.path.join(ATAC_AGGR_DIR, "{batch}_aggregation.csv")
        threads: 8
        resources:
            mem_gb = RESOURCES.get("mem_gb", 64),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ATAC_LOGS_DIR, "{batch}_atac_aggr.log")
        run:
            # Create aggregation CSV
            captures = atac_batch_to_samples[wildcards.batch]
            aggr_data = []
            for capture in captures:
                fragments = os.path.abspath(os.path.join(
                    ATAC_COUNT_DIR,
                    f"{wildcards.batch}_{capture}",
                    "outs",
                    "fragments.tsv.gz"
                ))
                singlecell = os.path.abspath(os.path.join(
                    ATAC_COUNT_DIR,
                    f"{wildcards.batch}_{capture}",
                    "outs",
                    "singlecell.csv"
                ))
                aggr_data.append({"library_id": f"{wildcards.batch}_{capture}", "fragments": fragments, "cells": singlecell})
            
            pd.DataFrame(aggr_data).to_csv(params.csv, index=False)
            
            # Only run aggregation if more than one capture
            if len(captures) > 1:
                shell(
                    f"""
                    cellranger-atac aggr \\
                        --id={wildcards.batch} \\
                        --csv={params.csv} \\
                        --reference={input.reference} \\
                        --normalize={params.normalize} \\
                        2>&1 > {log}
                    """
                )
                # Move output to final location
                if os.path.exists(wildcards.batch):
                    final_path = os.path.join(params.outdir, wildcards.batch)
                    if os.path.exists(final_path):
                        shutil.rmtree(final_path)
                    shutil.move(wildcards.batch, final_path)
            else:
                u.custom_logger.info(f"Batch {wildcards.batch} has only one capture ({captures[0]}). Skipping cellranger-atac aggr step.")


# ============================================================================
# CELL RANGER ARC
# ============================================================================

if config.get("cellranger_arc"):
    arc_cfg = parse_cellranger_config(config, "cellranger_arc", has_chemistry=False)
    
    ARC_REFERENCE = arc_cfg["reference"]
    ARC_LIBRARIES = arc_cfg["libraries"]
    ARC_NORMALIZE = arc_cfg["normalize"]
    ARC_LOGS_DIR = arc_cfg["logs_dir"]
    ARC_COUNT_DIR = arc_cfg["count_dir"]
    ARC_AGGR_DIR = arc_cfg["aggr_dir"]
    
    # Parse libraries file (ARC uses CSV column)
    arc_df = u.sanity_check_libraries_list_tsv(
        ARC_LIBRARIES,
        expected_columns={"batch", "capture", "CSV"},
        path_column="CSV",
        file_extension=".csv"
    )
    
    arc_captures = arc_df["capture"].unique().tolist()
    arc_batch_to_captures = {str(k): v for k, v in arc_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    u.custom_logger.info(f"Cell Ranger ARC: Found {len(arc_captures)} capture(s) across {len(arc_batch_to_captures)} batch(es)")
    
    
    rule cellranger_arc_count:
        """Run Cell Ranger ARC count for multiome data."""
        input:
            reference = ARC_REFERENCE,
            libraries_csv = lambda wc: arc_df[arc_df["capture"] == wc.capture]["CSV"].iloc[0]
        output:
            h5 = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
            summary = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "web_summary.html"),
            done = touch(os.path.join(ARC_LOGS_DIR, "{batch}_{capture}_arc_count.done"))
        params:
            outdir = ARC_COUNT_DIR
        threads: 8
        resources:
            mem_gb = RESOURCES.get("mem_gb", 64),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ARC_LOGS_DIR, "{batch}_{capture}_arc_count.log")
        run:
            output_id = f"{wildcards.batch}_{wildcards.capture}"
            shell(
                f"""
                cellranger-arc count \\
                    --id={output_id} \\
                    --reference={input.reference} \\
                    --libraries={input.libraries_csv} \\
                    2>&1 > {log}
                """
            )
            # Move output to final location
            if os.path.exists(output_id):
                final_path = os.path.join(params.outdir, output_id)
                if os.path.exists(final_path):
                    shutil.rmtree(final_path)
                shutil.move(output_id, final_path)
    
    
    rule cellranger_arc_aggr:
        """Aggregate multiple ARC captures."""
        input:
            captures = lambda wc: expand(
                os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
                batch=wc.batch,
                capture=arc_batch_to_captures[wc.batch]
            )
        output:
            h5 = os.path.join(ARC_AGGR_DIR, "{batch}", "outs", "filtered_feature_bc_matrix.h5"),
            done = touch(os.path.join(ARC_LOGS_DIR, "{batch}_arc_aggr.done"))
        params:
            outdir = ARC_AGGR_DIR,
            normalize = ARC_NORMALIZE,
            csv = os.path.join(ARC_AGGR_DIR, "{batch}_aggregation.csv"),
            reference = ARC_REFERENCE
        threads: 8
        resources:
            mem_gb = RESOURCES.get("mem_gb", 64),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(ARC_LOGS_DIR, "{batch}_arc_aggr.log")
        run:
            # Create aggregation CSV for cellranger-arc aggr
            captures = arc_batch_to_captures[wildcards.batch]
            aggr_data = []
            for capture in captures:
                library_id = f"{wildcards.batch}_{capture}"
                atac_fragments = os.path.abspath(os.path.join(
                    ARC_COUNT_DIR,
                    library_id,
                    "outs",
                    "atac_fragments.tsv.gz"
                ))
                per_barcode_metrics = os.path.abspath(os.path.join(
                    ARC_COUNT_DIR,
                    library_id,
                    "outs",
                    "per_barcode_metrics.csv"
                ))
                gex_molecule_info = os.path.abspath(os.path.join(
                    ARC_COUNT_DIR,
                    library_id,
                    "outs",
                    "gex_molecule_info.h5"
                ))
                aggr_data.append({
                    "library_id": library_id,
                    "atac_fragments": atac_fragments,
                    "per_barcode_metrics": per_barcode_metrics,
                    "gex_molecule_info": gex_molecule_info
                })
            
            pd.DataFrame(aggr_data).to_csv(params.csv, index=False)
            
            # Only run aggregation if more than one capture
            if len(captures) > 1:
                shell(
                    f"""
                    cellranger-arc aggr \\
                        --id={wildcards.batch} \\
                        --csv={params.csv} \\
                        --reference={params.reference} \\
                        --normalize={params.normalize} \\
                        2>&1 > {log}
                    """
                )
                # Move output to final location
                if os.path.exists(wildcards.batch):
                    final_path = os.path.join(params.outdir, wildcards.batch)
                    if os.path.exists(final_path):
                        shutil.rmtree(final_path)
                    shutil.move(wildcards.batch, final_path)
            else:
                u.custom_logger.info(f"Batch {wildcards.batch} has only one capture ({captures[0]}). Skipping cellranger-arc aggr step.")
