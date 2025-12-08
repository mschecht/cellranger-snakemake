"""Cell Ranger workflow rules for all modalities (GEX, ATAC, ARC)."""

import os
import sys
import pandas as pd
from pathlib import Path
from tempfile import gettempdir
from collections import defaultdict

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from utils import sanity_check_libraries_list_tsv, get_directories
from custom_logger import custom_logger


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
    output_dir = config.get("output_dir", "./output")
    
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
    GEX_NORMALIZE = gex_cfg["normalize"]
    GEX_LOGS_DIR = gex_cfg["logs_dir"]
    GEX_COUNT_DIR = gex_cfg["count_dir"]
    GEX_AGGR_DIR = gex_cfg["aggr_dir"]
    
    # Parse libraries file
    gex_df = sanity_check_libraries_list_tsv(
        GEX_LIBRARIES,
        expected_columns={"batch", "capture", "sample", "fastqs"},
        path_column="fastqs",
        file_extension=None
    )
    
    # Build sample mappings - convert batch to string for consistency
    gex_samples = gex_df["capture"].unique().tolist()
    gex_batch_to_samples = {str(k): v for k, v in gex_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    custom_logger.info(f"Cell Ranger GEX: Found {len(gex_samples)} sample(s) across {len(gex_batch_to_samples)} batch(es)")
    
    
    rule cellranger_gex_count:
        """Run Cell Ranger count for gene expression data."""
        input:
            reference = GEX_REFERENCE,
            fastqs = lambda wc: gex_df[gex_df["capture"] == wc.sample]["fastqs"].iloc[0]
        output:
            h5 = os.path.join(GEX_COUNT_DIR, "{sample}/outs/filtered_feature_bc_matrix.h5"),
            summary = os.path.join(GEX_COUNT_DIR, "{sample}/outs/web_summary.html"),
            done = touch(os.path.join(GEX_LOGS_DIR, "{sample}_gex_count.done"))
        params:
            outdir = GEX_COUNT_DIR,
            chemistry = GEX_CHEMISTRY,
            mem_gb = RESOURCES.get("mem_gb", 64),
            sample_name = lambda wc: gex_df[gex_df["capture"] == wc.sample]["sample"].iloc[0]
        threads: 8
        log:
            os.path.join(GEX_LOGS_DIR, "{sample}_gex_count.log")
        shell:
            """
            cd {params.outdir}
            cellranger count \\
                --id={wildcards.sample} \\
                --transcriptome={input.reference} \\
                --fastqs={input.fastqs} \\
                --sample={params.sample_name} \\
                --chemistry={params.chemistry} \\
                --localcores={threads} \\
                --localmem={params.mem_gb} \\
                2>&1 | tee {log}
            """
    
    
    rule cellranger_gex_aggr:
        """Aggregate multiple GEX samples."""
        input:
            samples = lambda wc: expand(
                os.path.join(GEX_COUNT_DIR, "{sample}/outs/filtered_feature_bc_matrix.h5"),
                sample=gex_batch_to_samples[wc.batch]
            )
        output:
            h5 = os.path.join(GEX_AGGR_DIR, "{batch}/outs/count/filtered_feature_bc_matrix.h5"),
            done = touch(os.path.join(GEX_LOGS_DIR, "{batch}_gex_aggr.done"))
        params:
            outdir = GEX_AGGR_DIR,
            normalize = GEX_NORMALIZE,
            mem_gb = RESOURCES.get("mem_gb", 64),
            csv = os.path.join(GEX_AGGR_DIR, "{batch}_aggregation.csv")
        threads: 8
        log:
            os.path.join(GEX_LOGS_DIR, "{batch}_gex_aggr.log")
        run:
            # Create aggregation CSV
            samples = gex_batch_to_samples[wildcards.batch]
            aggr_data = []
            for sample in samples:
                molecule_h5 = os.path.join(
                    GEX_COUNT_DIR,
                    f"{sample}/outs/molecule_info.h5"
                )
                aggr_data.append({"sample_id": sample, "molecule_h5": molecule_h5})
            
            pd.DataFrame(aggr_data).to_csv(params.csv, index=False)
            
            # Only run aggregation if more than one sample
            if len(samples) > 1:
                shell(
                    """
                    cd {params.outdir}
                    cellranger aggr \\
                        --id={wildcards.batch} \\
                        --csv={params.csv} \\
                        --normalize={params.normalize} \\
                        --localcores={threads} \\
                        --localmem={params.mem_gb} \\
                        2>&1 | tee {log}
                    """
                )
            else:
                custom_logger.info(f"Batch {wildcards.batch} has only one sample ({samples[0]}). Skipping cellranger aggr step.")


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
    atac_df = sanity_check_libraries_list_tsv(
        ATAC_LIBRARIES,
        expected_columns={"batch", "capture", "sample", "fastqs"},
        path_column="fastqs",
        file_extension=None
    )
    
    atac_samples = atac_df["capture"].unique().tolist()
    atac_batch_to_samples = {str(k): v for k, v in atac_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    custom_logger.info(f"Cell Ranger ATAC: Found {len(atac_samples)} sample(s) across {len(atac_batch_to_samples)} batch(es)")
    
    
    rule cellranger_atac_count:
        """Run Cell Ranger ATAC count."""
        input:
            reference = ATAC_REFERENCE,
            fastqs = lambda wc: atac_df[atac_df["capture"] == wc.sample]["fastqs"].iloc[0]
        output:
            h5 = os.path.join(ATAC_COUNT_DIR, "{sample}/outs/filtered_peak_bc_matrix.h5"),
            summary = os.path.join(ATAC_COUNT_DIR, "{sample}/outs/web_summary.html"),
            done = touch(os.path.join(ATAC_LOGS_DIR, "{sample}_atac_count.done"))
        params:
            outdir = ATAC_COUNT_DIR,
            mem_gb = RESOURCES.get("mem_gb", 64),
            sample_name = lambda wc: atac_df[atac_df["capture"] == wc.sample]["sample"].iloc[0]
        threads: 8
        log:
            os.path.join(ATAC_LOGS_DIR, "{sample}_atac_count.log")
        shell:
            """
            cd {params.outdir}
            cellranger-atac count \\
                --id={wildcards.sample} \\
                --reference={input.reference} \\
                --fastqs={input.fastqs} \\
                --sample={params.sample_name} \\
                --localcores={threads} \\
                --localmem={params.mem_gb} \\
                2>&1 | tee {log}
            """
    
    
    rule cellranger_atac_aggr:
        """Aggregate multiple ATAC samples."""
        input:
            samples = lambda wc: expand(
                os.path.join(ATAC_COUNT_DIR, "{sample}/outs/filtered_peak_bc_matrix.h5"),
                sample=atac_batch_to_samples[wc.batch]
            )
        output:
            h5 = os.path.join(ATAC_AGGR_DIR, "{batch}/outs/filtered_peak_bc_matrix.h5"),
            done = touch(os.path.join(ATAC_LOGS_DIR, "{batch}_atac_aggr.done"))
        params:
            outdir = ATAC_AGGR_DIR,
            normalize = ATAC_NORMALIZE,
            mem_gb = RESOURCES.get("mem_gb", 64),
            csv = os.path.join(ATAC_AGGR_DIR, "{batch}_aggregation.csv")
        threads: 8
        log:
            os.path.join(ATAC_LOGS_DIR, "{batch}_atac_aggr.log")
        run:
            # Create aggregation CSV
            samples = atac_batch_to_samples[wildcards.batch]
            aggr_data = []
            for sample in samples:
                fragments = os.path.join(
                    ATAC_COUNT_DIR,
                    f"{sample}/outs/fragments.tsv.gz"
                )
                aggr_data.append({"library_id": sample, "fragments": fragments, "cells": None})
            
            pd.DataFrame(aggr_data).to_csv(params.csv, index=False)
            
            # Only run aggregation if more than one sample
            if len(samples) > 1:
                shell(
                    """
                    cd {params.outdir}
                    cellranger-atac aggr \\
                        --id={wildcards.batch} \\
                        --csv={params.csv} \\
                        --normalize={params.normalize} \\
                        --localcores={threads} \\
                        --localmem={params.mem_gb} \\
                        2>&1 | tee {log}
                    """
                )
            else:
                custom_logger.info(f"Batch {wildcards.batch} has only one sample ({samples[0]}). Skipping cellranger-atac aggr step.")


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
    arc_df = sanity_check_libraries_list_tsv(
        ARC_LIBRARIES,
        expected_columns={"batch", "capture", "CSV"},
        path_column="CSV",
        file_extension=".csv"
    )
    
    arc_captures = arc_df["capture"].unique().tolist()
    arc_batch_to_captures = {str(k): v for k, v in arc_df.groupby("batch")["capture"].apply(list).to_dict().items()}
    
    custom_logger.info(f"Cell Ranger ARC: Found {len(arc_captures)} capture(s) across {len(arc_batch_to_captures)} batch(es)")
    
    
    rule cellranger_arc_count:
        """Run Cell Ranger ARC count for multiome data."""
        input:
            reference = ARC_REFERENCE,
            library_csv = lambda wc: arc_df[arc_df["capture"] == wc.capture]["CSV"].iloc[0]
        output:
            h5 = os.path.join(ARC_COUNT_DIR, "{capture}/outs/filtered_feature_bc_matrix.h5"),
            summary = os.path.join(ARC_COUNT_DIR, "{capture}/outs/web_summary.html"),
            done = touch(os.path.join(ARC_LOGS_DIR, "{capture}_arc_count.done"))
        params:
            outdir = ARC_COUNT_DIR,
            mem_gb = RESOURCES.get("mem_gb", 64)
        threads: 8
        log:
            os.path.join(ARC_LOGS_DIR, "{capture}_arc_count.log")
        shell:
            """
            cd {params.outdir}
            cellranger-arc count \\
                --id={wildcards.capture} \\
                --reference={input.reference} \\
                --libraries={input.library_csv} \\
                --localcores={threads} \\
                --localmem={params.mem_gb} \\
                2>&1 | tee {log}
            """
    
    
    rule cellranger_arc_aggr:
        """Aggregate multiple ARC captures."""
        input:
            captures = lambda wc: expand(
                os.path.join(ARC_COUNT_DIR, "{capture}/outs/filtered_feature_bc_matrix.h5"),
                capture=arc_batch_to_captures[wc.batch]
            )
        output:
            h5 = os.path.join(ARC_AGGR_DIR, "{batch}/outs/count/filtered_feature_bc_matrix.h5"),
            done = touch(os.path.join(ARC_LOGS_DIR, "{batch}_arc_aggr.done"))
        params:
            outdir = ARC_AGGR_DIR,
            normalize = ARC_NORMALIZE,
            mem_gb = RESOURCES.get("mem_gb", 64),
            csv = os.path.join(ARC_AGGR_DIR, "{batch}_aggregation.csv")
        threads: 8
        log:
            os.path.join(ARC_LOGS_DIR, "{batch}_arc_aggr.log")
        run:
            # Create aggregation CSV
            captures = arc_batch_to_captures[wildcards.batch]
            aggr_data = []
            for capture in captures:
                count_dir = os.path.join(ARC_COUNT_DIR, f"{capture}/outs")
                aggr_data.append({"library_id": capture, "atac_fragments": f"{count_dir}/atac_fragments.tsv.gz", 
                                 "per_barcode_metrics": f"{count_dir}/per_barcode_metrics.csv"})
            
            pd.DataFrame(aggr_data).to_csv(params.csv, index=False)
            
            # Only run aggregation if more than one capture
            if len(captures) > 1:
                shell(
                    """
                    cd {params.outdir}
                    cellranger-arc aggr \\
                        --id={wildcards.batch} \\
                        --csv={params.csv} \\
                        --normalize={params.normalize} \\
                        --localcores={threads} \\
                        --localmem={params.mem_gb} \\
                        2>&1 | tee {log}
                    """
                )
            else:
                custom_logger.info(f"Batch {wildcards.batch} has only one capture ({captures[0]}). Skipping cellranger-arc aggr step.")
