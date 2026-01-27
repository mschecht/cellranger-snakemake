"""Demultiplexing workflow rules for genetic demultiplexing methods."""

import os
import sys
from pathlib import Path

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger
from cellranger_snakemake.config_validator import parse_output_directories


# Get demultiplexing config
if config.get("demultiplexing"):
    DEMUX_CONFIG = config["demultiplexing"]
    DEMUX_METHOD = DEMUX_CONFIG["method"]
    DEMUX_PARAMS = DEMUX_CONFIG.get("parameters", {})
    DEMUX_SAMPLES = config.get("samples", [])
    
    custom_logger.info(f"Demultiplexing: Using {DEMUX_METHOD} method")


# ============================================================================
# DEMUXALOT
# ============================================================================

if config.get("demultiplexing") and DEMUX_METHOD == "demuxalot":

    rule demuxalot:
        """Run Demuxalot for genetic demultiplexing."""
        input:
            bam = "{sample}/outs/possorted_genome_bam.bam",
            vcf = DEMUX_PARAMS.get("vcf_file")
        output:
            best = "{sample}/demuxalot/{sample}.best",
            done = touch("{sample}/demuxalot/{sample}_demuxalot.done")
        params:
            field = DEMUX_PARAMS.get("field", "GT"),
            group_list = DEMUX_PARAMS.get("group_list", ""),
            alpha = DEMUX_PARAMS.get("alpha", [0.5]),
            out_prefix = "{sample}/demuxalot/{sample}"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/demuxalot/{sample}_demuxalot.log"
        shell:
            """
            popscle demuxalot \\
                --sam {input.bam} \\
                --vcf {input.vcf} \\
                --field {params.field} \\
                --out {params.out_prefix} \\
                2>&1 | tee {log}
            """


# ============================================================================
# VIREO
# ============================================================================

if config.get("demultiplexing") and DEMUX_METHOD == "vireo":
    
    # Parse vireo config
    VIREO_CONFIG = DEMUX_CONFIG.get("vireo", {})
    CELLSNP_CONFIG = VIREO_CONFIG.get("cellsnp", {})
    VIREO_DONORS = VIREO_CONFIG.get("donors")
    
    # cellsnp-lite parameters
    CELLSNP_VCF = CELLSNP_CONFIG.get("vcf")
    CELLSNP_THREADS = CELLSNP_CONFIG.get("threads", 4)
    CELLSNP_MIN_MAF = CELLSNP_CONFIG.get("min_maf", 0.0)
    CELLSNP_MIN_COUNT = CELLSNP_CONFIG.get("min_count", 1)
    CELLSNP_UMI_TAG = CELLSNP_CONFIG.get("umi_tag", "Auto")
    CELLSNP_CELL_TAG = CELLSNP_CONFIG.get("cell_tag", "CB")
    CELLSNP_GZIP = CELLSNP_CONFIG.get("gzip", True)
    
    custom_logger.info(f"Vireo: {VIREO_DONORS} donors, using cellsnp-lite with VCF: {CELLSNP_VCF}")
    
    # Get paths from cellranger GEX output
    if config.get("cellranger_gex"):
        from cellranger_snakemake import utils as u
        gex_cfg_vireo = {
            "reference": config["cellranger_gex"]["reference"],
            "libraries": config["cellranger_gex"]["libraries"],
            "normalize": config["cellranger_gex"].get("normalize", "none"),
            "create_bam": config["cellranger_gex"].get("create-bam", False),
        }
        GEX_LIBRARIES_VIREO = gex_cfg_vireo["libraries"]
        GEX_COUNT_DIR_VIREO = os.path.join(config.get("output_dir", "output"), "01_CELLRANGERGEX_COUNT")
        
        gex_df_vireo = u.sanity_check_libraries_list_tsv(
            GEX_LIBRARIES_VIREO,
            expected_columns={"batch", "capture", "sample", "fastqs"},
            path_column="fastqs",
            file_extension=None
        )
    
    # Get output directory from centralized configuration
    OUTPUT_DIRS = parse_output_directories(config)
    DEMUX_OUTPUT_DIR = OUTPUT_DIRS["demultiplexing_dir"]
    
    rule cellsnp_lite:
        """Run cellsnp-lite for SNP calling from BAM."""
        input:
            gex_done = os.path.join(config.get("output_dir", "output"), "00_LOGS", "{batch}_{capture}_gex_count.done"),
            bam = os.path.join(GEX_COUNT_DIR_VIREO, "{batch}_{capture}", "outs", "possorted_genome_bam.bam"),
            barcodes = os.path.join(GEX_COUNT_DIR_VIREO, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
        output:
            base_vcf = os.path.join(DEMUX_OUTPUT_DIR, "cellsnp_output_{batch}_{capture}", "cellSNP.base.vcf.gz"),
            samples = os.path.join(DEMUX_OUTPUT_DIR, "cellsnp_output_{batch}_{capture}", "cellSNP.samples.tsv"),
            done = touch(os.path.join(config.get("output_dir", "output"), "00_LOGS", "cellsnp_output_{batch}_{capture}.done"))
        params:
            vcf_ref = CELLSNP_VCF,
            outdir = os.path.join(DEMUX_OUTPUT_DIR, "cellsnp_output_{batch}_{capture}"),
            min_maf = CELLSNP_MIN_MAF,
            min_count = CELLSNP_MIN_COUNT,
            umi_tag = CELLSNP_UMI_TAG,
            cell_tag = CELLSNP_CELL_TAG
        threads: CELLSNP_THREADS
        log:
            os.path.join(DEMUX_OUTPUT_DIR, "cellsnp_output_{batch}_{capture}", "cellsnp_lite.log")
        shell:
            """
            cellsnp-lite \\
                -s {input.bam} \\
                -b {input.barcodes} \\
                -O {params.outdir} \\
                -R {params.vcf_ref} \\
                -p {threads} \\
                --minMAF {params.min_maf} \\
                --minCOUNT {params.min_count} \\
                --UMItag {params.umi_tag} \\
                --cellTAG {params.cell_tag} \\
                --gzip \\
                2>&1 | tee {log}
            """
    
    
    rule vireo:
        """Run Vireo for donor deconvolution using cellsnp-lite output."""
        input:
            cellsnp_done = rules.cellsnp_lite.output.done
        output:
            donor_ids = os.path.join(DEMUX_OUTPUT_DIR, "vireo_output_{batch}_{capture}", "donor_ids.tsv"),
            summary = os.path.join(DEMUX_OUTPUT_DIR, "vireo_output_{batch}_{capture}", "summary.tsv"),
            done = touch(os.path.join(OUTPUT_DIRS["logs_dir"], "vireo_output_{batch}_{capture}.done"))
        params:
            n_donor = VIREO_DONORS,
            cellsnp_dir = os.path.join(DEMUX_OUTPUT_DIR, "cellsnp_output_{batch}_{capture}"),
            outdir = os.path.join(DEMUX_OUTPUT_DIR, "vireo_output_{batch}_{capture}")
        threads: 1
        log:
            os.path.join(DEMUX_OUTPUT_DIR, "vireo_output_{batch}_{capture}", "vireo.log")
        shell:
            """
            vireo \\
                -c {params.cellsnp_dir} \\
                -N {params.n_donor} \\
                -o {params.outdir} \\
                2>&1 | tee {log}
            """
