"""Demultiplexing workflow rules for genetic demultiplexing methods."""

import os
import sys
from pathlib import Path

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger


# Get demultiplexing config
if config.get("demultiplexing"):
    DEMUX_CONFIG = config["demultiplexing"]
    DEMUX_METHOD = DEMUX_CONFIG["method"]
    DEMUX_PARAMS = DEMUX_CONFIG.get("parameters", {})
    DEMUX_SAMPLES = config.get("samples", [])
    
    custom_logger.info(f"Demultiplexing: Using {DEMUX_METHOD} method")


# ============================================================================
# DEMUXLET
# ============================================================================

if config.get("demultiplexing") and DEMUX_METHOD == "demuxlet":
    
    rule demuxlet:
        """Run Demuxlet for genetic demultiplexing."""
        input:
            bam = "{sample}/outs/possorted_genome_bam.bam",
            vcf = DEMUX_PARAMS.get("vcf_file")
        output:
            best = "{sample}/demuxlet/{sample}.best",
            done = touch("{sample}/demuxlet/{sample}_demuxlet.done")
        params:
            field = DEMUX_PARAMS.get("field", "GT"),
            group_list = DEMUX_PARAMS.get("group_list", ""),
            alpha = DEMUX_PARAMS.get("alpha", [0.5]),
            out_prefix = "{sample}/demuxlet/{sample}"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/demuxlet/{sample}_demuxlet.log"
        shell:
            """
            popscle demuxlet \\
                --sam {input.bam} \\
                --vcf {input.vcf} \\
                --field {params.field} \\
                --out {params.out_prefix} \\
                2>&1 | tee {log}
            """


# ============================================================================
# FREEMUXLET
# ============================================================================

if config.get("demultiplexing") and DEMUX_METHOD == "freemuxlet":
    
    rule freemuxlet:
        """Run Freemuxlet for reference-free genetic demultiplexing."""
        input:
            bam = "{sample}/outs/possorted_genome_bam.bam",
            plp = "{sample}/freemuxlet/{sample}.plp.gz"
        output:
            clust = "{sample}/freemuxlet/{sample}.clust1.samples.gz",
            done = touch("{sample}/freemuxlet/{sample}_freemuxlet.done")
        params:
            nsample = DEMUX_PARAMS.get("n_sample"),
            out_prefix = "{sample}/freemuxlet/{sample}"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/freemuxlet/{sample}_freemuxlet.log"
        shell:
            """
            popscle freemuxlet \\
                --plp {input.plp} \\
                --out {params.out_prefix} \\
                --nsample {params.nsample} \\
                2>&1 | tee {log}
            """


# ============================================================================
# SOUPORCELL
# ============================================================================

if config.get("demultiplexing") and DEMUX_METHOD == "souporcell":
    
    rule souporcell:
        """Run Souporcell for demultiplexing and doublet detection."""
        input:
            bam = "{sample}/outs/possorted_genome_bam.bam",
            barcodes = "{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
            fasta = DEMUX_PARAMS.get("fasta")
        output:
            clusters = "{sample}/souporcell/clusters.tsv",
            done = touch("{sample}/souporcell/{sample}_souporcell.done")
        params:
            clusters = DEMUX_PARAMS.get("n_clusters"),
            known_genotypes = DEMUX_PARAMS.get("known_genotypes", ""),
            outdir = "{sample}/souporcell"
        threads: RESOURCES.get("threads", 8)
        log:
            "{sample}/souporcell/{sample}_souporcell.log"
        shell:
            """
            souporcell_pipeline.py \\
                -i {input.bam} \\
                -b {input.barcodes} \\
                -f {input.fasta} \\
                -t {threads} \\
                -o {params.outdir} \\
                -k {params.clusters} \\
                2>&1 | tee {log}
            """


# ============================================================================
# VIREO
# ============================================================================

if config.get("demultiplexing") and DEMUX_METHOD == "vireo":
    
    rule vireo:
        """Run Vireo for donor deconvolution."""
        input:
            celldata = "{sample}/vireo/cellSNP.cells.vcf.gz",
            donors_vcf = DEMUX_PARAMS.get("donor_vcf", "")
        output:
            donor_ids = "{sample}/vireo/donor_ids.tsv",
            done = touch("{sample}/vireo/{sample}_vireo.done")
        params:
            n_donor = DEMUX_PARAMS.get("n_donor"),
            outdir = "{sample}/vireo"
        threads: RESOURCES.get("threads", 4)
        log:
            "{sample}/vireo/{sample}_vireo.log"
        shell:
            """
            vireo \\
                -c {input.celldata} \\
                -N {params.n_donor} \\
                -o {params.outdir} \\
                2>&1 | tee {log}
            """
