"""Demultiplexing workflow rules for genetic demultiplexing methods."""

import os
import sys
import gzip

import pandas as pd

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))

from pathlib import Path
from custom_logger import custom_logger
from cellranger_snakemake.config_validator import parse_output_directories


# Get demultiplexing config
if config.get("demultiplexing"):
    DEMUX_CONFIG = config["demultiplexing"]
    DEMUX_METHOD = DEMUX_CONFIG["method"]
    DEMUX_PARAMS = DEMUX_CONFIG.get("parameters", {})
    DEMUX_SAMPLES = config.get("samples", [])

    # Shared output directories for all demux methods
    OUTPUT_DIRS = parse_output_directories(config)
    DEMUX_OUTPUT_DIR = OUTPUT_DIRS["demultiplexing_dir"]

    # Determine count directory based on enabled modality
    # Demux works with BAM files from any modality
    if config.get("cellranger_gex"):
        COUNT_DIR = os.path.join(config.get("output_dir", "output"), "01_CELLRANGERGEX_COUNT")
        MODALITY = "gex"
    elif config.get("cellranger_atac"):
        COUNT_DIR = os.path.join(config.get("output_dir", "output"), "01_CELLRANGERATAC_COUNT")
        MODALITY = "atac"
    elif config.get("cellranger_arc"):
        COUNT_DIR = os.path.join(config.get("output_dir", "output"), "01_CELLRANGERARC_COUNT")
        MODALITY = "arc"
    else:
        raise ValueError("Demultiplexing requires cellranger_gex, cellranger_atac, or cellranger_arc to be enabled")

    # Set modality-specific file paths for BAM and barcodes
    if MODALITY == "gex":
        BAM_FILE = "possorted_genome_bam.bam"
        BARCODES_FILE = os.path.join("filtered_feature_bc_matrix", "barcodes.tsv.gz")
    elif MODALITY == "atac":
        BAM_FILE = "possorted_bam.bam"
        BARCODES_FILE = os.path.join("filtered_peak_bc_matrix", "barcodes.tsv.gz")
    elif MODALITY == "arc":
        BAM_FILE = "gex_possorted_bam.bam"
        BARCODES_FILE = os.path.join("filtered_feature_bc_matrix", "barcodes.tsv.gz")

    custom_logger.info(f"Demultiplexing: Using {DEMUX_METHOD} method on {MODALITY.upper()} data")


# ============================================================================
# DEMUXALOT
# ============================================================================

if config.get("demultiplexing") and DEMUX_METHOD == "demuxalot":

    DEMUXALOT_CONFIG = DEMUX_CONFIG.get("demuxalot", {})
    DEMUXALOT_VCF = DEMUXALOT_CONFIG.get("vcf")
    DEMUXALOT_GENOME_NAMES_FILE = DEMUXALOT_CONFIG.get("genome_names")
    DEMUXALOT_REFINE = DEMUXALOT_CONFIG.get("refine")
    DEMUXALOT_CELLTAG = DEMUXALOT_CONFIG.get("celltag")
    DEMUXALOT_UMITAG = DEMUXALOT_CONFIG.get("umitag")

    # Read genome names from file into a list
    genome_names_df = pd.read_csv(DEMUXALOT_GENOME_NAMES_FILE, header=None, names=['genome'])
    DEMUXALOT_GENOME_NAMES_LIST = genome_names_df['genome'].tolist()

    custom_logger.info(f"Demuxalot: Using VCF file: {DEMUXALOT_VCF}")
    custom_logger.info(f"Demuxalot: Using genome names file: {DEMUXALOT_GENOME_NAMES_FILE}")
    custom_logger.info(f"Demuxalot: Genome names: {DEMUXALOT_GENOME_NAMES_LIST}")

    rule demuxalot:
        """Run Demuxalot for genetic demultiplexing."""
        input:
            count_done = os.path.join(config.get("output_dir", "output"), "00_LOGS", "{batch}_{capture}_" + MODALITY + "_count.done"),
            bam = os.path.join(COUNT_DIR, "{batch}_{capture}", "outs", BAM_FILE),
            barcodes = os.path.join(COUNT_DIR, "{batch}_{capture}", "outs", BARCODES_FILE)
        output:
            assign = os.path.join(DEMUX_OUTPUT_DIR, "{batch}_{capture}", "demuxalot", "{batch}_{capture}_assignments.tsv.gz"),
            probs = os.path.join(DEMUX_OUTPUT_DIR, "{batch}_{capture}", "demuxalot", "{batch}_{capture}_posterior_probabilities.tsv.gz"),
            done = touch(os.path.join(config.get("output_dir", "output"), "00_LOGS", "demuxalot_output_{batch}_{capture}.done"))
        params:
            vcf_ref = DEMUXALOT_VCF,
            refine = DEMUXALOT_REFINE or "True",
            celltag = DEMUXALOT_CELLTAG or "CB",
            umitag = DEMUXALOT_UMITAG or "UB",
            out_dir = DEMUX_OUTPUT_DIR,
            sample_names = DEMUXALOT_GENOME_NAMES_LIST
        threads: 4
        log:
            os.path.join(DEMUX_OUTPUT_DIR, "{batch}_{capture}", "demuxalot", "{batch}_{capture}_demuxalot.log"),
        run:
            # below is a simple demuxalot implementation grabbed 
            # from here: https://github.com/arogozhnikov/demuxalot/blob/master/examples/1-plain_demultiplexing.py
            
            from demuxalot import Demultiplexer, BarcodeHandler, ProbabilisticGenotypes, count_snps

            # Load sample names
            genotypes = ProbabilisticGenotypes(genotype_names=DEMUXALOT_GENOME_NAMES_LIST)
            genotypes.add_vcf(DEMUXALOT_VCF)

            print(f'Loaded genotypes: {genotypes}')

            barcode_handler = BarcodeHandler.from_file(input.barcodes)
            print(f'Loaded barcodes: {barcode_handler}')

            snps = count_snps(
                bamfile_location=input.bam,
                chromosome2positions=genotypes.get_chromosome2positions(),
                barcode_handler=barcode_handler,
            )

            print('Collected SNPs: ')
            for chromosome, snps_in_chromosome in snps.items():
                print(f'Chromosome {chromosome}, {snps_in_chromosome.n_snp_calls} calls in {snps_in_chromosome.n_molecules} mols')

            # returns two dataframes with likelihoods and posterior probabilities
            log_likelihoods, posterior_probabilities = Demultiplexer.learn_genotypes(
                snps,
                genotypes=genotypes,
                barcode_handler=barcode_handler,
                doublet_prior=0.25,
            )

            # Assign most probable donor or doublet per cell
            assignments = posterior_probabilities.idxmax(axis=1)

            # Confidence score
            max_probs = posterior_probabilities.max(axis=1)

            # Make result DataFrame
            results_df = pd.DataFrame({
                "assignment": assignments,
                "max_prob": max_probs,
            })

            print('Result:')
            print(results_df)
            results_df.to_csv(output.assign, sep="\t", index=True)
            posterior_probabilities.to_csv(output.probs, sep="\t", index=True)


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
    
    rule cellsnp_lite:
        """Run cellsnp-lite for SNP calling from BAM."""
        input:
            count_done = os.path.join(config.get("output_dir", "output"), "00_LOGS", "{batch}_{capture}_" + MODALITY + "_count.done"),
            bam = os.path.join(COUNT_DIR, "{batch}_{capture}", "outs", BAM_FILE),
            barcodes = os.path.join(COUNT_DIR, "{batch}_{capture}", "outs", BARCODES_FILE)
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
