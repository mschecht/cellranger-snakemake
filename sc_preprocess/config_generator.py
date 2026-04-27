"""Config generator for single-cell preprocessing pipeline."""

_GEX_TEMPLATE = """\
project_name: my_project  # REQUIRED
output_dir: output  # REQUIRED

resources:
  mem_gb: 32
  tmpdir: ""  # temp directory for large file operations

directories_suffix: none  # suffix appended to output directory names; "none" to disable

cellranger_gex:
  enabled: true
  reference: /path/to/cellranger/reference  # REQUIRED
  libraries: /path/to/libraries.tsv  # REQUIRED: TSV with columns: batch, capture, sample, fastqs
  chemistry: auto  # options: auto, SC3Pv1, SC3Pv2, SC3Pv3, fiveprime, SC5P-PE, SC5P-R2, ARC-v1
  normalize: none  # options: none, mapped, depth
  create-bam: false
  threads: 10
  mem_gb: 64
  runtime_minutes: 720
  jobmode: null  # options: local, slurm, sge, lsf, or path to .template file
  mempercore: null  # cluster only: GB RAM per core (--mempercore)
  maxjobs: null  # cluster only: max concurrent jobs (--maxjobs)
  jobinterval: null  # cluster only: ms between job submissions (--jobinterval)
  anndata_threads: 1
  anndata_mem_gb: 16

doublet_detection:
  enabled: false
  method: scrublet  # only supported method
  threads: 1
  mem_gb: 16
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    threshold: null  # doublet score cutoff; auto-determined if null
    n_neighbors: null  # KNN neighbors; auto-set to 0.5*sqrt(n_obs) if null
    random_state: 0

demultiplexing:
  enabled: false
  method: demuxalot  # options: demuxalot, vireo — only the selected method is used
  demuxalot:  # genotype-based; GEX only
    vcf: /path/to/genotypes.vcf  # REQUIRED
    genome_names: /path/to/genome_names.txt  # REQUIRED
    refine: false  # run genotype refinement step
    celltag: CB
    umitag: UB
  vireo:  # SNP-based; works with all modalities
    donors: 4  # REQUIRED: number of donors
    cellsnp:
      vcf: /path/to/variants.vcf  # REQUIRED
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
"""

_ATAC_TEMPLATE = """\
project_name: my_project  # REQUIRED
output_dir: output  # REQUIRED

resources:
  mem_gb: 32
  tmpdir: ""  # temp directory for large file operations

directories_suffix: none  # suffix appended to output directory names; "none" to disable

cellranger_atac:
  enabled: true
  reference: /path/to/cellranger-atac/reference  # REQUIRED
  libraries: /path/to/libraries.tsv  # REQUIRED: TSV with columns: batch, capture, sample, fastqs
  chemistry: auto  # options: auto, ARC-v1
  normalize: none  # options: none, depth
  threads: 10
  mem_gb: 64
  runtime_minutes: 720
  jobmode: null  # options: local, slurm, sge, lsf, or path to .template file
  mempercore: null  # cluster only: GB RAM per core (--mempercore)
  maxjobs: null  # cluster only: max concurrent jobs (--maxjobs)
  jobinterval: null  # cluster only: ms between job submissions (--jobinterval)
  anndata_threads: 1
  anndata_mem_gb: 32  # SnapATAC2 fragment sorting requires extra memory

doublet_detection:
  enabled: false
  method: scrublet  # only supported method
  threads: 1
  mem_gb: 16
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    threshold: null  # doublet score cutoff; auto-determined if null
    n_neighbors: null  # KNN neighbors; auto-set to 0.5*sqrt(n_obs) if null
    random_state: 0

demultiplexing:
  enabled: false
  method: vireo  # options: demuxalot, vireo — only the selected method is used
  demuxalot:  # genotype-based; GEX only
    vcf: /path/to/genotypes.vcf  # REQUIRED
    genome_names: /path/to/genome_names.txt  # REQUIRED
    refine: false  # run genotype refinement step
    celltag: CB
    umitag: UB
  vireo:  # SNP-based; works with all modalities
    donors: 4  # REQUIRED: number of donors
    cellsnp:
      vcf: /path/to/variants.vcf  # REQUIRED
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
"""

_ARC_TEMPLATE = """\
project_name: my_project  # REQUIRED
output_dir: output  # REQUIRED

resources:
  mem_gb: 32
  tmpdir: ""  # temp directory for large file operations

directories_suffix: none  # suffix appended to output directory names; "none" to disable

cellranger_arc:
  enabled: true
  reference: /path/to/cellranger-arc/reference  # REQUIRED
  libraries: /path/to/libraries.tsv  # REQUIRED: TSV with columns: batch, capture, CSV (per-capture CSV path)
  normalize: none  # options: none, depth
  threads: 10
  mem_gb: 64
  runtime_minutes: 720
  jobmode: null  # options: local, slurm, sge, lsf, or path to .template file
  mempercore: null  # cluster only: GB RAM per core (--mempercore)
  maxjobs: null  # cluster only: max concurrent jobs (--maxjobs)
  jobinterval: null  # cluster only: ms between job submissions (--jobinterval)
  anndata_threads: 1
  anndata_mem_gb: 16

doublet_detection:
  enabled: false
  method: scrublet  # only supported method
  threads: 1
  mem_gb: 16
  scrublet:
    filter_cells_min_genes: 100
    filter_genes_min_cells: 3
    expected_doublet_rate: 0.06
    min_gene_variability_pctl: 85.0
    n_prin_comps: 30
    sim_doublet_ratio: 2.0
    threshold: null  # doublet score cutoff; auto-determined if null
    n_neighbors: null  # KNN neighbors; auto-set to 0.5*sqrt(n_obs) if null
    random_state: 0

demultiplexing:
  enabled: false
  method: vireo  # options: demuxalot, vireo — only the selected method is used
  demuxalot:  # genotype-based; GEX only
    vcf: /path/to/genotypes.vcf  # REQUIRED
    genome_names: /path/to/genome_names.txt  # REQUIRED
    refine: false  # run genotype refinement step
    celltag: CB
    umitag: UB
  vireo:  # SNP-based; works with all modalities
    donors: 4  # REQUIRED: number of donors
    cellsnp:
      vcf: /path/to/variants.vcf  # REQUIRED
      threads: 4
      min_maf: 0.0
      min_count: 1
      umi_tag: Auto
      cell_tag: CB
      gzip: true
"""

_TEMPLATES = {
    "gex": _GEX_TEMPLATE,
    "atac": _ATAC_TEMPLATE,
    "arc": _ARC_TEMPLATE,
}


def generate_config_yaml(modality: str) -> str:
    """Return a fully-commented YAML config string for the given modality."""
    return _TEMPLATES[modality]
