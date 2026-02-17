# TLDR

The tool `cellranger-snakemake` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) wrapper for [10x Cell Ranger workflows](https://www.10xgenomics.com/software) supporting [GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest), [ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), and [ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) single-cell data processing.

# Description

Reproducibility and scalability are essential components of contemporary [FAIR](https://www.nature.com/articles/sdata201618) (Findable, Accessible, Interoperable, and Reproducible) single-cell ‘omics data analysis, yet preprocessing steps lack workflow infrastructure needed to standardize large-scale and collaborative studies. 10x Genomics’ [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) is critical software for preprocessing raw single-cell ‘omics modalities, but executing it reproducibly across hundreds or thousands of samples remains cumbersome, error-prone, and computationally ineﬃcient. We present `cellranger-snakemake`, a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow wrapper that automates, scales, and standardizes Cell Ranger preprocessing for Gene Expression ([GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest)), Chromatin accessibility ([ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest)), and multiome ([ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest)) data. The workflow supports flexible input specifications, integrated logging, and portable configuration files, making it straightforward
to deploy in high-performance computing or cloud environments. By combining Snakemake's reproducible workflow management with Cell Ranger, `cellranger-snakemake` improves reproducibility, reduces user error, and accelerates downstream single-cell 'omics.

---

## Features

- **Per-capture objects**: AnnData/MuData objects created immediately after Cell Ranger count (one per capture), then aggregated into batch-level objects
- **Traceability**: All objects include `batch_id`, `capture_id`, and globally unique `cell_id` metadata
- **Python-only tools**: Doublet detection (Scrublet, SOLO), cell type annotation (celltypist, scANVI, decoupler), demultiplexing (Vireo, demuxalot)
- **Metadata enrichment**: Analysis results automatically merged into final batch objects
- **Multi-modality**: Supports GEX, ATAC, and ARC (multiome) workflows

---

## Documentation

For full documentation, visit **[cellranger-snakemake.readthedocs.io](https://cellranger-snakemake.readthedocs.io/en/latest/)**.

- [Installation](https://cellranger-snakemake.readthedocs.io/en/latest/installation.html)
- [Quick Start](https://cellranger-snakemake.readthedocs.io/en/latest/quickstart.html)
- [Tutorial](https://cellranger-snakemake.readthedocs.io/en/latest/tutorial.html)

# Developer notes

## Project Guidelines

See [CLAUDE.md](CLAUDE.md) for project architecture, coding conventions, and development guidelines. Key topics include:
- Per-capture object creation patterns
- Traceability metadata requirements
- Modality-specific rules (GEX, ATAC, ARC)
- Import ordering conventions

## Updating the Conda environment

If you add or update any Conda dependency, re-export the environment files:
```bash
# Full environment (exact pins, Linux only)
conda env export | grep -v "^prefix:" > environment.yaml

# Portable environment (no build hashes)
conda env export --no-builds | grep -v "^prefix:" > environment_portable.yaml
```

> **Important:** After exporting, manually remove the `cellranger-snakemake` line from the `pip:` section in both files. The package should be installed separately via `pip install -e .`, not pinned in the environment.