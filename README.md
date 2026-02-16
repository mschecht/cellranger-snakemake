# TLDR

The tool `cellranger-snakemake` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) wrapper for [10x Cell Ranger workflows](https://www.10xgenomics.com/software) supporting [GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest), [ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), and [ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) single-cell data processing.

# Description

Reproducibility and scalability are essential components of contemporary [FAIR](https://www.nature.com/articles/sdata201618) (Findable, Accessible, Interoperable, and Reproducible) single-cell ‘omics data analysis, yet preprocessing steps lack workflow infrastructure needed to standardize large-scale and collaborative studies. 10x Genomics’ [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) is critical software for preprocessing raw single-cell ‘omics modalities, but executing it reproducibly across hundreds or thousands of samples remains cumbersome, error-prone, and computationally ineﬃcient. We present `cellranger-snakemake`, a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow wrapper that automates, scales, and standardizes Cell Ranger preprocessing for Gene Expression ([GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest)), Chromatin accessibility ([ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest)), and multiome ([ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest)) data. The workflow supports flexible input specifications, integrated logging, and portable configuration files, making it straightforward
to deploy in high-performance computing or cloud environments. By combining Snakemake's reproducible workflow management with Cell Ranger, `cellranger-snakemake` improves reproducibility, reduces user error, and accelerates downstream single-cell 'omics.

---

## 🚨 What's New in v3.0.0 (2026-02-16)

**Major refactoring with breaking changes!** This release introduces a Python-only, per-capture object architecture:

- **Per-capture objects**: AnnData/MuData objects now created immediately after Cell Ranger count (one per capture) instead of after aggregation
- **Batch aggregation**: New phase merges per-capture objects into batch-level objects with full traceability metadata
- **Python-only tools**: All R-based analysis tools (DoubletFinder, Azimuth, SingleR, ScType) replaced with Python equivalents (SOLO, scANVI, celltypist, decoupler)
- **New directory structure**: Per-capture objects in `03_ANNDATA/`, batch objects in `04_BATCH_OBJECTS/`, downstream steps renumbered
- **Required metadata**: All objects include `batch_id`, `capture_id`, and globally unique `cell_id` for traceability

**⚠️ v2.x configs and objects are incompatible with v3.0.0.** See [CHANGELOG.md](CHANGELOG.md) for full migration notes and feature details.

---

## Documentation

For full documentation, visit **[cellranger-snakemake.readthedocs.io](https://cellranger-snakemake.readthedocs.io/en/latest/)**.

- [Installation](https://cellranger-snakemake.readthedocs.io/en/latest/installation.html)
- [Quick Start](https://cellranger-snakemake.readthedocs.io/en/latest/quickstart.html)
- [Tutorial](https://cellranger-snakemake.readthedocs.io/en/latest/tutorial.html)

# Quick Start

1. **Activate the environment**
```bash
conda activate snakemake8
```

2. **Verify Cell Ranger installation (optional but recommended)**
```bash
snakemake-run-cellranger check-versions
```

3. **Create a config file:**
```bash
snakemake-run-cellranger init-config --output pipeline_config.yaml
```

4. **Run a workflow:**
```bash
snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 8
```

# Tutorial

[Check out the wiki for an in-depth tutorial with publicly available practice data!
](https://github.com/mschecht/cellranger-snakemake/wiki)

# Developer notes

## Project Guidelines

See [CLAUDE.md](CLAUDE.md) for project architecture, coding conventions, and development guidelines. Key topics include:
- Per-capture object creation patterns
- Traceability metadata requirements
- Modality-specific rules (GEX, ATAC, ARC)
- Import ordering conventions

## Updating the Conda environment

If you add or update any Conda dependency, re-export the environment and commit the change:
```bash
conda env export > environment.yaml
```