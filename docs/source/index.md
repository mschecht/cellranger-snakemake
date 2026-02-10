# single-cell-snakemake documentation

## TLDR

The tool `cellranger-snakemake` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) wrapper for [10x Cell Ranger workflows](https://www.10xgenomics.com/software) supporting [GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest), [ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest), and [ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest) single-cell data processing.

## Description

Reproducibility and scalability are essential components of contemporary [FAIR](https://www.nature.com/articles/sdata201618) (Findable, Accessible, Interoperable, and Reproducible) single-cell 'omics data analysis, yet preprocessing steps lack workflow infrastructure needed to standardize large-scale and collaborative studies. 10x Genomics' [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) is critical software for preprocessing raw single-cell 'omics modalities, but executing it reproducibly across hundreds or thousands of samples remains cumbersome, error-prone, and computationally inefficient. We present `cellranger-snakemake`, a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow wrapper that automates, scales, and standardizes Cell Ranger preprocessing for Gene Expression ([GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest)), Chromatin accessibility ([ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest)), and multiome ([ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest)) data. The workflow supports flexible input specifications, integrated logging, and portable configuration files, making it straightforward to deploy in high-performance computing or cloud environments. By combining Snakemake's reproducible workflow management with Cell Ranger, `cellranger-snakemake` improves reproducibility, reduces user error, and accelerates downstream single-cell 'omics.

```{toctree}
:maxdepth: 2
:caption: Contents:

installation
quickstart
tutorial
development
```
