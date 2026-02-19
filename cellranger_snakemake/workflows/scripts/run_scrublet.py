"""Run Scrublet doublet detection on per-capture AnnData object."""

import re
import sys

import scanpy as sc
import pandas as pd

# Access snakemake object
h5ad_path = snakemake.input.h5ad
output_tsv = snakemake.output.tsv
log_file = snakemake.log[0]

# Scrublet parameters
expected_doublet_rate = snakemake.params.expected_doublet_rate
min_counts = snakemake.params.min_counts
min_cells = snakemake.params.min_cells
min_gene_variability_pctl = snakemake.params.min_gene_variability_pctl
n_prin_comps = snakemake.params.n_prin_comps

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    print(f"\nRunning Scrublet doublet detection with parameters:")
    print(f"  expected_doublet_rate: {expected_doublet_rate}")
    print(f"  min_counts: {min_counts}")
    print(f"  min_cells: {min_cells}")
    print(f"  min_gene_variability_pctl: {min_gene_variability_pctl}")
    print(f"  n_prin_comps: {n_prin_comps}")

    try:
        sc.pp.scrublet(
            adata,
            expected_doublet_rate=expected_doublet_rate,
            n_prin_comps=n_prin_comps,
        )
    except ValueError as e:
        # Scrublet filters to highly variable genes internally, so the effective
        # feature count can be much smaller than adata.n_vars. Parse the actual
        # limit from the error message and retry.
        match = re.search(r'min\(n_samples, n_features\)=(\d+)', str(e))
        if not match:
            raise
        max_comps = int(match.group(1)) - 1
        print(f"  n_prin_comps reduced to {max_comps} (after internal gene filtering)")
        sc.pp.scrublet(
            adata,
            expected_doublet_rate=expected_doublet_rate,
            n_prin_comps=max_comps,
        )

    n_doublets = adata.obs['predicted_doublet'].sum()
    print(f"Detected {n_doublets} doublets out of {adata.n_obs} cells")
    print(f"Doublet rate: {n_doublets / adata.n_obs:.2%}")

    output_df = pd.DataFrame({
        'barcode': adata.obs.index,
        'scrublet_score': adata.obs['doublet_score'],
        'scrublet_predicted_doublet': adata.obs['predicted_doublet'],
    })

    print(f"\nWriting results to: {output_tsv}")
    output_df.to_csv(output_tsv, sep='\t', index=False, compression='gzip')
    print("✓ Scrublet doublet detection complete")

except Exception as e:
    print(f"ERROR running Scrublet: {e}", file=sys.stderr)
    raise

finally:
    sys.stdout.close()
