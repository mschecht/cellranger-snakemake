"""Run Scrublet doublet detection on per-capture AnnData object."""

import scanpy as sc
import scrublet as scr
import pandas as pd
import sys

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

    # Initialize Scrublet
    print(f"\nInitializing Scrublet with parameters:")
    print(f"  expected_doublet_rate: {expected_doublet_rate}")
    print(f"  min_counts: {min_counts}")
    print(f"  min_cells: {min_cells}")
    print(f"  min_gene_variability_pctl: {min_gene_variability_pctl}")
    print(f"  n_prin_comps: {n_prin_comps}")

    scrub = scr.Scrublet(
        adata.X,
        expected_doublet_rate=expected_doublet_rate
    )

    # Run doublet detection
    print("\nRunning Scrublet doublet detection...")
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=min_counts,
        min_cells=min_cells,
        min_gene_variability_pctl=min_gene_variability_pctl,
        n_prin_comps=n_prin_comps
    )

    print(f"Detected {predicted_doublets.sum()} doublets out of {len(predicted_doublets)} cells")
    print(f"Doublet rate: {predicted_doublets.sum() / len(predicted_doublets):.2%}")

    # Create output dataframe with barcode index
    output_df = pd.DataFrame({
        'barcode': adata.obs.index,
        'scrublet_score': doublet_scores,
        'scrublet_predicted_doublet': predicted_doublets
    })

    # Write to compressed TSV
    print(f"\nWriting results to: {output_tsv}")
    output_df.to_csv(output_tsv, sep='\t', index=False, compression='gzip')
    print(f"✓ Scrublet doublet detection complete")

except Exception as e:
    print(f"ERROR running Scrublet: {e}", file=sys.stderr)
    raise

finally:
    sys.stdout.close()
