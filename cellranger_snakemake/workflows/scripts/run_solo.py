"""Run SOLO (scvi-tools) doublet detection on per-capture AnnData object."""

import scanpy as sc
import scvi
import pandas as pd
import sys

# Access snakemake object
h5ad_path = snakemake.input.h5ad
output_tsv = snakemake.output.tsv
log_file = snakemake.log[0]

# SOLO parameters
n_hidden = snakemake.params.n_hidden
n_latent = snakemake.params.n_latent
n_layers = snakemake.params.n_layers
learning_rate = snakemake.params.learning_rate
max_epochs = snakemake.params.max_epochs

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    # Preprocessing for SOLO
    print("\nPreprocessing data for SOLO...")

    # Basic QC filtering (minimal - assume cellranger did most)
    sc.pp.filter_genes(adata, min_cells=3)

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    print(f"After preprocessing: {adata.n_obs} cells, {adata.n_vars} genes")

    # Setup AnnData for scvi
    print("\nSetting up AnnData for scvi-tools...")
    scvi.model.SCVI.setup_anndata(adata)

    # Train scVI model
    print(f"\nTraining scVI model with parameters:")
    print(f"  n_hidden: {n_hidden}")
    print(f"  n_latent: {n_latent}")
    print(f"  n_layers: {n_layers}")

    vae = scvi.model.SCVI(
        adata,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers
    )

    print(f"Training for {max_epochs} epochs with learning_rate={learning_rate}...")
    vae.train(max_epochs=max_epochs, train_size=0.9, early_stopping=True)

    # Train SOLO model
    print("\nTraining SOLO doublet detection model...")
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train(max_epochs=max_epochs, train_size=0.9)

    # Get predictions
    print("\nGenerating doublet predictions...")
    predictions = solo.predict()

    # Extract doublet scores and classifications
    # SOLO outputs columns: 'doublet', 'singlet', 'prediction'
    doublet_scores = predictions['doublet'].values
    predicted_doublets = (predictions['prediction'] == 'doublet').values

    print(f"Detected {predicted_doublets.sum()} doublets out of {len(predicted_doublets)} cells")
    print(f"Doublet rate: {predicted_doublets.sum() / len(predicted_doublets):.2%}")

    # Create output dataframe with original barcode index
    # Note: adata index may have changed due to filtering, so we need to map back
    output_df = pd.DataFrame({
        'barcode': adata.obs.index,
        'solo_score': doublet_scores,
        'solo_predicted_doublet': predicted_doublets
    })

    # Write to compressed TSV
    print(f"\nWriting results to: {output_tsv}")
    output_df.to_csv(output_tsv, sep='\t', index=False, compression='gzip')
    print(f"✓ SOLO doublet detection complete")

except Exception as e:
    print(f"ERROR running SOLO: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
