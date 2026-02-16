"""Run scANVI (scvi-tools) reference-based cell type annotation on per-capture AnnData."""

import scanpy as sc
import scvi
import pandas as pd
import sys
import numpy as np

# Access snakemake object
h5ad_path = snakemake.input.h5ad
reference_path = snakemake.input.reference
output_tsv = snakemake.output.tsv
log_file = snakemake.log[0]

# scANVI parameters
label_key = snakemake.params.label_key
n_hidden = snakemake.params.n_hidden
n_latent = snakemake.params.n_latent
n_layers = snakemake.params.n_layers
max_epochs = snakemake.params.max_epochs

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading query AnnData from: {h5ad_path}")
    adata_query = sc.read_h5ad(h5ad_path)
    print(f"Query: {adata_query.n_obs} cells, {adata_query.n_vars} genes")

    print(f"\nReading reference AnnData from: {reference_path}")
    adata_ref = sc.read_h5ad(reference_path)
    print(f"Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")

    # Verify label_key exists in reference
    if label_key not in adata_ref.obs.columns:
        raise ValueError(f"Label key '{label_key}' not found in reference.obs. Available columns: {adata_ref.obs.columns.tolist()}")

    print(f"\nReference cell type distribution:")
    print(adata_ref.obs[label_key].value_counts())

    # Find common genes between reference and query
    print("\nFinding common genes...")
    common_genes = np.intersect1d(adata_ref.var_names, adata_query.var_names)
    print(f"Common genes: {len(common_genes)}")

    if len(common_genes) < 500:
        raise ValueError(f"Too few common genes ({len(common_genes)}) between reference and query. Ensure compatible data sources.")

    # Subset to common genes
    adata_ref = adata_ref[:, common_genes].copy()
    adata_query = adata_query[:, common_genes].copy()

    # Preprocessing
    print("\nPreprocessing data...")

    # Basic QC
    sc.pp.filter_genes(adata_ref, min_cells=3)
    sc.pp.filter_genes(adata_query, min_cells=3)

    # Normalize and log-transform
    for adata in [adata_ref, adata_query]:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Find highly variable genes on reference
    sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000, subset=True)

    # Subset query to same genes
    adata_query = adata_query[:, adata_ref.var_names].copy()

    print(f"After preprocessing: {adata_ref.n_obs} ref cells, {adata_query.n_obs} query cells, {adata_ref.n_vars} genes")

    # Concatenate reference and query for joint processing
    print("\nConcatenating reference and query...")
    adata_ref.obs['_batch'] = 'reference'
    adata_query.obs['_batch'] = 'query'
    adata_query.obs[label_key] = 'Unknown'  # Placeholder for query

    adata_full = sc.concat([adata_ref, adata_query], join='inner')
    print(f"Combined: {adata_full.n_obs} cells")

    # Setup AnnData for scvi
    print("\nSetting up AnnData for scvi-tools...")
    scvi.model.SCVI.setup_anndata(
        adata_full,
        batch_key='_batch',
        labels_key=label_key
    )

    # Train scVI model
    print(f"\nTraining scVI model with parameters:")
    print(f"  n_hidden: {n_hidden}")
    print(f"  n_latent: {n_latent}")
    print(f"  n_layers: {n_layers}")

    vae = scvi.model.SCVI(
        adata_full,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers
    )

    print(f"Training scVI for {max_epochs} epochs...")
    vae.train(max_epochs=max_epochs, train_size=0.9, early_stopping=True)

    # Train scANVI model
    print("\nTraining scANVI model for label transfer...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        vae,
        unlabeled_category="Unknown",
        labels_key=label_key
    )

    print(f"Training scANVI for {max_epochs} epochs...")
    scanvi_model.train(max_epochs=max_epochs, train_size=0.9)

    # Get predictions for query cells
    print("\nPredicting cell types for query cells...")
    predictions = scanvi_model.predict()

    # Filter to only query cells
    query_mask = adata_full.obs['_batch'] == 'query'
    query_predictions = predictions[query_mask]

    # Get prediction probabilities (confidence)
    pred_probs = scanvi_model.predict(soft=True)
    query_probs = pred_probs[query_mask]
    confidence_scores = query_probs.max(axis=1)

    print(f"\nAnnotation summary:")
    unique_labels, counts = pd.Series(query_predictions).value_counts().index, pd.Series(query_predictions).value_counts().values
    for label, count in zip(unique_labels[:10], counts[:10]):
        print(f"  {label}: {count} cells")

    # Create output dataframe
    output_df = pd.DataFrame({
        'barcode': adata_query.obs.index,
        'scanvi_predicted_label': query_predictions,
        'scanvi_confidence': confidence_scores
    })

    # Write to compressed TSV
    print(f"\nWriting results to: {output_tsv}")
    output_df.to_csv(output_tsv, sep='\t', index=False, compression='gzip')
    print(f"✓ scANVI annotation complete")

except Exception as e:
    print(f"ERROR running scANVI: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
