"""Run Celltypist cell type annotation on per-capture AnnData object."""

import scanpy as sc
import celltypist
import pandas as pd
import sys

# Access snakemake object
h5ad_path = snakemake.input.h5ad
output_tsv = snakemake.output.tsv
log_file = snakemake.log[0]

# Celltypist parameters
model = snakemake.params.model
majority_voting = snakemake.params.majority_voting
over_clustering = snakemake.params.over_clustering
min_prop = snakemake.params.min_prop

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    # Preprocessing for celltypist (minimal - needs log-normalized data)
    print("\nPreprocessing data for celltypist...")

    # Check if data is already normalized
    if 'log1p' not in adata.uns:
        print("Normalizing and log-transforming data...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("Data appears to be log-normalized already")

    # Load model
    print(f"\nLoading celltypist model: {model}")
    try:
        # Try loading as pre-trained model name
        model_obj = celltypist.models.Model.load(model=model)
    except:
        # Otherwise load from file path
        model_obj = celltypist.models.Model.load(model)

    # Run celltypist
    print(f"\nRunning celltypist annotation...")
    print(f"  majority_voting: {majority_voting}")
    if majority_voting:
        print(f"  over_clustering: {over_clustering}")
        print(f"  min_prop: {min_prop}")

    # Annotate
    predictions = celltypist.annotate(
        adata,
        model=model_obj,
        majority_voting=majority_voting,
        over_clustering=over_clustering if majority_voting else None,
        min_prop=min_prop if majority_voting else 0.0
    )

    # Extract predictions
    if majority_voting and over_clustering:
        predicted_labels = predictions.predicted_labels['majority_voting'].values
    else:
        predicted_labels = predictions.predicted_labels['predicted_labels'].values

    # Get confidence scores (probability of predicted label)
    conf_matrix = predictions.probability_matrix
    confidence_scores = conf_matrix.max(axis=1).values

    print(f"\nAnnotation summary:")
    unique_labels, counts = pd.Series(predicted_labels).value_counts().index, pd.Series(predicted_labels).value_counts().values
    for label, count in zip(unique_labels[:10], counts[:10]):  # Show top 10
        print(f"  {label}: {count} cells")
    if len(unique_labels) > 10:
        print(f"  ... and {len(unique_labels) - 10} more cell types")

    # Create output dataframe
    output_df = pd.DataFrame({
        'barcode': adata.obs.index,
        'celltypist_predicted_label': predicted_labels,
        'celltypist_confidence': confidence_scores
    })

    # Write to compressed TSV
    print(f"\nWriting results to: {output_tsv}")
    output_df.to_csv(output_tsv, sep='\t', index=False, compression='gzip')
    print(f"✓ Celltypist annotation complete")

except Exception as e:
    print(f"ERROR running celltypist: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
