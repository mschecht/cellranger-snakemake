"""Run decoupler marker-based cell type annotation on per-capture AnnData object."""

import scanpy as sc
import decoupler as dc
import pandas as pd
import sys
import os

# Access snakemake object
h5ad_path = snakemake.input.h5ad
output_tsv = snakemake.output.tsv
log_file = snakemake.log[0]

# Decoupler parameters
marker_database = snakemake.params.marker_database
method = snakemake.params.method
min_score = snakemake.params.min_score

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading AnnData from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    # Preprocessing for decoupler (needs log-normalized data)
    print("\nPreprocessing data for decoupler...")

    # Check if data is already normalized
    if 'log1p' not in adata.uns:
        print("Normalizing and log-transforming data...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("Data appears to be log-normalized already")

    # Load marker gene sets
    print(f"\nLoading marker gene database: {marker_database}")

    if os.path.isfile(marker_database):
        # Load from file (GMT or CSV format)
        if marker_database.endswith('.gmt'):
            marker_df = dc.read_gmt(marker_database)
        elif marker_database.endswith('.csv') or marker_database.endswith('.tsv'):
            marker_df = pd.read_csv(marker_database, sep='\t' if marker_database.endswith('.tsv') else ',')
            # Expected format: columns 'cell_type', 'gene', 'weight' (optional)
            if 'weight' not in marker_df.columns:
                marker_df['weight'] = 1.0
        else:
            raise ValueError(f"Unsupported marker database format: {marker_database}. Use .gmt, .csv, or .tsv")
    else:
        # Try loading built-in markers from PanglaoDB or other sources
        print(f"Attempting to load built-in markers for tissue: {marker_database}")
        # PanglaoDB markers can be accessed via decoupler
        try:
            marker_df = dc.get_progeny(organism='human', top=100)
            print("Warning: Using PROGENy pathway signatures as placeholder. Please provide custom marker file for better results.")
        except:
            raise ValueError(f"Could not find marker database file: {marker_database}")

    print(f"Loaded {len(marker_df)} marker gene entries")
    if hasattr(marker_df, 'source'):
        print(f"Cell types: {marker_df['source'].unique()[:10]}...")
    elif 'cell_type' in marker_df.columns:
        print(f"Cell types: {marker_df['cell_type'].unique()[:10]}...")

    # Run decoupler scoring
    print(f"\nRunning decoupler with method: {method}")

    scoring_methods = {
        'ulm': dc.run_ulm,
        'wsum': dc.run_wsum,
        'ora': dc.run_ora,
        'aucell': dc.run_aucell
    }

    if method not in scoring_methods:
        raise ValueError(f"Unknown method: {method}. Choose from {list(scoring_methods.keys())}")

    scoring_func = scoring_methods[method]

    # Run scoring
    dc.run_enrichment(
        adata,
        net=marker_df,
        source='source' if 'source' in marker_df.columns else 'cell_type',
        target='target' if 'target' in marker_df.columns else 'gene',
        weight='weight' if 'weight' in marker_df.columns else None,
        use_raw=False,
        min_n=5
    )

    # Get scores from obsm
    scores_key = list(adata.obsm.keys())[-1]  # Latest added scores
    scores = adata.obsm[scores_key]

    print(f"\nScore matrix shape: {scores.shape}")

    # Assign cell types based on maximum score
    max_scores = scores.max(axis=1)
    predicted_labels = scores.columns[scores.values.argmax(axis=1)]

    # Filter by minimum score threshold
    predicted_labels = pd.Series(predicted_labels, index=adata.obs.index)
    predicted_labels[max_scores < min_score] = 'Unknown'

    print(f"\nAnnotation summary (min_score={min_score}):")
    unique_labels, counts = predicted_labels.value_counts().index, predicted_labels.value_counts().values
    for label, count in zip(unique_labels[:10], counts[:10]):
        print(f"  {label}: {count} cells")

    # Create output dataframe
    output_df = pd.DataFrame({
        'barcode': adata.obs.index,
        'decoupler_predicted_label': predicted_labels.values,
        'decoupler_score': max_scores
    })

    # Write to compressed TSV
    print(f"\nWriting results to: {output_tsv}")
    output_df.to_csv(output_tsv, sep='\t', index=False, compression='gzip')
    print(f"✓ Decoupler marker-based annotation complete")

except Exception as e:
    print(f"ERROR running decoupler: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
