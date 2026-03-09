"""Create AnnData object from GEX Cell Ranger count output with traceability metadata."""

import scanpy as sc
import sys

# Access snakemake object
h5_path = snakemake.input.h5
batch_id = snakemake.params.batch
capture_id = snakemake.params.capture
output_h5ad = snakemake.output.h5ad
log_file = snakemake.log[0]

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading h5 file from: {h5_path}")
    adata = sc.read_10x_h5(h5_path)
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    # Add traceability metadata
    adata.obs['batch_id'] = str(batch_id)
    adata.obs['capture_id'] = str(capture_id)

    # Create unique cell ID: batch_capture_barcode
    adata.obs['cell_id'] = (
        adata.obs['batch_id'].astype(str) + '_' +
        adata.obs['capture_id'].astype(str) + '_' +
        adata.obs.index.astype(str)
    )

    # Verify cell_id uniqueness
    if not adata.obs['cell_id'].is_unique:
        raise ValueError("cell_id is not unique! This should never happen.")

    print(f"\nCell metadata columns: {adata.obs.columns.tolist()}")
    print(f"Batch: {batch_id}, Capture: {capture_id}")
    print(f"First few cells:\n{adata.obs.head()}")

    # Store raw counts
    adata.raw = adata

    # Compute QC metrics
    qc_vars = []
    if adata.var_names.str.startswith(('MT-', 'MT.')).any():
        adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'MT.'))
        qc_vars.append('mt')
    if adata.var_names.str.startswith(('RPS', 'RPL')).any():
        adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
        qc_vars.append('ribo')
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, inplace=True)

    # Write to disk
    print(f"\nWriting AnnData to: {output_h5ad}")
    adata.write_h5ad(output_h5ad)
    print(f"✓ GEX AnnData creation complete: {adata.n_obs} cells")

except Exception as e:
    print(f"ERROR creating GEX AnnData: {e}", file=sys.stderr)
    raise

finally:
    sys.stdout.close()
