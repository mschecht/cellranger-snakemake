"""Create MuData object from ARC Cell Ranger count output with GEX and ATAC modalities."""

import muon as mu
import scanpy as sc
import sys

# Access snakemake object
h5_path = snakemake.input.h5
fragments_path = snakemake.input.fragments
batch_id = snakemake.params.batch
capture_id = snakemake.params.capture
output_h5mu = snakemake.output.h5mu
log_file = snakemake.log[0]

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading ARC multiome data from: {h5_path}")

    # Read the multiome h5 file which contains both GEX and ATAC
    adata_full = sc.read_10x_h5(h5_path)
    print(f"Loaded {adata_full.n_obs} cells and {adata_full.n_vars} features")

    # Separate GEX and ATAC modalities based on feature_types
    gex_mask = adata_full.var['feature_types'] == 'Gene Expression'
    atac_mask = adata_full.var['feature_types'] == 'Peaks'

    adata_gex = adata_full[:, gex_mask].copy()
    adata_atac = adata_full[:, atac_mask].copy()

    print(f"GEX: {adata_gex.n_obs} cells, {adata_gex.n_vars} genes")
    print(f"ATAC: {adata_atac.n_obs} cells, {adata_atac.n_vars} peaks")

    # Store fragment file path in ATAC modality
    adata_atac.uns['fragments_file'] = fragments_path

    # Add traceability metadata to both modalities
    for adata_mod in [adata_gex, adata_atac]:
        adata_mod.obs['batch_id'] = str(batch_id)
        adata_mod.obs['capture_id'] = str(capture_id)

        # Create unique cell ID
        adata_mod.obs['cell_id'] = (
            adata_mod.obs['batch_id'].astype(str) + '_' +
            adata_mod.obs['capture_id'].astype(str) + '_' +
            adata_mod.obs.index.astype(str)
        )

        # Verify uniqueness
        if not adata_mod.obs['cell_id'].is_unique:
            raise ValueError("cell_id is not unique!")

    print(f"\nGEX metadata columns: {adata_gex.obs.columns.tolist()}")
    print(f"Batch: {batch_id}, Capture: {capture_id}")
    print(f"First few GEX cells:\n{adata_gex.obs.head()}")

    # Create MuData object with GEX and ATAC modalities
    mdata = mu.MuData({'gex': adata_gex, 'atac': adata_atac})
    print(f"\nCreated MuData with modalities: {list(mdata.mod.keys())}")

    # Write to disk
    print(f"\nWriting MuData to: {output_h5mu}")
    mdata.write(output_h5mu)
    print(f"✓ ARC MuData creation complete: {mdata.n_obs} cells")

except Exception as e:
    print(f"ERROR creating ARC MuData: {e}", file=sys.stderr)
    raise

finally:
    sys.stdout.close()
