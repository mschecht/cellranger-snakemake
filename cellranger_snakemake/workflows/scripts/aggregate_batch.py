"""Aggregate per-capture objects into batch-level objects."""

import scanpy as sc
import sys
import os

# Access snakemake object
batch_id = snakemake.params.batch
modality = snakemake.params.modality
output_file = snakemake.output.h5ad if modality in ["gex", "atac"] else snakemake.output.h5mu
log_file = snakemake.log[0]

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Aggregating batch {batch_id} for {modality.upper()} modality")

    if modality == "arc":
        import muon as mu
        input_files = snakemake.input.h5mus
        print(f"Found {len(input_files)} per-capture MuData files to aggregate:")
    else:
        input_files = snakemake.input.h5ads
        print(f"Found {len(input_files)} per-capture AnnData files to aggregate:")

    for f in input_files:
        print(f"  - {f}")

    # Read all objects
    print("\nReading per-capture objects...")
    objects = []
    for i, file_path in enumerate(input_files):
        print(f"  [{i+1}/{len(input_files)}] Reading: {os.path.basename(file_path)}")

        if modality == "arc":
            obj = mu.read(file_path)
        else:
            obj = sc.read_h5ad(file_path)

        objects.append(obj)
        print(f"      Loaded {obj.n_obs} cells")

    # Concatenate objects
    print(f"\nConcatenating {len(objects)} objects...")

    if modality == "arc":
        # For MuData, use muon.concat
        batch_obj = mu.concat(objects, axis=0, join='outer', merge='same')
    else:
        # For AnnData (GEX/ATAC), use anndata.concat
        import anndata as ad
        batch_obj = ad.concat(objects, axis=0, join='outer', merge='same')

    print(f"✓ Concatenated to {batch_obj.n_obs} cells total")

    # Verify cell_id uniqueness
    if modality == "arc":
        # For MuData, check each modality
        for mod_name, mod_data in batch_obj.mod.items():
            if not mod_data.obs['cell_id'].is_unique:
                raise ValueError(f"cell_id is not unique in {mod_name} modality after aggregation!")
        print(f"✓ cell_id is unique in all modalities")
    else:
        # For AnnData
        if not batch_obj.obs['cell_id'].is_unique:
            raise ValueError("cell_id is not unique after aggregation!")
        print(f"✓ cell_id is unique")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Batch {batch_id} {modality.upper()} aggregation summary:")
    print(f"{'='*60}")
    print(f"Total cells: {batch_obj.n_obs}")

    if modality == "arc":
        print(f"Modalities: {list(batch_obj.mod.keys())}")
        for mod_name, mod_data in batch_obj.mod.items():
            print(f"  {mod_name}: {mod_data.n_obs} cells x {mod_data.n_vars} features")
    else:
        print(f"Features: {batch_obj.n_vars}")

    print(f"Metadata columns: {batch_obj.obs.columns.tolist()}")

    # Show batch/capture distribution
    print(f"\nCapture distribution:")
    print(batch_obj.obs['capture_id'].value_counts())

    print(f"\n{'='*60}")

    # Write batch object
    print(f"\nWriting batch-level object to: {output_file}")

    if modality == "arc":
        batch_obj.write(output_file)
    else:
        batch_obj.write_h5ad(output_file)

    print(f"✓ Batch aggregation complete!")

except Exception as e:
    print(f"ERROR in batch aggregation: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
