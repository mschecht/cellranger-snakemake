"""Merge analysis metadata into batch-level AnnData/MuData objects."""

import os
import sys
import glob
import pandas as pd
import scanpy as sc

# Access snakemake object
batch_id = snakemake.params.batch
modality = snakemake.params.modality
input_object = snakemake.input.batch_object
output_object = snakemake.output.enriched_object
log_file = snakemake.log[0]

# Directories for metadata files
demux_dir = snakemake.params.get("demux_dir")
doublet_dir = snakemake.params.get("doublet_dir")
annotation_dir = snakemake.params.get("annotation_dir")

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"{'='*70}")
    print(f"Enriching batch {batch_id} {modality.upper()} object with metadata")
    print(f"{'='*70}\n")

    # Read batch object
    print(f"Reading batch object: {input_object}")
    if modality == "arc":
        import muon as mu
        adata = mu.read(input_object)
        is_mudata = True
    else:
        adata = sc.read_h5ad(input_object)
        is_mudata = False

    print(f"✓ Loaded {adata.n_obs} cells\n")

    metadata_added = []

    # ========================================================================
    # Merge Demultiplexing Metadata
    # ========================================================================
    if demux_dir and os.path.exists(demux_dir):
        print("Searching for demultiplexing metadata...")

        # Find all per-capture demux results
        # Pattern: {demux_dir}/{batch}_{capture}/{method}/{batch}_{capture}_assignments.tsv.gz
        demux_files = glob.glob(os.path.join(demux_dir, f"{batch_id}_*", "*", f"{batch_id}_*_assignments.tsv.gz"))

        if demux_files:
            print(f"Found {len(demux_files)} demux file(s)")
            all_demux = []

            for demux_file in demux_files:
                print(f"  - {os.path.basename(demux_file)}")

                # Extract capture ID from filename: {batch}_{capture}_assignments.tsv.gz
                basename = os.path.basename(demux_file)
                # Remove prefix and suffix more carefully to avoid multiple replacements
                if basename.startswith(f"{batch_id}_"):
                    capture_id = basename[len(f"{batch_id}_"):]  # Remove prefix
                    capture_id = capture_id.replace("_assignments.tsv.gz", "")  # Remove suffix
                else:
                    raise ValueError(f"Unexpected demux filename format: {basename}")

                # Read demux results
                df = pd.read_csv(demux_file, sep="\t", index_col=0)

                # Create cell_id column for matching: {batch}_{capture}_{barcode}
                df['cell_id'] = df.index.map(lambda barcode: f"{batch_id}_{capture_id}_{barcode}")

                # Add method name from directory
                method = os.path.basename(os.path.dirname(demux_file))
                df.columns = [f"demux_{method}_{col}" if col != 'cell_id' else col for col in df.columns]

                all_demux.append(df)

            # Merge all demux results (concatenate vertically - different captures)
            demux_metadata = pd.concat(all_demux, axis=0)

            # Set cell_id as index for joining
            demux_metadata = demux_metadata.set_index('cell_id')

            # Join with object on cell_id
            if is_mudata:
                # For MuData, add to each modality
                for mod_name, mod_data in adata.mod.items():
                    temp_df = mod_data.obs.set_index('cell_id')
                    temp_df = temp_df.join(demux_metadata, how='left')
                    temp_df = temp_df.reset_index(drop=False)
                    mod_data.obs = temp_df
                    matched = temp_df[demux_metadata.columns].notna().any(axis=1).sum()
                    print(f"✓ Added demux metadata to {mod_name} ({matched} cells matched)")
            else:
                # Create temporary df with cell_id as index
                temp_df = adata.obs.set_index('cell_id')
                temp_df = temp_df.join(demux_metadata, how='left')
                temp_df = temp_df.reset_index(drop=False)
                adata.obs = temp_df
                matched = temp_df[demux_metadata.columns].notna().any(axis=1).sum()
                print(f"✓ Added demux metadata ({matched} cells matched)")

            metadata_added.append("demultiplexing")
        else:
            print("  No demux files found")

    print()

    # ========================================================================
    # Merge Doublet Detection Metadata
    # ========================================================================
    if doublet_dir and os.path.exists(doublet_dir):
        print("Searching for doublet detection metadata...")

        # Find all doublet results for this batch
        # Pattern: {doublet_dir}/{batch}/{method}/{batch}_scores.tsv.gz or similar
        doublet_files = glob.glob(os.path.join(doublet_dir, f"{batch_id}", "*", "*.tsv.gz"))

        if doublet_files:
            print(f"Found {len(doublet_files)} doublet file(s)")
            all_doublet = []

            for doublet_file in doublet_files:
                print(f"  - {os.path.basename(doublet_file)}")
                df = pd.read_csv(doublet_file, sep="\t", index_col=0)

                # Add method name from directory
                method = os.path.basename(os.path.dirname(doublet_file))
                df.columns = [f"doublet_{method}_{col}" for col in df.columns]
                all_doublet.append(df)

            # Merge all doublet results (concatenate vertically if per-capture, horizontally if per-method)
            doublet_metadata = pd.concat(all_doublet, axis=0)

            # Join with object
            if is_mudata:
                for mod_name, mod_data in adata.mod.items():
                    common_barcodes = mod_data.obs.index.intersection(doublet_metadata.index)
                    if len(common_barcodes) > 0:
                        mod_data.obs = mod_data.obs.join(doublet_metadata, how='left')
                        print(f"✓ Added doublet metadata to {mod_name} ({len(common_barcodes)} cells matched)")
            else:
                common_barcodes = adata.obs.index.intersection(doublet_metadata.index)
                if len(common_barcodes) > 0:
                    adata.obs = adata.obs.join(doublet_metadata, how='left')
                    print(f"✓ Added doublet metadata ({len(common_barcodes)} cells matched)")

            metadata_added.append("doublet_detection")
        else:
            print("  No doublet files found")

    print()

    # ========================================================================
    # Merge Cell Type Annotation Metadata
    # ========================================================================
    if annotation_dir and os.path.exists(annotation_dir):
        print("Searching for cell type annotation metadata...")

        # Find all annotation results for this batch
        annotation_files = glob.glob(os.path.join(annotation_dir, f"{batch_id}", "*", "*.tsv.gz"))

        if annotation_files:
            print(f"Found {len(annotation_files)} annotation file(s)")
            all_annotation = []

            for annotation_file in annotation_files:
                print(f"  - {os.path.basename(annotation_file)}")
                df = pd.read_csv(annotation_file, sep="\t", index_col=0)

                # Add method name from directory
                method = os.path.basename(os.path.dirname(annotation_file))
                df.columns = [f"celltype_{method}_{col}" for col in df.columns]
                all_annotation.append(df)

            # Merge all annotation results (concatenate vertically if per-capture, horizontally if per-method)
            annotation_metadata = pd.concat(all_annotation, axis=0)

            # Join with object
            if is_mudata:
                for mod_name, mod_data in adata.mod.items():
                    common_barcodes = mod_data.obs.index.intersection(annotation_metadata.index)
                    if len(common_barcodes) > 0:
                        mod_data.obs = mod_data.obs.join(annotation_metadata, how='left')
                        print(f"✓ Added annotation metadata to {mod_name} ({len(common_barcodes)} cells matched)")
            else:
                common_barcodes = adata.obs.index.intersection(annotation_metadata.index)
                if len(common_barcodes) > 0:
                    adata.obs = adata.obs.join(annotation_metadata, how='left')
                    print(f"✓ Added annotation metadata ({len(common_barcodes)} cells matched)")

            metadata_added.append("celltype_annotation")
        else:
            print("  No annotation files found")

    print()

    # ========================================================================
    # Summary
    # ========================================================================
    print(f"{'='*70}")
    print(f"Metadata enrichment summary:")
    print(f"{'='*70}")
    print(f"Batch: {batch_id}")
    print(f"Modality: {modality.upper()}")
    print(f"Cells: {adata.n_obs}")

    if is_mudata:
        print(f"\nMetadata columns per modality:")
        for mod_name, mod_data in adata.mod.items():
            print(f"  {mod_name}: {mod_data.obs.columns.tolist()}")
    else:
        print(f"\nMetadata columns: {adata.obs.columns.tolist()}")

    print(f"\nMetadata added: {', '.join(metadata_added) if metadata_added else 'None (no analysis completed yet)'}")
    print(f"{'='*70}\n")

    # Write enriched object
    print(f"Writing enriched object to: {output_object}")
    if is_mudata:
        adata.write(output_object)
    else:
        adata.write_h5ad(output_object)

    print("✓ Metadata enrichment complete!\n")

except Exception as e:
    print(f"ERROR in metadata enrichment: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    raise

finally:
    sys.stdout.close()
