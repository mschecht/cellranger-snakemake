"""Create AnnData object from ATAC Cell Ranger count output using SnapATAC2."""

import snapatac2 as snap
import subprocess
import sys
import os
import tempfile

# Access snakemake object
fragments_path = snakemake.input.fragments
peak_matrix_path = snakemake.input.peak_matrix
batch_id = snakemake.params.batch
capture_id = snakemake.params.capture
output_h5ad = snakemake.output.h5ad
log_file = snakemake.log[0]

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Importing ATAC fragments from: {fragments_path}")

    # Sort fragments by barcode for efficient import
    # SnapATAC2 is much faster with sorted fragments
    sorted_fragments = fragments_path.replace('.tsv.gz', '.sorted_by_barcode.tsv.gz')

    if not os.path.exists(sorted_fragments):
        print(f"Sorting fragments by barcode for efficient import...")

        # Create temp directory for sorting
        temp_dir = os.path.dirname(fragments_path)

        # Decompress, sort by barcode (col 4), then chromosome and position, recompress
        sort_cmd = f"""
        zcat {fragments_path} | \
        awk 'NR==1 {{print; next}} {{print | "sort -k4,4V -k1,1V -k2,2n -T {temp_dir}"}}' | \
        bgzip > {sorted_fragments}
        """

        print(f"Running sort command...")
        subprocess.run(sort_cmd, shell=True, check=True, executable='/bin/bash')
        print(f"✓ Fragments sorted and saved to: {sorted_fragments}")
    else:
        print(f"Using existing sorted fragments: {sorted_fragments}")

    # Import fragments using SnapATAC2 (computes QC metrics)
    # Using memory mode to allow metadata additions
    adata = snap.pp.import_fragments(
        sorted_fragments,
        chrom_sizes=snap.genome.hg38,  # Default to hg38, can be parameterized later
        sorted_by_barcode=True,  # Now we can use sorted mode for efficiency
        min_num_fragments=1  # Relaxed for testing; default is 200
    )
    print(f"Loaded {adata.n_obs} cells with fragment-level data")

    # Add traceability metadata
    adata.obs['batch_id'] = str(batch_id)
    adata.obs['capture_id'] = str(capture_id)

    # Create unique cell ID
    adata.obs['cell_id'] = (
        adata.obs['batch_id'].astype(str) + '_' +
        adata.obs['capture_id'].astype(str) + '_' +
        adata.obs.index.astype(str)
    )

    # Verify uniqueness
    if not adata.obs['cell_id'].is_unique:
        raise ValueError("cell_id is not unique!")

    print(f"\nATAC QC metrics: {[col for col in adata.obs.columns if col.startswith(('n_', 'frac_'))]}")
    print(f"Cell metadata columns: {adata.obs.columns.tolist()}")
    print(f"Batch: {batch_id}, Capture: {capture_id}")
    print(f"First few cells:\n{adata.obs.head()}")

    # Write to disk
    print(f"\nWriting ATAC AnnData to: {output_h5ad}")
    adata.write(output_h5ad)
    print(f"✓ ATAC AnnData creation complete: {adata.n_obs} cells")

except Exception as e:
    print(f"ERROR creating ATAC AnnData: {e}", file=sys.stderr)
    raise

finally:
    sys.stdout.close()
