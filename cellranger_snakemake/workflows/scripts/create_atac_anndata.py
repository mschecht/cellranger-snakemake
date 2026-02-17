"""Create AnnData object from ATAC Cell Ranger count output using SnapATAC2."""

import gzip
import subprocess
import sys
import os

import snapatac2 as snap

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


def get_chrom_sizes_from_reference(fragments_path):
    """Extract chromosome sizes from reference path in fragments file header."""
    print("Extracting reference path from fragments file...")

    # Read fragments file header to get reference path
    reference_path = None
    with gzip.open(fragments_path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('# reference_path='):
                reference_path = line.strip().split('=', 1)[1]
                break

    if not reference_path:
        raise ValueError("Could not find reference_path in fragments file header")

    print(f"Reference path: {reference_path}")

    # Read chromosome sizes from reference .fai file
    fai_file = os.path.join(reference_path, "fasta", "genome.fa.fai")

    if not os.path.exists(fai_file):
        raise FileNotFoundError(f"Reference .fai file not found: {fai_file}")

    print(f"Reading chromosome sizes from: {fai_file}")

    chrom_sizes = {}
    with open(fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom_name = fields[0]
            chrom_length = int(fields[1])
            chrom_sizes[chrom_name] = chrom_length

    print(f"✓ Found {len(chrom_sizes)} chromosomes: {list(chrom_sizes.keys())}")

    return chrom_sizes


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

    # Get chromosome sizes from reference
    chrom_sizes = get_chrom_sizes_from_reference(fragments_path)

    # Import fragments using SnapATAC2 (computes QC metrics)
    adata = snap.pp.import_fragments(
        sorted_fragments,
        chrom_sizes=chrom_sizes,
        sorted_by_barcode=True,  # Use sorted mode for efficiency
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
