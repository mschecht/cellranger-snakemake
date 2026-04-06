"""Create MuData object from ARC Cell Ranger count output with GEX and ATAC modalities."""

import sys

import muon as mu
import scanpy as sc

# Access snakemake object
mtx_dir = snakemake.params.mtx_dir
fragments_path = snakemake.input.fragments
batch_id = snakemake.params.batch
capture_id = snakemake.params.capture
output_h5mu = snakemake.output.h5mu
log_file = snakemake.log[0]

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

try:
    print(f"Reading ARC multiome data from: {mtx_dir}")

    # Reads both GEX and ATAC, returning a MuData with modalities named "rna" and "atac".
    # extended=False skips auto-loading peak annotation (mixed types cause h5 write errors)
    # and fragment file auto-location (we set fragments_path manually below).
    mdata = mu.read_10x_mtx(mtx_dir, extended=False)

    print(f"Loaded modalities: {list(mdata.mod.keys())}")
    print(f"rna: {mdata['rna'].n_obs} cells, {mdata['rna'].n_vars} genes")
    print(f"atac: {mdata['atac'].n_obs} cells, {mdata['atac'].n_vars} peaks")
    
    # Rename "rna" -> "gex" to match pipeline conventions
    print(f'Renaming "rna" modality to "gex"')
    mdata.mod["gex"] = mdata.mod.pop("rna")

    print(mdata)

    # Store fragment file path in ATAC modality
    print(f"Storing fragment file path in ATAC modality: {fragments_path}")
    mdata["atac"].uns["fragments_file"] = fragments_path

    # Add traceability metadata to both modalities
    for mod_name in ["gex", "atac"]:
        adata = mdata[mod_name]
        adata.obs["batch_id"] = str(batch_id)
        adata.obs["capture_id"] = str(capture_id)
        adata.obs["cell_id"] = (
            str(batch_id) + "_" + str(capture_id) + "_" + adata.obs_names.astype(str)
        )
        if not adata.obs["cell_id"].is_unique:
            raise ValueError(f"cell_id is not unique in {mod_name} modality!")
        adata.obs["barcode"] = adata.obs_names
        adata.obs_names = adata.obs["cell_id"]

    # Sync top-level obs index with the updated modality obs_names
    mdata.update()

    # Promote traceability metadata to MuData top-level obs
    mdata.obs["batch_id"] = mdata["gex"].obs["batch_id"]
    mdata.obs["capture_id"] = mdata["gex"].obs["capture_id"]
    mdata.obs["cell_id"] = mdata["gex"].obs["cell_id"]

    # Compute QC metrics
    for mod_name, adata in [("gex", mdata["gex"]), ("atac", mdata["atac"])]:
        qc_vars = []
        if mod_name == "gex":
            if adata.var_names.str.startswith(("MT-", "MT.")).any():
                adata.var["mt"] = adata.var_names.str.startswith(("MT-", "MT."))
                qc_vars.append("mt")
            if adata.var_names.str.startswith(("RPS", "RPL")).any():
                adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
                qc_vars.append("ribo")
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, inplace=True)

    print(f"\nWriting MuData to: {output_h5mu}")
    mdata.write(output_h5mu)
    print(f"✓ ARC MuData creation complete: {mdata.n_obs} cells")
    print(f"  gex: {mdata['gex'].n_obs} x {mdata['gex'].n_vars}")
    print(f"  atac: {mdata['atac'].n_obs} x {mdata['atac'].n_vars}")

except Exception as e:
    print(f"ERROR creating ARC MuData: {e}", file=sys.stderr)
    raise

finally:
    sys.stdout.close()
