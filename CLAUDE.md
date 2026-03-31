# cellranger-snakemake Project Guidelines

## Project Overview

This is a Snakemake pipeline for processing single-cell RNA-seq, ATAC-seq, and multiome (ARC) data from 10X Genomics Cell Ranger outputs. The pipeline creates AnnData/MuData objects with per-capture traceability and supports demultiplexing, doublet detection, and cell type annotation.

**Version:** 0.1.0-dev (pre-release)
**Architecture:** Per-capture objects with batch aggregation
**Language:** Python

---

## Core Principles

### 1. Simplicity Over Complexity
- Keep solutions minimal and focused
- Don't add features beyond what's requested
- Avoid premature optimization or abstraction
- Three similar lines are better than a premature abstraction

### 2. Python-Only Stack
- **NEVER** use R-based tools (e.g., DoubletFinder, Azimuth, SingleR)
- Always use Python equivalents: scanpy, snapatac2, muon, scvi-tools, celltypist, decoupler
- This is a **hard requirement**

### 3. Traceability is Critical
- **ALWAYS** ensure every cell has: `batch_id`, `capture_id`, `cell_id`
- `cell_id` **MUST** be unique across entire dataset
- Format: `{batch}_{capture}_{barcode}`
- **NEVER** skip uniqueness verification: `assert adata.obs['cell_id'].is_unique`

---

## Architecture Requirements

### Per-Capture Pattern
- Objects are created **per-capture**, not per-batch
- Naming: `{batch}_{capture}.h5ad` or `{batch}_{capture}.h5mu`
- Location: `03_ANNDATA/`
- Each capture gets its own object immediately after Cell Ranger count

### Batch Aggregation
- Batch-level objects are created by merging per-capture objects
- Location: `04_BATCH_OBJECTS/{batch}_{modality}.h5ad`
- Preserve all per-capture metadata during merging
- Always verify `cell_id` uniqueness after aggregation

### Directory Structure (Fixed)
```
00_LOGS/                    - Log files and done flags
01_CELLRANGER{GEX|ATAC|ARC}_COUNT/  - Per-capture Cell Ranger outputs
02_CELLRANGER{GEX|ATAC|ARC}_AGGR/   - Cell Ranger batch aggregation
03_ANNDATA/                 - Per-capture .h5ad/.h5mu objects
04_BATCH_OBJECTS/           - Batch-level aggregated objects
05_DEMULTIPLEXING/          - Demux outputs (per-capture)
06_DOUBLET_DETECTION/       - Doublet outputs (per-capture)
07_CELLTYPE_ANNOTATION/     - Annotation outputs (per-capture)
08_FINAL/                   - Future use
```

**Do NOT** renumber directories without explicit approval.

---

## Modality-Specific Rules

### GEX (Gene Expression)
- **Library:** scanpy
- **Input:** `filtered_feature_bc_matrix.h5`
- **Read function:** `sc.read_10x_h5()`
- **Output:** `.h5ad` (AnnData)

### ATAC (Chromatin Accessibility)
- **Library:** snapatac2
- **Input:** `fragments.tsv.gz`
- **Read function:** `snap.pp.import_fragments()`
- **IMPORTANT:** Always sort fragments by barcode first
  - Creates cached file: `fragments.sorted_by_barcode.tsv.gz`
  - Use `sorted_by_barcode=True` for efficiency
- **QC Parameter:** Adjust `min_num_fragments` (default 200 is too high for test data)
- **Output:** `.h5ad` (AnnData with ATAC-specific QC metrics)

### ARC (Multiome: GEX + ATAC)
- **Library:** muon
- **Input:** `filtered_feature_bc_matrix.h5` (contains both modalities)
- **Process:** Read → split by `feature_types` → create MuData
- **IMPORTANT:** Add traceability to **both** GEX and ATAC modalities
- **Output:** `.h5mu` (MuData)

---

## Code Style & Conventions

### Python
- Follow PEP 8
- Use snake_case for functions and variables
- Type hints preferred but not required
- Clear, descriptive variable names

#### Import Order (by character length)
Imports must be organized in three groups, sorted by character length within each group:

1. **Standard library** imports
2. **Third-party** imports (e.g., scanpy, muon, pandas)
3. **Local** imports (from cellranger_snakemake)

**Example:**
```python
import sys

import muon as mu
import scanpy as sc

from cellranger_snakemake.config_validator import parse_output_directories
```

**Rules:**
- Blank line between groups
- Within each group: sort by total line length (shortest first)
- `import X` comes before `import X as Y` if same base length
- `from X import Y` statements last in their group

### Snakemake Rules
- Use per-capture wildcards: `{batch}_{capture}`
- Include descriptive docstrings for each rule
- Always use `.done` files for checkpoints
- Log files: `{batch}_{capture}_{step}.log`
- **Always include `threads:` and `resources:` blocks** — without them SLURM defaults to 1 GB and jobs will be OOM-killed on real data
  - Add `threads` and `mem_gb` to the step's top-level Pydantic config class
  - Use `from tempfile import gettempdir` (must be explicitly imported in each `.smk` file)
  - `RESOURCES` is a global from `main.smk` available in all rule files

### File Naming
- Per-capture objects: `{batch}_{capture}.h5ad` or `.h5mu`
- Batch objects: `{batch}_{modality}.h5ad` or `.h5mu`
- Log files: `{batch}_{capture}_{step}.log`
- Done flags: `{batch}_{capture}_{step}.done`

---

## Required Tools & Libraries

### Core
- `snakemake>=8.0`
- `python>=3.10`

### Single-Cell Analysis
- `scanpy>=1.10` (GEX) — installed: 1.12
- `anndata>=0.10` — installed: 0.12.10
- `muon>=0.1.6` (ARC) — installed: 0.1.7
- `snapatac2>=2.6` (ATAC) — installed: 2.9.0

### Analysis Tools
- `scvi-tools>=1.1` (SOLO, scANVI)
- `celltypist>=1.6` (cell type annotation)
- `decoupler>=1.7` (marker-based annotation)
- `scrublet>=0.2.3` (doublet detection)

### Other
- `pandas`, `numpy`, `leidenalg`

---

## Testing Guidelines

### Test Data Locations
- GEX: `tests/00_TEST_DATA_GEX/` (2 captures: L001, L002)
- ATAC: `tests/00_TEST_DATA_ATAC/` (1 capture: L001)
- ARC: `tests/00_TEST_DATA_ARC/` (1 capture: L001)

### Generate Test Data
```bash
snakemake-run-cellranger generate-test-data GEX --output-dir tests/00_TEST_DATA_GEX
snakemake-run-cellranger generate-test-data ATAC --output-dir tests/00_TEST_DATA_ATAC
snakemake-run-cellranger generate-test-data ARC --output-dir tests/00_TEST_DATA_ARC
```

### Run Tests

#### GEX
```bash
# Dry run
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dry-run

# DAG visualization
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dag | dot -Tpng > dag_gex.png

# Local run
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1

# Run with HPC profile
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml \
                             --cores all \
                             --snakemake-args --profile tests/00_TEST_DATA_GEX/HPC_profiles
```

#### ATAC
```bash
# Dry run
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ATAC/test_config_atac.yaml --cores 1 --dry-run

# DAG visualization
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ATAC/test_config_atac.yaml --cores 1 --dag | dot -Tpng > dag_atac.png

# Local run
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ATAC/test_config_atac.yaml --cores 1

# Run with HPC profile
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ATAC/test_config_atac.yaml \
                             --cores all \
                             --snakemake-args --profile tests/00_TEST_DATA_GEX/HPC_profiles
```

#### ARC
```bash
# Dry run
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ARC/test_config_arc.yaml --cores 1 --dry-run

# DAG visualization
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ARC/test_config_arc.yaml --cores 1 --dag | dot -Tpng > dag_arc.png

# Local run
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ARC/test_config_arc.yaml --cores 1

# Run with HPC profile
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ARC/test_config_arc.yaml \
                             --cores all \
                             --snakemake-args --profile tests/00_TEST_DATA_GEX/HPC_profiles
```

### Testing Checklist
Before marking a feature complete:
- [ ] Dry-run succeeds (`--dry-run`)
- [ ] Full pipeline runs without errors
- [ ] Traceability metadata verified (`cell_id` unique)
- [ ] Output files created in correct locations
- [ ] Log files show expected output

### Verification Pattern
```python
import scanpy as sc

adata = sc.read_h5ad("output.h5ad")

# Required checks
assert 'batch_id' in adata.obs.columns
assert 'capture_id' in adata.obs.columns
assert 'cell_id' in adata.obs.columns
assert adata.obs['cell_id'].is_unique, "cell_id must be unique!"

print("✓ Verification passed")
```

---

## Common Issues (Don't Repeat)

### ATAC: SnapATAC2 API
- **Wrong:** `snap.read(file=...)`
- **Right:** `snap.pp.import_fragments(sorted_by_barcode=True)`

### ATAC: QC Filters
- Default `min_num_fragments=200` filters out test data
- Use `min_num_fragments=1` for testing
- Increase for production data

### Demux: Multi-Modality
- Use generic `COUNT_DIR` variable, not `GEX_COUNT_DIR`
- Detect modality: check `cellranger_gex`, `cellranger_atac`, `cellranger_arc`

### Batch Aggregation: Duplicate Warnings
- Warning about duplicate `obs_names` is expected (same barcodes in different captures)
- **Important:** `cell_id` MUST be unique (verify this!)

---

## Package API Reference

Exact API signatures for the versions pinned in this project. Use these to avoid confusing function names, argument names, or behaviors across packages.

---

### muon 0.1.7

**Import:** `import muon as mu`

#### I/O

```python
mu.read_10x_h5(filename: PathLike, extended: bool = True, *args, **kwargs) -> MuData
```
Reads a 10X HDF5 file (`.h5`). With `extended=True` (default) auto-locates peak annotations and fragment files. Use `extended=False` to skip — required when mixed feature types would cause write errors or when setting `fragments_file` manually.

```python
mu.read_10x_mtx(path: PathLike, extended: bool = True, *args, **kwargs) -> MuData
```
Reads a 10X MTX directory (`filtered_feature_bc_matrix/`). Wraps `sc.read_10x_mtx()`.
**CRITICAL — default modality names created by this function:**
| Cell Ranger `feature_types` value | Modality key in MuData |
|------------------------------------|------------------------|
| `"Gene Expression"`                | `"rna"` ← NOT `"gex"` |
| `"Peaks"`                          | `"atac"`               |
| `"Antibody Capture"`               | `"prot"`               |

This pipeline renames `"rna"` → `"gex"` immediately after reading:
```python
mdata.mod["gex"] = mdata.mod.pop("rna")
```

```python
mu.read(filename: PathLike) -> MuData
```
Reads a saved `.h5mu` file.

#### MuData Class

```python
mu.MuData(
    data: AnnData | Mapping[str, AnnData] | MuData | None = None,
    feature_types_names: Mapping[str, str] | None = {"Gene Expression": "rna", "Peaks": "atac", "Antibody Capture": "prot"},
    **kwargs,
)
```
Constructor. When passing a `dict[str, AnnData]`, the dict keys become modality names.

**Key attributes:**
- `mdata.mod` — `Mapping[str, AnnData]`, read-only dict of modalities
- `mdata["gex"]` — shorthand for `mdata.mod["gex"]`, returns the AnnData for that modality
- `mdata.obs` — top-level `pd.DataFrame` (shared observation metadata across modalities)
- `mdata.var` — top-level `pd.DataFrame` (all features concatenated across modalities)
- `mdata.n_obs` / `mdata.n_vars` — cell/feature counts at MuData level
- `mdata.uns` — unstructured metadata dict

**Key methods:**
```python
mdata.update()          # Sync global obs/var from modality-level data; call after modifying modalities
mdata.write(filename)   # Write to .h5mu (alias: mdata.write_h5mu(filename))
mdata.copy()            # Return an independent copy
mdata.pull_obs(key)     # Pull a column from modality obs into top-level mdata.obs
mdata.push_obs(key)     # Push a column from top-level mdata.obs into modality obs
```

**After `mu.MuData(dict)` or modality rename, always call `mdata.update()`** to re-sync the global obs index.

**Top-level obs after `mdata.update()`:** muon prefixes per-modality obs columns with the modality name (e.g., `"gex:batch_id"`). To expose columns without prefix, set them directly:
```python
mdata.obs["batch_id"] = mdata["gex"].obs["batch_id"]
```

#### Concatenation (ARC batch aggregation pattern)

muon does not have a top-level `mu.concat()`. For ARC, concatenate each modality separately with `anndata.concat()`, then reconstruct MuData:
```python
merged_mods = {
    mod_name: ad.concat([obj.mod[mod_name] for obj in objects], axis=0, join="outer", merge="same")
    for mod_name in objects[0].mod.keys()
}
batch_obj = mu.MuData(merged_mods)
batch_obj.update()
```

---

### scanpy ≥1.10 (GEX)

**Import:** `import scanpy as sc`

```python
sc.read_10x_h5(filename, genome=None, gex_only=True, backup_url=None) -> AnnData
```
- `gex_only=True` (default): returns only Gene Expression features
- `gex_only=False`: returns all feature types (used internally by muon)

```python
sc.read_h5ad(filename, backed=None) -> AnnData
```
Reads a saved `.h5ad` file.

```python
sc.pp.calculate_qc_metrics(adata, qc_vars=[], percent_top=(50,), log1p=True, inplace=False)
```
Computes per-cell (`obs`) and per-gene (`var`) QC metrics. Common `qc_vars`: `["mt", "ribo"]`.

---

### snapatac2 2.9.0 (ATAC)

**Import:** `import snapatac2 as snap`

```python
snap.pp.import_fragments(
    fragment_file: Path | list[Path],  # path to fragments.tsv.gz
    chrom_sizes: Genome | dict[str, int],  # chromosome sizes; or snap.genome.* presets
    *,
    is_paired: bool = True,            # True for paired-end (10X default)
    file: Path | list[Path] | None = None,  # output .h5ad path (None = in-memory)
    min_num_fragments: int = 200,      # REDUCE TO 1 FOR TEST DATA
    sorted_by_barcode: bool = True,    # True after pre-sorting (2.9.0 default is True)
    whitelist: Path | list[str] | None = None,
    chrM: list[str] = ["chrM", "M"],
    chunk_size: int = 2000,
    tempdir: Path | None = None,
    backend: Literal["hdf5"] = "hdf5",
    n_jobs: int = 8,
) -> AnnData
```
**WRONG:** `snap.read(file=...)` — this does not exist.
**RIGHT:** `snap.pp.import_fragments(fragment_file=..., chrom_sizes=...)` — always use this.

Note: In 2.9.0 `sorted_by_barcode` defaults to `True` (unlike older versions). The pipeline pre-sorts and caches as `fragments.sorted_by_barcode.tsv.gz`, then passes `sorted_by_barcode=True`.

---

### anndata ≥0.10

**Import:** `import anndata as ad`

```python
ad.concat(
    adatas: Collection[AnnData] | Mapping[str, AnnData],
    axis: Literal[0, 1] = 0,          # 0=obs (cells), 1=var (features)
    join: Literal["inner", "outer"] = "inner",
    merge: str | None = None,          # how to merge .uns; "same" keeps keys with identical values
    uns_merge: str | None = None,
    label: str | None = None,          # column name added to .obs for source key
    keys: Collection | None = None,    # values for label column
    index_unique: str | None = None,   # separator to make obs_names unique; None keeps originals
    fill_value = None,                 # fill for outer join missing values
    pairwise: bool = False,
) -> AnnData
```
Duplicate `obs_names` warning during batch aggregation is **expected** (same barcodes in different captures). `cell_id` uniqueness is what matters — always verify after concat.

---

## What NOT to Do

### ❌ Never
- Add R-based analysis tools
- Skip traceability metadata
- Hardcode GEX-specific paths in multi-modality code
- Create objects at batch-level instead of per-capture (unless aggregating)
- Skip `cell_id` uniqueness verification
- Over-engineer simple solutions
- Add features not explicitly requested

### ⚠️ Ask First
- Major architectural changes
- Renumbering directory structure
- Changing file naming conventions
- Adding new dependencies
- Modifying traceability metadata format

---

## Git Workflow

### Commits
- Use conventional commits format: `feat:`, `fix:`, `refactor:`, etc.
- Include co-author tag:
  ```
  Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
  ```

### Branches
- Work on feature branches when appropriate
- Keep `main` stable

---

## Documentation

### Code Comments
- Only add comments where logic isn't self-evident
- Don't add docstrings to unchanged code
- Focus comments on "why" not "what"

### Log Messages
- Use clear, informative print statements in scripts
- Include progress indicators for long operations
- Log file locations and key statistics

---

## Performance Considerations

### ATAC Data
- Always sort fragments by barcode before import
- Cache sorted files (don't re-sort on every run)
- Use `sorted_by_barcode=True` for efficiency

### Large Datasets
- Consider using backed mode for very large objects
- Monitor memory usage with multi-capture aggregation

---