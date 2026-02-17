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
- `scanpy>=1.10` (GEX)
- `anndata>=0.10`
- `muon>=0.1.6` (ARC)
- `snapatac2>=2.6` (ATAC)

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
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ARC/test_config_gex.yaml --cores 1

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
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_ARC/test_config_atac.yaml --cores 1

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