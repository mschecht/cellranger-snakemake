# Development Guide

This guide walks you through contributing to the cellranger-snakemake pipeline. Whether you're adding a new analysis method or an entirely new pipeline step, this document provides the workflow, patterns, and testing strategies you need.

## Key Concepts

**Step**: A preprocessing stage in single-cell analysis, such as demultiplexing, doublet detection, or cell type annotation. Each step has its own Snakemake rule file (`workflows/rules/<step>.smk`) and Pydantic schema (`schemas/<step>.py`).

**Method**: A software tool that implements a step. For example, the demultiplexing step supports two methods: `vireo` and `demuxalot`. Adding a new method means integrating another tool into an existing step.

## Developer Workflow

Follow these steps when making changes:

1. **Orient** - Understand the [project architecture](#project-architecture) and how the pipeline resolves what to run
2. **Develop** - Choose your task: [add a method](#adding-a-method-to-an-existing-step) or [add a step](#adding-a-new-pipeline-step)
3. **Test** - Validate with [DAG checks](#dag-validation) and [integration tests](#integration-tests)
4. **Document** - Update [documentation](#building-and-editing-documentation) as needed

---

## Project Architecture

Understanding the file layout is essential before making changes:

```
cellranger_snakemake/
├── cli.py                          # CLI entry points
├── config_generator.py             # Interactive config builder
├── config_validator.py             # PIPELINE_DIRECTORIES, validation
├── schemas/                        # Pydantic models for config validation
│   ├── base.py                     # BaseStepConfig (all steps inherit this)
│   ├── config.py                   # PipelineConfig (unified schema)
│   ├── cellranger.py               # GEX/ATAC/ARC Cell Ranger configs
│   ├── demultiplexing.py           # DemuxalotConfig, VireoConfig
│   ├── doublet_detection.py        # ScrubletConfig, SoloConfig
│   └── annotation.py               # CelltypistConfig, ScANVIConfig, DecouplerMarkerConfig
├── workflows/
│   ├── main.smk                    # Master workflow, rule all, includes
│   ├── rules/                      # One .smk file per pipeline step
│   │   ├── cellranger.smk          # Cell Ranger count/aggregation
│   │   ├── object_creation.smk     # Per-capture AnnData/MuData creation
│   │   ├── batch_aggregation.smk   # Batch-level aggregation + metadata enrichment
│   │   ├── demultiplexing.smk      # Demuxalot/Vireo
│   │   ├── doublet_detection.smk   # Scrublet/SOLO
│   │   └── celltype_annotation.smk # Celltypist/scANVI/Decoupler
│   └── scripts/
│       ├── build_targets.py        # Generates target files for rule all
│       ├── parse_config.py         # Extracts enabled steps, methods, etc.
│       ├── create_gex_anndata.py   # GEX per-capture object creation
│       ├── create_atac_anndata.py  # ATAC per-capture object creation
│       ├── create_arc_mudata.py    # ARC per-capture MuData creation
│       ├── aggregate_batch.py      # Batch aggregation (per-capture → batch)
│       ├── merge_metadata.py       # Merge analysis metadata into batch objects
│       ├── run_scrublet.py         # Scrublet doublet detection
│       ├── run_solo.py             # SOLO doublet detection
│       ├── run_celltypist.py       # Celltypist annotation
│       ├── run_scanvi.py           # scANVI annotation
│       └── run_decoupler_markers.py # Decoupler marker-based annotation
tests/
├── test.sh                         # Integration test script
├── 00_TEST_DATA_GEX/               # Test configs and library lists
└── ...
docs/source/                        # Read the Docs documentation (this file)
```

### How the pipeline resolves what to run

1. **Config** (`pipeline_config.yaml`) declares which steps are `enabled: true`
2. **`parse_config.py`** → `get_enabled_steps()` reads the config and returns a list of enabled step names
3. **`main.smk`** conditionally includes `.smk` rule files based on enabled steps
4. **`build_targets.py`** → `build_all_targets()` generates the list of expected `.done` files for `rule all`
5. **Each rule** produces a `.done` marker file in `{output_dir}/00_LOGS/` that matches what `build_targets.py` expects

If the target filename from `build_targets.py` doesn't match the `done` output in the rule, Snakemake will raise a `MissingInputException`.

---

## Development Workflow

### Quick Reference Checklists

**Adding a method to an existing step** (e.g., a new demultiplexing tool):

1. Add method config schema in `schemas/<step>.py` with `tool_meta`
2. Register the method in the parent config class
3. Add the rule in `workflows/rules/<step>.smk`
4. Add target generation in `workflows/scripts/build_targets.py`
5. Add to config generator in `config_generator.py`
6. [Test](#testing)
7. [Document](#building-and-editing-documentation)

**Adding a new pipeline step** (e.g., cell type annotation):

1. Create the Pydantic schema in `schemas/<new_step>.py`
2. Register the output directory in `config_validator.py`
3. Register the step in `parse_config.py`
4. Add target generation in `build_targets.py`
5. Create the rule file in `workflows/rules/<new_step>.smk`
6. Include the rule file in `main.smk`
7. Create a dummy rule and test the DAG
8. Implement the rule
9. Add to config generator
10. [Write tests](#testing)
11. [Document](#building-and-editing-documentation)

---

### Adding a Method to an Existing Step

This example shows how [Vireo](https://github.com/single-cell-genetics/vireo) was added alongside demuxalot for demultiplexing. Use this as a template for adding new methods.

#### Step 1: Add method config schema

In `cellranger_snakemake/schemas/demultiplexing.py`:

```python
from typing import ClassVar
from .base import ToolMeta

class CellSNPConfig(BaseModel):
    """cellsnp-lite parameters for SNP calling."""

    vcf: str = Field(description="Path to VCF reference file with known variants")
    threads: int = Field(default=4, ge=1, description="Number of threads for cellsnp-lite")
    min_maf: float = Field(default=0.0, ge=0.0, le=1.0, description="Minimum minor allele frequency")
    min_count: int = Field(default=1, ge=0, description="Minimum UMI count")

    class Config:
        extra = "forbid"


class VireoConfig(BaseModel):
    """Vireo demultiplexing parameters (requires cellsnp-lite preprocessing)."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="vireoSNP",
        url="https://github.com/single-cell-genetics/vireo",
    )

    cellsnp: CellSNPConfig = Field(description="cellsnp-lite configuration for SNP calling")
    donors: int = Field(description="Number of donors to demultiplex")

    class Config:
        extra = "forbid"
```

Every method schema must include a `tool_meta` class variable. This enables `show-params` to display the installed tool version and a link to the source. Use `ClassVar` so Pydantic treats it as a class attribute (not a config field). For tools that aren't Python packages (e.g., Cell Ranger), set `shell_version_cmd`:

```python
tool_meta: ClassVar[ToolMeta] = ToolMeta(
    package="vireo",
    url="https://github.com/single-cell-genetics/vireo/releases/tag/v0.2.3",
    shell_version_cmd="cellranger --version",
)
```

#### Step 2: Register the method in the parent config class

```python
class DemultiplexingConfig(BaseStepConfig):
    method: Literal["demuxalot", "vireo"] = Field(...)

    vireo: Optional[VireoConfig] = None  # Add this

    @model_validator(mode='after')
    def validate_method_params(self):
        method_configs = {
            "demuxalot": self.demuxalot,
            "vireo": self.vireo,  # Add this
        }
        # ... rest stays the same
```

#### Step 3: Add the rule

In `cellranger_snakemake/workflows/rules/demultiplexing.smk`:

```python
if config.get("demultiplexing") and DEMUX_METHOD == "vireo":

    # Parse vireo config
    VIREO_CONFIG = DEMUX_CONFIG.get("vireo", {})
    CELLSNP_CONFIG = VIREO_CONFIG.get("cellsnp", {})
    VIREO_DONORS = VIREO_CONFIG.get("donors")

    # cellsnp-lite parameters
    CELLSNP_VCF = CELLSNP_CONFIG.get("vcf")
    CELLSNP_THREADS = CELLSNP_CONFIG.get("threads", 4)

    rule cellsnp_lite:
        """Run cellsnp-lite for SNP calling from BAM."""
        input:
            gex_done = os.path.join(config.get("output_dir", "output"), "00_LOGS", "{batch}_{capture}_gex_count.done"),
            bam = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "possorted_genome_bam.bam"),
            barcodes = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
        output:
            base_vcf = os.path.join(DEMUX_OUTPUT_DIR, "cellsnp_output_{batch}_{capture}", "cellSNP.base.vcf.gz"),
            done = touch(os.path.join(config.get("output_dir", "output"), "00_LOGS", "cellsnp_output_{batch}_{capture}.done"))
        # ...

    rule vireo:
        """Run Vireo for donor deconvolution using cellsnp-lite output."""
        input:
            cellsnp_done = rules.cellsnp_lite.output.done
        output:
            donor_ids = os.path.join(DEMUX_OUTPUT_DIR, "vireo_output_{batch}_{capture}", "donor_ids.tsv"),
            done = touch(os.path.join(OUTPUT_DIRS["logs_dir"], "vireo_output_{batch}_{capture}.done"))
        # ...
```

**Key conventions:**
- The `.done` filename pattern **must** match what `build_targets.py` generates: `{method}_output_{batch}_{capture}.done`
- `.done` files always go in `{output_dir}/00_LOGS/`
- Use `touch()` for `.done` outputs
- Shared variables (output dirs, GEX count dir) go in the top-level `if config.get(...)` block
- Method-specific config goes in the method-level `if` block

#### Step 4: Add target generation

In `cellranger_snakemake/workflows/scripts/build_targets.py`:

```python
if method == "vireo":
    for batch in batches:
        for capture in captures:
            outputs.append(os.path.join(logs_dir, f"vireo_output_{batch}_{capture}.done"))
```

#### Step 5: Add to config generator

Update `cellranger_snakemake/config_generator.py` so `init-config` can produce the new method's parameters interactively.

#### Step 6: Test

See [Testing](#testing) below.

---

### Adding a New Pipeline Step

This example shows how cell type annotation was added as a pipeline step. Use this as a template for adding new steps.

#### Step 1: Create the Pydantic schema

Create `cellranger_snakemake/schemas/annotation.py`:

```python
"""Cell type annotation configuration schemas."""

from typing import ClassVar, Literal, Optional
from pydantic import BaseModel, Field, model_validator
from .base import BaseStepConfig, ToolMeta


class CelltypistConfig(BaseModel):
    """Celltypist annotation parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="celltypist",
        url="https://github.com/Teichlab/celltypist",
    )

    model: str = Field(description="Path to celltypist model file or model name")
    majority_voting: bool = Field(default=False, description="Use majority voting")

    class Config:
        extra = "forbid"


class CelltypeAnnotationConfig(BaseStepConfig):
    """Cell type annotation step configuration."""

    method: Literal["celltypist", "azimuth", "singler", "sctype"] = Field(
        description="Cell type annotation method to use"
    )
    celltypist: Optional[CelltypistConfig] = None
    # ... other method configs ...

    @model_validator(mode='after')
    def validate_method_params(self):
        # Validate that the selected method has its config block
        ...
```

#### Step 2: Register the output directory

In `cellranger_snakemake/config_validator.py`, add to `PIPELINE_DIRECTORIES`:

```python
PIPELINE_DIRECTORIES = [
    ("logs", "00_LOGS"),
    # ... existing entries ...
    ("celltype_annotation", "05_CELLTYPE_ANNOTATION"),
]
```

#### Step 3: Register the step in parse_config.py

In `cellranger_snakemake/workflows/scripts/parse_config.py`, add the step name to the `get_enabled_steps` list:

```python
for step in ["cellranger_gex", "cellranger_atac", "cellranger_arc",
             "demultiplexing", "doublet_detection", "celltype_annotation"]:
```

#### Step 4: Add target generation in build_targets.py

In `cellranger_snakemake/workflows/scripts/build_targets.py`:

1. Add a call in `build_all_targets()`:

```python
if "celltype_annotation" in enabled_steps:
    targets.extend(get_annotation_outputs(config))
```

2. Add the target function:

```python
def get_annotation_outputs(config):
    if not config.get("celltype_annotation"):
        return []

    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    annot_config = config["celltype_annotation"]
    method = annot_config["method"]

    outputs = []
    if config.get("cellranger_gex"):
        df = pd.read_csv(config["cellranger_gex"]["libraries"], sep="\t")
        batches = df['batch'].unique().tolist()
        captures = df['capture'].unique().tolist()

        for batch in batches:
            for capture in captures:
                outputs.append(os.path.join(logs_dir, f"{method}_output_{batch}_{capture}.done"))

    return outputs
```

#### Step 5: Create the rule file

Create `cellranger_snakemake/workflows/rules/celltype_annotation.smk`:

```python
"""Cell type annotation workflow rules."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger

if config.get("celltype_annotation"):
    ANNOT_CONFIG = config["celltype_annotation"]
    ANNOT_METHOD = ANNOT_CONFIG["method"]

    custom_logger.info(f"Cell Type Annotation: Using {ANNOT_METHOD} method")

# ============================================================================
# CELLTYPIST
# ============================================================================

if config.get("celltype_annotation") and ANNOT_METHOD == "celltypist":

    rule celltypist:
        """Run Celltypist for cell type annotation."""
        input:
            h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
        output:
            predictions = "{sample}/celltypist/predicted_labels.csv",
            done = touch("{sample}/celltypist/{sample}_celltypist.done")
        params:
            model = ANNOT_CONFIG.get("celltypist", {}).get("model", "Immune_All_Low.pkl")
        script:
            "../scripts/run_celltypist.py"
```

#### Step 6: Include the rule file in main.smk

In `cellranger_snakemake/workflows/main.smk`:

```python
if "celltype_annotation" in ENABLED_STEPS:
    include: "rules/celltype_annotation.smk"
```

#### Step 7: Create a dummy rule and test the DAG

Before implementing the actual tool logic, write a dummy `shell` block:

```python
shell:
    """
    echo "Placeholder for celltypist"
    touch {output.predictions}
    """
```

Then verify the DAG resolves correctly:

```bash
snakemake-run-cellranger run --config-file your_test_config.yaml --cores 1 --dry-run
snakemake-run-cellranger run --config-file your_test_config.yaml --cores 1 --dag | dot -Tpng > dag.png
```

Check that:
- Your new rule appears in the DAG
- It depends on the correct upstream rules (e.g., `cellranger_gex_count`)
- `rule all` connects to your new rule's `.done` output
- No `MissingInputException` errors

#### Step 8: Implement the rule

Replace the dummy shell with the actual tool invocation. Use either `shell:` for command-line tools or `run:` for Python-based tools.

#### Step 9: Add to config generator

Update `cellranger_snakemake/config_generator.py` so that `snakemake-run-cellranger init-config` can interactively generate config for the new step.

#### Step 10: Write tests

See [Testing](#testing) below.

#### Step 11: Document

Update this file and the tutorial if the new step is relevant to the standard user workflow.

---

## Testing

Validate your changes before merging.

### DAG validation

Always verify the DAG first. This catches target/output mismatches without executing any rules:

```bash
# Dry run - checks all inputs/outputs resolve
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dry-run

# Visual DAG - confirm rule dependencies look correct
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dag | dot -Tpng > dag.png
```

### Integration tests

Integration tests run the full pipeline on test data:

```bash
# Run integration tests
bash tests/test.sh

# Or run a specific workflow manually
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1
```

### Test checklist for a new step

When adding a new step, verify all of the following before merging:

- [ ] **Config validation**: `snakemake-run-cellranger validate-config --config-file your_config.yaml` succeeds
- [ ] **Dry run**: `snakemake-run-cellranger run --config-file ... --cores 1 --dry-run` shows your rule
- [ ] **DAG**: Your rule appears with correct dependencies in the DAG visualization
- [ ] **Dummy execution**: Pipeline completes with placeholder `touch` commands
- [ ] **Real execution**: Pipeline completes with the actual tool on test data
- [ ] **Config generator**: `snakemake-run-cellranger init-config` includes the new step

---

## Common Pitfalls

| Problem | Cause | Fix |
|---|---|---|
| `MissingInputException` for `.done` files | Target filename in `build_targets.py` doesn't match rule `output.done` | Ensure both use the exact same pattern, e.g., `{method}_output_{batch}_{capture}.done` |
| `NameError` for shared variables | Variable defined in one method block but used in another | Move shared variables (output dirs, GEX count dir) to the common `if config.get("step"):` block |
| Pydantic `ValidationError` on valid config | Schema too strict or missing `Optional` on method-specific config | Method configs should be `Optional[...] = None` with a `model_validator` to enforce presence based on `method` |
| Rule not in DAG | Step not added to `get_enabled_steps()` or `main.smk` includes | Check `parse_config.py` and `main.smk` both reference your step name |
| Config generator skips new step | Step not added to `config_generator.py` | Add interactive prompts for the new step's parameters |

---

## Building and Editing Documentation

This project uses [Sphinx](https://www.sphinx-doc.org/) with MyST-Parser for markdown support and is deployed on [Read the Docs](https://about.readthedocs.com/). Documentation is hosted at: https://cellranger-snakemake.readthedocs.io/

### Documentation structure

```
docs/
├── source/
│   ├── conf.py           # Sphinx configuration
│   ├── index.md          # Landing page with toctree
│   ├── installation.md
│   ├── quickstart.md
│   ├── tutorial.md
│   └── development.md    # This file
├── requirements.txt      # Sphinx dependencies
└── Makefile              # Build commands
```

### Building docs locally

Render the documentation locally in your web browser:

```bash
cd docs

# Build HTML
make html

# Serve locally
python3 -m http.server 8000 -d build/html
```

Then open http://localhost:8000 in your browser.

If you are using VS Code on a remote session, you can render the docs in the IDE itself.

### Live reload during editing

For automatic rebuilds when you save changes, use `sphinx-autobuild`:

```bash
pip install sphinx-autobuild
sphinx-autobuild source build/html --port 8000
```

This watches `source/` for changes, rebuilds automatically, and refreshes your browser.

### Adding a new page

1. Create a new `.md` file in `docs/source/`
2. Add it to the `toctree` in `docs/source/index.md`:

```markdown
```{toctree}
:maxdepth: 2
:caption: Contents

installation
quickstart
tutorial
your-new-page
development
```
