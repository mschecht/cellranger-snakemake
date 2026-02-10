# Development Guide

This guide covers two audiences: **users** who want to modify existing workflow parameters, and **developers** who want to add entirely new pipeline steps or rules.

## Project Architecture

Understanding the file layout is essential before making changes:

```
cellranger_snakemake/
├── cli.py                          # CLI entry points
├── config_generator.py             # Interactive config builder
├── config_validator.py             # PIPELINE_DIRECTORIES, validation
├── schemas/                        # Pydantic models for config validation
│   ├── base.py                     # BaseStepConfig (all steps inherit this)
│   ├── demultiplexing.py           # Example: DemuxalotConfig, VireoConfig
│   ├── doublet_detection.py
│   └── annotation.py
├── workflows/
│   ├── main.smk                    # Master workflow, rule all, includes
│   ├── rules/                      # One .smk file per pipeline step
│   │   ├── cellranger.smk
│   │   ├── demultiplexing.smk
│   │   ├── doublet_detection.smk
│   │   └── celltype_annotation.smk
│   └── scripts/
│       ├── build_targets.py        # Generates target files for rule all
│       └── parse_config.py         # Extracts enabled steps, methods, etc.
tests/
├── test.sh                         # Integration test script
├── 00_TEST_DATA_GEX/               # Test configs and library lists
└── ...
docs/source/                        # Read the Docs documentation (this file)
docs/source/development.md          # Dev documentation (this file)
```

### How the pipeline resolves what to run

1. **Config** (`pipeline_config.yaml`) declares which steps are `enabled: true`
2. **`parse_config.py`** → `get_enabled_steps()` reads the config and returns a list of enabled step names
3. **`main.smk`** conditionally includes `.smk` rule files based on enabled steps
4. **`build_targets.py`** → `build_all_targets()` generates the list of expected `.done` files for `rule all`
5. **Each rule** produces a `.done` marker file in `{output_dir}/00_LOGS/` that matches what `build_targets.py` expects

If the target filename from `build_targets.py` doesn't match the `done` output in the rule, Snakemake will raise a `MissingInputException`.

---

## For Users: Modifying Existing Steps

### Changing parameters for an existing method

All parameters are defined in the YAML config. To see what's available:

```bash
# Show all parameters for a specific method
snakemake-run-cellranger show-params

# Generate a default config with all options
snakemake-run-cellranger init-config --get-default-config
```

Edit the generated YAML and validate:

```bash
snakemake-run-cellranger validate-config --config-file your_config.yaml
```

### Adding a new method to an existing step

For example, adding a new demultiplexing method called `souporcell` alongside the existing `vireo` and `demuxalot`:

1. **Add method config schema** in `cellranger_snakemake/schemas/demultiplexing.py`:

```python
from typing import ClassVar
from .base import ToolMeta

class SouporcellConfig(BaseModel):
    """Souporcell demultiplexing parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="souporcell",
        url="https://github.com/wheaton5/souporcell",
    )

    num_clusters: int = Field(description="Expected number of clusters/donors")
    # ... additional params

    class Config:
        extra = "forbid"
```

Every method schema must include a `tool_meta` class variable. This enables `show-params` to display the installed tool version and a link to the source. Use `ClassVar` so Pydantic treats it as a class attribute (not a config field). For tools that aren't Python packages (e.g., Cell Ranger), set `shell_version_cmd`:

```python
tool_meta: ClassVar[ToolMeta] = ToolMeta(
    package="cellranger",
    url="https://www.10xgenomics.com/support/software/cell-ranger/latest",
    shell_version_cmd="cellranger --version",
)
```

2. **Register the method** in the parent config class in the same file:

```python
class DemultiplexingConfig(BaseStepConfig):
    method: Literal["demuxalot", "vireo", "souporcell"] = Field(...)

    souporcell: Optional[SouporcellConfig] = None  # Add this

    @model_validator(mode='after')
    def validate_method_params(self):
        method_configs = {
            "demuxalot": self.demuxalot,
            "vireo": self.vireo,
            "souporcell": self.souporcell,  # Add this
        }
        # ... rest stays the same
```

3. **Add the rule** in `cellranger_snakemake/workflows/rules/demultiplexing.smk`:

```python
if config.get("demultiplexing") and DEMUX_METHOD == "souporcell":
    # Method-specific config only
    SOUPORCELL_CONFIG = DEMUX_CONFIG.get("souporcell", {})
    # ...

    rule souporcell:
        input:
            # Use shared GEX_COUNT_DIR for BAM/barcodes
            bam = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "possorted_genome_bam.bam"),
            barcodes = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
        output:
            # .done file MUST match the pattern in build_targets.py
            done = touch(os.path.join(config.get("output_dir", "output"), "00_LOGS", "souporcell_output_{batch}_{capture}.done"))
        # ...
```

4. **Add target generation** in `cellranger_snakemake/workflows/scripts/build_targets.py`:

```python
if method == "souporcell":
    for batch in batches:
        for capture in captures:
            outputs.append(os.path.join(logs_dir, f"souporcell_output_{batch}_{capture}.done"))
```

5. **Add to config generator** in `cellranger_snakemake/config_generator.py` so `init-config` can produce the new method's parameters interactively.

6. **Test** (see [Testing](#testing) below).

---

## For Developers: Adding a New Pipeline Step

This is the complete checklist for adding a new step (e.g., `ambient_rna_removal`). Follow these in order.

### Step 1: Create the Pydantic schema

Create `cellranger_snakemake/schemas/ambient_rna.py`:

```python
"""Ambient RNA removal configuration schemas."""

from typing import ClassVar, Literal, Optional
from pydantic import BaseModel, Field
from .base import BaseStepConfig, ToolMeta


class CellBenderConfig(BaseModel):
    """CellBender ambient RNA removal parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellbender",
        url="https://github.com/broadinstitute/CellBender",
    )

    epochs: int = Field(default=150, description="Number of training epochs")
    fpr: float = Field(default=0.01, description="False positive rate")

    class Config:
        extra = "forbid"


class AmbientRNAConfig(BaseStepConfig):
    """Ambient RNA removal step configuration."""

    method: Literal["cellbender"] = Field(description="Method to use")
    cellbender: Optional[CellBenderConfig] = None

    @model_validator(mode='after')
    def validate_method_params(self):
        # Same validation pattern as other steps
        ...
```

### Step 2: Register the output directory

In `cellranger_snakemake/config_validator.py`, add to `PIPELINE_DIRECTORIES`:

```python
PIPELINE_DIRECTORIES = [
    ("logs", "00_LOGS"),
    # ... existing entries ...
    ("ambient_rna", "06_AMBIENT_RNA"),  # Add new entry
]
```

### Step 3: Register the step in parse_config.py

In `cellranger_snakemake/workflows/scripts/parse_config.py`, add the step name to the `get_enabled_steps` list:

```python
for step in ["cellranger_gex", "cellranger_atac", "cellranger_arc",
             "demultiplexing", "doublet_detection", "celltype_annotation",
             "ambient_rna_removal"]:  # Add here
```

### Step 4: Add target generation in build_targets.py

In `cellranger_snakemake/workflows/scripts/build_targets.py`:

1. Add a call in `build_all_targets()`:

```python
if "ambient_rna_removal" in enabled_steps:
    targets.extend(get_ambient_rna_outputs(config))
```

2. Add the target function:

```python
def get_ambient_rna_outputs(config):
    if not config.get("ambient_rna_removal"):
        return []

    output_dirs = parse_output_directories(config)
    logs_dir = output_dirs["logs_dir"]
    method = config["ambient_rna_removal"]["method"]

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

### Step 5: Create the rule file

Create `cellranger_snakemake/workflows/rules/ambient_rna.smk`. Follow the existing pattern:

```python
"""Ambient RNA removal workflow rules."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger
from cellranger_snakemake.config_validator import parse_output_directories

if config.get("ambient_rna_removal"):
    AMBIENT_CONFIG = config["ambient_rna_removal"]
    AMBIENT_METHOD = AMBIENT_CONFIG["method"]

    # Shared setup
    OUTPUT_DIRS = parse_output_directories(config)
    AMBIENT_OUTPUT_DIR = OUTPUT_DIRS["ambient_rna_dir"]

    if config.get("cellranger_gex"):
        GEX_COUNT_DIR = os.path.join(config.get("output_dir", "output"), "01_CELLRANGERGEX_COUNT")

    custom_logger.info(f"Ambient RNA removal: Using {AMBIENT_METHOD} method")

# ============================================================================
# CELLBENDER
# ============================================================================

if config.get("ambient_rna_removal") and AMBIENT_METHOD == "cellbender":

    CELLBENDER_CONFIG = AMBIENT_CONFIG.get("cellbender", {})
    # ... method-specific params ...

    rule cellbender:
        """Run CellBender for ambient RNA removal."""
        input:
            gex_done = os.path.join(config.get("output_dir", "output"), "00_LOGS", "{batch}_{capture}_gex_count.done"),
            raw_h5 = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "raw_feature_bc_matrix.h5")
        output:
            filtered = os.path.join(AMBIENT_OUTPUT_DIR, "{batch}_{capture}", "cellbender_filtered.h5"),
            done = touch(os.path.join(config.get("output_dir", "output"), "00_LOGS", "cellbender_output_{batch}_{capture}.done"))
        # ...
```

**Key conventions:**
- The `.done` filename pattern **must** match what `build_targets.py` generates: `{method}_output_{batch}_{capture}.done`
- `.done` files always go in `{output_dir}/00_LOGS/`
- Use `touch()` for `.done` outputs
- Shared variables (output dirs, GEX count dir) go in the top-level `if config.get(...)` block
- Method-specific config goes in the method-level `if` block

### Step 6: Include the rule file in main.smk

In `cellranger_snakemake/workflows/main.smk`:

```python
if "ambient_rna_removal" in ENABLED_STEPS:
    include: "rules/ambient_rna.smk"
```

### Step 7: Create a dummy rule and test the DAG

Before implementing the actual tool logic, write a dummy `shell` block:

```python
shell:
    """
    echo "Placeholder for cellbender"
    touch {output.filtered}
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

### Step 8: Implement the rule

Replace the dummy shell with the actual tool invocation. Use either `shell:` for command-line tools or `run:` for Python-based tools.

### Step 9: Add to config generator

Update `cellranger_snakemake/config_generator.py` so that `snakemake-run-cellranger init-config` can interactively generate config for the new step.

### Step 10: Write tests

See [Testing](#testing) below.

### Step 11: Document

Update this file and the tutorial if the new step is relevant to the standard user workflow.

---

## Testing

### DAG validation (quick smoke test)

Always verify the DAG before running anything. This catches target/output mismatches without executing any rules:

```bash
# Dry run - checks all inputs/outputs resolve
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dry-run

# Visual DAG - confirm rule dependencies look correct
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dag | dot -Tpng > dag.png
```

### Integration tests

Integration tests run the full pipeline on test data. Use the automated script or run manually:

```bash
# Run all integration tests
bash tests/test.sh

# Or run a specific workflow manually
snakemake-run-cellranger run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1
```

### Unit tests with pytest

Unit tests validate individual Python components (schemas, config parsing, target generation) without running Snakemake. Tests live in `tests/` and use pytest.

#### Config validation tests

Test that your Pydantic schema accepts valid config and rejects invalid config:

```python
# tests/test_schemas.py
import pytest
from cellranger_snakemake.schemas.demultiplexing import (
    DemultiplexingConfig,
    DemuxalotConfig,
)
from pydantic import ValidationError


def test_valid_demuxalot_config():
    """Valid demuxalot config should pass validation."""
    config = DemultiplexingConfig(
        enabled=True,
        method="demuxalot",
        demuxalot=DemuxalotConfig(
            vcf="/path/to/file.vcf.gz",
            genome_names="/path/to/names.txt",
            refine=True,
        ),
    )
    assert config.method == "demuxalot"
    assert config.demuxalot.refine is True


def test_missing_method_params():
    """Selecting a method without providing its config block should fail."""
    with pytest.raises(ValidationError):
        DemultiplexingConfig(
            enabled=True,
            method="demuxalot",
            # Missing demuxalot config block
        )


def test_extra_fields_rejected():
    """Unknown parameters should be rejected."""
    with pytest.raises(ValidationError):
        DemuxalotConfig(
            vcf="/path/to/file.vcf.gz",
            genome_names="/path/to/names.txt",
            refine=True,
            unknown_param="bad",
        )
```

#### Target generation tests

Test that `build_targets.py` produces the correct `.done` file paths:

```python
# tests/test_build_targets.py
import pytest
from cellranger_snakemake.workflows.scripts.build_targets import get_demux_outputs
from unittest.mock import patch
import pandas as pd


@pytest.fixture
def base_config(tmp_path):
    """Create a minimal config with a libraries file."""
    libraries = tmp_path / "libraries.tsv"
    libraries.write_text("batch\tcapture\tsample\tfastqs\n1\tL001\tS1\t/fastqs\n1\tL002\tS2\t/fastqs\n")
    return {
        "output_dir": str(tmp_path / "output"),
        "cellranger_gex": {
            "enabled": True,
            "libraries": str(libraries),
            "reference": "/ref",
        },
    }


def test_demuxalot_targets(base_config):
    """Demuxalot should produce one .done per batch-capture combo."""
    base_config["demultiplexing"] = {"enabled": True, "method": "demuxalot"}
    outputs = get_demux_outputs(base_config)
    assert len(outputs) == 2
    assert any("demuxalot_output_1_L001.done" in o for o in outputs)
    assert any("demuxalot_output_1_L002.done" in o for o in outputs)


def test_vireo_targets(base_config):
    """Vireo should produce one .done per batch-capture combo."""
    base_config["demultiplexing"] = {"enabled": True, "method": "vireo"}
    outputs = get_demux_outputs(base_config)
    assert len(outputs) == 2
    assert any("vireo_output_1_L001.done" in o for o in outputs)


def test_disabled_demux_returns_empty():
    """Disabled demultiplexing should return no targets."""
    assert get_demux_outputs({}) == []
```

#### Config validator tests

```python
# tests/test_config_validator.py
from cellranger_snakemake.config_validator import parse_output_directories


def test_output_directories_contain_all_steps():
    """All pipeline steps should have a directory mapping."""
    config = {"output_dir": "test_output"}
    dirs = parse_output_directories(config)
    assert "logs_dir" in dirs
    assert "demultiplexing_dir" in dirs
    assert dirs["logs_dir"] == "test_output/00_LOGS"
```

#### Running unit tests

```bash
# Run all unit tests
pytest tests/ -v

# Run a specific test file
pytest tests/test_schemas.py -v

# Run with coverage
pytest tests/ --cov=cellranger_snakemake --cov-report=term-missing
```

### Test checklist for a new step

When adding a new step, verify all of the following before merging:

- [ ] **Schema validation**: `pytest tests/test_schemas.py` passes
- [ ] **Target generation**: `pytest tests/test_build_targets.py` passes
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

This project uses [Sphinx](https://www.sphinx-doc.org/) with MyST-Parser for markdown support and is deployed on [Read the docs](https://about.readthedocs.com/). Documentation rendered documentation is hosted on here: https://cellranger-snakemake.readthedocs.io/.

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

This is how you can render the documentation locally on your web browser: 

```bash
cd docs

# Build HTML
make html

# Serve locally
python3 -m http.server 8000 -d build/html
```

Then open http://localhost:8000 in your browser.

If you are using VS Code on a remote sesions, it is super conventient to render the docs in the IDE itself! 

### Live reload during editing

For automatic rebuilds when you save changes, use `sphinx-autobuild`:

```bash
pip install sphinx-autobuild
sphinx-autobuild source build/html --port 8000
```

This watches `source/` for changes, rebuilds automatically, and refreshes your browser (EXTREMELY CONVENTIENT).

### Adding a new page

1. Create a new `.md` file in `docs/source/`
2. Add it to the `toctree` in `docs/source/index.md`: