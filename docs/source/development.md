# Development Guide

This guide walks you through contributing to the sc-preprocess pipeline. Whether you're adding a new analysis method or an entirely new single-cell preprocessing step, this document provides the workflow, patterns, and testing strategies you need.

## Getting Started

Follow the [Installation](installation.md) instructions (steps 1-5), using `pip install -e .` for an editable install.

After pulling new changes from the repository:
```bash
git pull
pip install -e .   
```

### Building the documentation

To build and preview the docs locally, install the documentation dependencies:

```bash
pip install -e ".[docs]"
```

Then serve with auto-reload:

```bash
sphinx-autobuild docs/source docs/build/html
```

### Updating dependencies

The project uses three types of environment files:

| File | Purpose | Editable? |
|------|---------|-----------|
| `environment.yaml` | Minimal install spec (hand-maintained) | Yes — edit this when adding/changing deps |
| `environment_lock_linux.yaml` | Exact pins for Linux reproducibility | No — regenerate with `conda env export` |
| `environment_lock_portable.yaml` | Exact pins without build hashes | No — regenerate with `conda env export --no-builds` |

**To add a new dependency:**

> **Do NOT** run `conda env export > environment.yaml` — this overwrites the minimal spec with a 240+ line platform-specific dump.

1. Add it to `environment.yaml` (conda section or pip section)
2. Recreate your environment: `conda env remove -n snakemake8 && conda env create -f environment.yaml`
3. Optionally regenerate lock files:
   ```bash
   conda env export | grep -v "^prefix:" > environment_lock_linux.yaml
   conda env export --no-builds | grep -v "^prefix:" > environment_lock_portable.yaml
   ```

**Per-rule conda environments:** If you add a tool with dependency conflicts that cannot be resolved in the main environment, Snakemake supports isolated per-rule conda environments via the `conda:` directive. Place a `.yaml` file in `workflows/envs/` and reference it in the rule:

```python
rule my_rule:
    ...
    conda:
        "../envs/my_tool.yaml"
    script:
        "../scripts/my_script.py"
```

> **Note for RHEL8/HPC users:** Snakemake's `conda:` directive activates isolated environments in a subprocess that does not inherit `LD_LIBRARY_PATH`. On systems where the system `libstdc++` is older than GLIBCXX_3.4.29 (e.g., RHEL8), this causes import failures for any modern conda-forge package. In that case, add the tool to `environment.yaml` instead.

## Key Concepts

**Preprocessing Step**: A preprocessing stage in single-cell analysis, such as demultiplexing or doublet detection. Each step has its own Snakemake rule file (`workflows/rules/<step>.smk`) and Pydantic schema (`schemas/<step>.py`) to encourage expansion to different tool options that solve the same preprocessing step e.g. [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite) + [vireo](https://github.com/single-cell-genetics/vireo) or [demuxalot](https://github.com/arogozhnikov/demuxalot).

**Method**: A software tool that implements a **preprocessing step**. For example, the demultiplexing step supports two methods: [vireo](https://github.com/single-cell-genetics/vireo) and [demuxalot](https://github.com/arogozhnikov/demuxalot). Adding a new method means integrating another tool into an existing prepcrocessing step.

## Developer Workflow

Follow these steps when making changes:

1. **Orient** - Understand the [project architecture](#project-architecture) and how the pipeline resolves what to run
2. **Develop** - Choose your task: [add a method](#adding-a-method-to-an-existing-step) or [add a step](#adding-a-new-pipeline-step)
3. **Test** - Validate with [DAG checks](#dag-validation) and [integration tests](#integration-tests)
4. **Document** - Update [documentation](#building-and-editing-documentation)!

---

## Project Architecture

Understanding the file layout is essential before making changes:

```
sc_preprocess/   
├── cli.py                          # CLI entry points
├── config_generator.py             # Interactive config builder
├── config_validator.py             # PIPELINE_DIRECTORIES, validation
├── schemas/                        # Pydantic models for config validation
│   ├── base.py                     # BaseStepConfig (all steps inherit this)
│   ├── config.py                   # PipelineConfig (unified schema)
│   ├── cellranger.py               # GEX/ATAC/ARC Cell Ranger configs
│   ├── demultiplexing.py           # DemuxalotConfig, VireoConfig
│   └── doublet_detection.py        # ScrubletConfig
├── workflows/
│   ├── main.smk                    # Master workflow, rule all, includes
│   ├── rules/                      # One .smk file per pipeline step
│   │   ├── cellranger.smk          # Cell Ranger count/aggregation
│   │   ├── object_creation.smk     # Per-capture AnnData/MuData creation
│   │   ├── batch_aggregation.smk   # Batch-level aggregation + metadata enrichment
│   │   ├── demultiplexing.smk      # Demuxalot/Vireo
│   │   └── doublet_detection.smk   # Scrublet
│   └── scripts/
│       ├── build_targets.py        # Generates target files for rule all
│       ├── parse_config.py         # Extracts enabled steps, methods, etc.
│       ├── create_gex_anndata.py   # GEX per-capture object creation
│       ├── create_atac_anndata.py  # ATAC per-capture object creation
│       ├── create_arc_mudata.py    # ARC per-capture MuData creation
│       ├── aggregate_batch.py      # Batch aggregation (per-capture → batch)
│       ├── merge_metadata.py       # Merge analysis metadata into batch objects
│       └── run_scrublet.py         # Scrublet doublet detection
tests/
├── test.sh                         # Integration test script
├── 00_TEST_DATA_GEX/               # Test configs and library lists
└── ...
docs/source/                        # Read the Docs documentation (this file)
```

### How the pipeline resolves what to run

1. **Config** (`pipeline_config.yaml`) declares which steps are `enabled: true`
2. `parse_config.py` → `get_enabled_steps()` reads the config and returns a list of enabled step names
3. `main.smk` conditionally includes `.smk` rule files based on enabled steps
4. `build_targets.py` → `build_all_targets()` generates the list of expected `.done` files for `rule all`
5. **Each rule** produces a `.done` marker file in `{output_dir}/00_LOGS/` that matches what `build_targets.py` expects

If the target filename from `build_targets.py` doesn't match the `done` output in the rule, Snakemake will raise a `MissingInputException`.

### Pipeline phases

The pipeline executes in phases, each building on the previous:

```
Phase 1: Cell Ranger count (per-capture)
  → 01_CELLRANGER{GEX|ATAC|ARC}_COUNT/{batch}_{capture}/

Phase 2: Cell Ranger aggregation (per-batch)
  → 02_CELLRANGER{GEX|ATAC|ARC}_AGGR/{batch}/

Phase 3: Per-capture object creation (AnnData/MuData with traceability metadata)
  → 03_ANNDATA/{batch}_{capture}.h5ad|h5mu

Phase 4: Batch aggregation (merge per-capture objects)
  → 04_BATCH_OBJECTS/{batch}_{modality}.h5ad|h5mu

Phase 5-6: Per-capture analysis (demux, doublet — run in parallel)
  → 05_DEMULTIPLEXING/, 06_DOUBLET_DETECTION/

Phase 7: Metadata enrichment (merge analysis results into batch objects)
  → 07_FINAL/{batch}_{modality}.h5ad|h5mu
```

**Metadata enrichment** is the final phase. The rules live in `batch_aggregation.smk` (alongside the aggregation rules) and the logic is in `scripts/merge_metadata.py`. For each batch, enrichment:

1. Reads the batch object from `04_BATCH_OBJECTS/`
2. Searches `05_DEMULTIPLEXING/` and `06_DOUBLET_DETECTION/` for per-capture result files
3. Joins analysis metadata onto `adata.obs` using `cell_id` as the key
4. Writes the enriched object to `07_FINAL/`

Enrichment is always the last step. It runs for every enabled modality regardless of which analysis steps are enabled — if no analysis metadata is found, it copies the batch object as-is.

Target generation for enrichment is handled by `get_enriched_object_outputs()` in `build_targets.py`, producing done files like `{batch}_gex_enrichment.done`.

---

## Development Workflow

### Quick Reference Checklists

**Adding a method to an existing preprocessing step** (e.g., a new demultiplexing tool):

1. Add method config schema in `schemas/<step>.py` with `tool_meta`
2. Register the method in the parent config class
3. Add the rule in `workflows/rules/<step>.smk` — **include a `threads:` and `resources:` block** (see [Resource parameters](#resource-parameters))
4. Add target generation in `workflows/scripts/build_targets.py`
5. Add to config generator in `config_generator.py`
6. [Test](#testing)
7. [Document](#building-and-editing-documentation)

**Adding a new pipeline step** (e.g., a new QC filter or analysis method):

1. Create the Pydantic schema in `schemas/<new_step>.py`
2. Register the output directory in `config_validator.py`
3. Register the step in `parse_config.py`
4. Add target generation in `build_targets.py`
5. Create the rule file in `workflows/rules/<new_step>.smk` — **include a `threads:` and `resources:` block** (see [Resource parameters](#resource-parameters))
6. Include the rule file in `main.smk`
7. Create a dummy rule and test the DAG
8. Implement the rule
9. Add to config generator
10. [Write tests](#testing)
11. [Document](#building-and-editing-documentation)

---

### Adding a Method to an Existing Preprocessing Step

This example shows how [Vireo](https://github.com/single-cell-genetics/vireo) was added alongside [demuxalot](https://github.com/arogozhnikov/demuxalot) for demultiplexing. Use this as a template for adding new methods.

#### Step 1: Add method config schema

In `sc_preprocess/schemas/demultiplexing.py`:

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

#### Checkpoint: Verify schema registration

After completing Steps 1-2, verify the new method is visible to the CLI:

```bash
# Confirm the method appears in the registry
sc-preprocess list-methods

# Confirm schema fields and tool_meta are correct
sc-preprocess show-params --step demultiplexing --method vireo
```

#### Step 3: Add the rule

In `sc_preprocess/workflows/rules/demultiplexing.smk`:

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
- **Always include `threads:` and `resources:` blocks** — without them SLURM defaults to 1 GB and jobs will OOM on real data:

```python
from tempfile import gettempdir  # must be imported explicitly in each .smk file

rule my_rule:
    ...
    threads: config["my_step"].get("threads", 1)
    resources:
        mem_mb = config["my_step"].get("mem_gb", 16) * 1024,  # SLURM uses MB
        runtime = config["my_step"].get("runtime_minutes", 720),  # passed as --time to SLURM
        tmpdir = RESOURCES.get("tmpdir") or gettempdir()
```

Add `threads`, `mem_gb`, and `runtime_minutes` to the **parent config class** (not the method sub-config) in `schemas/<step>.py`:

```python
threads: int = Field(default=1, ge=1, description="Number of threads")
mem_gb: int = Field(default=16, ge=1, description="Memory in GB")
runtime_minutes: int = Field(default=720, gt=0, description="Maximum runtime in minutes for the SLURM job")
```
- Method-specific config goes in the method-level `if` block

#### Step 4: Add target generation

In `sc_preprocess/workflows/scripts/build_targets.py`:

```python
if method == "vireo":
    for batch in batches:
        for capture in captures:
            outputs.append(os.path.join(logs_dir, f"vireo_output_{batch}_{capture}.done"))
```

#### Checkpoint: Verify rule and targets

After completing Steps 3-4, create a test config that uses the new method and verify the DAG resolves:

```bash
# Validate the config parses correctly
sc-preprocess validate-config --config-file your_test_config.yaml

# Dry run - confirm the new rule appears and targets match
sc-preprocess run --config-file your_test_config.yaml --cores 1 --dry-run

# Visual check - confirm rule dependencies look correct
sc-preprocess run --config-file your_test_config.yaml --cores 1 --dag | dot -Tpng > dag.png
```

If you get a `MissingInputException`, the `.done` filename in `build_targets.py` doesn't match the rule output. See [Common Pitfalls](#common-pitfalls).

#### Step 5: Add to config generator

Update `sc_preprocess/config_generator.py` so `init-config` can produce the new method's parameters interactively.

#### Step 6: Test

See [Testing](#testing) below.

---

### Adding a New Pipeline Step

Follow these steps to add a new preprocessing step (e.g., a quality-control filter, RNA velocity, etc.):

#### Step 1: Create the Pydantic schema

Create `sc_preprocess/schemas/<new_step>.py`:

```python
"""<Step name> configuration schemas."""

from typing import ClassVar, Literal, Optional
from pydantic import BaseModel, Field, model_validator
from .base import BaseStepConfig, ToolMeta


class MyMethodConfig(BaseModel):
    """Parameters for my method."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="my-package",
        url="https://github.com/my/package",
    )

    my_param: str = Field(description="An example parameter")

    class Config:
        extra = "forbid"


class MyStepConfig(BaseStepConfig):
    """My new step configuration."""

    method: Literal["my_method"] = Field(description="Method to use")
    my_method: Optional[MyMethodConfig] = None

    @model_validator(mode='after')
    def validate_method_params(self):
        # Validate that the selected method has its config block
        ...
```

#### Checkpoint: Verify schema registration

After Step 1, verify the new step and its methods are visible:

```bash
sc-preprocess list-methods
sc-preprocess show-params --step my_step --method my_method
```

#### Step 2: Register the output directory

In `sc_preprocess/config_validator.py`, add to `PIPELINE_DIRECTORIES`:

```python
("my_step", "0N_MY_STEP"),
```

#### Step 3: Register the step in parse_config.py

Add the step name to the `get_enabled_steps` list in `schemas/config.py`.

#### Step 4: Add target generation in build_targets.py

Add a call in `build_all_targets()` and a `get_<step>_outputs()` function following the pattern of `get_doublet_outputs()`.

#### Checkpoint: Verify config validation

```bash
sc-preprocess validate-config --config-file your_test_config.yaml
```

#### Step 5: Create the rule file

Create `sc_preprocess/workflows/rules/<new_step>.smk`. Always include `threads:` and `resources:` — without them SLURM defaults to 1 GB and jobs will OOM on real data. Import `gettempdir` explicitly:

```python
from tempfile import gettempdir
```

#### Step 6: Include the rule file in main.smk

```python
if "my_step" in ENABLED_STEPS:
    include: "rules/my_step.smk"
```

#### Step 7: Create a dummy rule and test the DAG

Before implementing the actual tool logic, write a dummy `shell` block:

```python
shell:
    """
    echo "Placeholder for my_method"
    touch {output.predictions}
    """
```

Then verify the DAG resolves correctly:

```bash
sc-preprocess run --config-file your_test_config.yaml --cores 1 --dry-run
sc-preprocess run --config-file your_test_config.yaml --cores 1 --dag | dot -Tpng > dag.png
```

Check that:
- Your new rule appears in the DAG
- It depends on the correct upstream rules (e.g., `cellranger_gex_count`)
- `rule all` connects to your new rule's `.done` output
- No `MissingInputException` errors

#### Step 8: Implement the rule

Replace the dummy shell with the actual tool invocation. Use either `shell:` for command-line tools or `run:` for Python-based tools.

#### Step 9: Add to config generator

Update `sc_preprocess/config_generator.py` so that `sc-preprocess init-config` can interactively generate config for the new step.

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
sc-preprocess run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dry-run

# Visual DAG - confirm rule dependencies look correct
sc-preprocess run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1 --dag | dot -Tpng > dag.png
```

### Integration tests

Integration tests run the full pipeline on test data:

```bash
# Run integration tests
bash tests/test.sh

# Or run a specific workflow manually
sc-preprocess run --config-file tests/00_TEST_DATA_GEX/test_config_gex.yaml --cores 1
```

### Test checklist for a new step

When adding a new step, verify all of the following before merging:

- [ ] **Config validation**: `sc-preprocess validate-config --config-file your_config.yaml` succeeds
- [ ] **Dry run**: `sc-preprocess run --config-file ... --cores 1 --dry-run` shows your rule
- [ ] **DAG**: Your rule appears with correct dependencies in the DAG visualization
- [ ] **Dummy execution**: Pipeline completes with placeholder `touch` commands
- [ ] **Real execution**: Pipeline completes with the actual tool on test data
- [ ] **Config generator**: `sc-preprocess init-config` includes the new step

### End-to-end test walkthrough

For a step-by-step walkthrough of running the full pipeline on test data (GEX, ATAC, and ARC), see:

```{toctree}
:maxdepth: 1

tutorial
```

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

This project uses [Sphinx](https://www.sphinx-doc.org/) with MyST-Parser for markdown support and is deployed on [Read the Docs](https://about.readthedocs.com/). Documentation is hosted at: https://sc-preprocess.readthedocs.io/

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

Here is how you can render the documentation locally with live reloading in your web browser:

```bash
sc-preprocess render-docs
```

If you would like to do it manually, here you go: 

```bash
cd docs

# install html rendering software
pip install sphinx
pip install sphinx-copybutton
pip install myst-parser
pip install sphinx-autobuild

# Build HTML
make html

# Serve locally
python3 -m http.server 8000 -d build/html
```

Then open http://localhost:8000 in your browser.

If you are using VS Code on a remote session, you can render the docs in the IDE itself.

For automatic rebuilds when you save changes, use `sphinx-autobuild`:

```bash
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
