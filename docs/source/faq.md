# FAQ

## How do `.done` files track finished steps?

Each pipeline step writes an empty flag file to `00_LOGS/` when it completes successfully — for example, `1_L001_cellranger_gex.done`. Snakemake uses these files as targets: if the `.done` file already exists, the step is skipped on re-run.

To force a step to re-run, delete its `.done` file, or use `--forcerun` (see below).

## How do I re-run the workflow from a specific step?

Pass `--forcerun <rule_name>` to Snakemake via `--snakemake-args`. For example, to restart from `create_gex_anndata`:

```bash
sc-preprocess run --config-file pipeline_config.yaml \
                             --cores 1 \
                             --snakemake-args --forcerun create_gex_anndata
```

To find available rule names, visualize the pipeline DAG:

```bash
sc-preprocess run --config-file pipeline_config.yaml --cores 1 --dag | dot -Tpng > dag.png
```

## What happens when a Cell Ranger cluster-mode job is killed?

When `jobmode` is set to `slurm` (or another scheduler), Cell Ranger's Martian runtime saves checkpoints after each completed stage. If the job is killed mid-run, Cell Ranger leaves a `_lock` file in the output directory. You read more about Cell Ranger cluster mode [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode).

The pipeline automatically removes this `_lock` file when resubmitted, so the job resumes from the last completed stage rather than starting over.
