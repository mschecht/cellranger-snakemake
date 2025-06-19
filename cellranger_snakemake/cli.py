#!/usr/bin/env python

import sys
import argparse

from cellranger_snakemake.utils.custom_logger import custom_logger
from cellranger_snakemake.workflows.workflow import Workflow
from cellranger_snakemake.config_templates import ARC_CONFIG, ATAC_CONFIG, GEX_CONFIG

__version__ = "1.0.0"
__description__ = "Snakemake wrapper for Cell Ranger ARC workflows"

WORKFLOW_CONFIGS = {
    "ARC": ARC_CONFIG,
    "ATAC": ATAC_CONFIG,
    "GEX": GEX_CONFIG,
}

def main():
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Print default config to stdout
  snakemake-run-cellranger --workflow ARC --get-default-config

  # Run workflow (dry run)
  snakemake-run-cellranger --workflow ARC --dry-run

  # Run workflow with 4 cores
  snakemake-run-cellranger --workflow ARC --run --cores 4

  # Run with custom config file
  snakemake-run-cellranger --workflow ARC --run --config-file my_config.yaml
        """
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    parser.add_argument(
        "--workflow",
        required=True,
        choices=["ARC", "ATAC", "GEX"],
        help="Select which Cell Ranger workflow to use."
    )

    action_group = parser.add_mutually_exclusive_group(required=False)  # no longer required

    action_group.add_argument(
        "--get-default-config",
        dest="get_default_config",
        action="store",
        nargs="?",
        const=True,
        default=False,
        metavar="FILE",
        help="Write default configuration YAML to FILE (or NAME_default_config.yaml if none provided)"
    )
    action_group.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run the workflow"
    )
    action_group.add_argument(
        "--dag-to-pdf",
        action="store_true",
        help="Export PDF of workflow DAG"
    )
    parser.add_argument("--config-file", help="Config file (default: config.yaml)")
    parser.add_argument("--snakefile", default="Snakefile", help="Snakefile to use (default: Snakefile)")
    parser.add_argument("--cores", type=int, help="Number of CPU cores to use")
    parser.add_argument("--output-dir", default=".", help="Output dir for default files")
    parser.add_argument("--additional-params", default="", help="Extra Snakemake args")

    args = parser.parse_args()

    # If not using --get-default-config, we need a path to a config file
    if not args.get_default_config:
        if not args.config_file:
            parser.error("`--config-file` is required unless using `--get-default-config`")

    config = WORKFLOW_CONFIGS[args.workflow]
    wf = Workflow(name=args.workflow, get_default_config=args.get_default_config, config_file=args.config_file)

    try:
        if args.get_default_config:
            wf.write_default_config()
            wf.generate_config_readme()
        else:
            wf.run(
                config_file=args.config_file,
                dry_run=args.dry_run,
                dag=args.dag_to_pdf,
                cores=args.cores,
                additional_args=args.additional_params
            )
    except Exception as e:
        custom_logger.error(f"Failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
