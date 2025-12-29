#!/usr/bin/env python

import sys
import os
import shlex
import argparse
import subprocess

from pathlib import Path

from cellranger_snakemake.utils.custom_logger import custom_logger
from cellranger_snakemake.config_generator import ConfigGenerator
from cellranger_snakemake.config_validator import ConfigValidator
from cellranger_snakemake.utils.version_check import CellRangerVersionChecker
from cellranger_snakemake.utils.test_data_generator import TestDataGenerator

__version__ = "2.0.0"
__description__ = "Snakemake wrapper for single-cell preprocessing pipelines"

def main():
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate unified config interactively
  snakemake-run-cellranger init-config

  # Quick config validation (optional - 'run' auto-validates)
  snakemake-run-cellranger validate-config --config-file pipeline_config.yaml

  # Check installed Cell Ranger versions
  snakemake-run-cellranger check-versions
  snakemake-run-cellranger check-versions --workflow GEX

  # Run the pipeline (dry-run, auto-validates config first)
  snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 1 --dry-run

  # Run the pipeline (real run with 8 cores, auto-validates config first)
  snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 8

  # Generate DAG visualization
  snakemake-run-cellranger run --config-file pipeline_config.yaml --cores 1 --dag | dot -Tpng > dag.png

  # Show available parameters for a Cell Ranger modality
  snakemake-run-cellranger show-params --step cellranger --method gex

  # Show available parameters for a processing method
  snakemake-run-cellranger show-params --step doublet_detection --method scrublet

  # List all available methods
  snakemake-run-cellranger list-methods

  # Generate example config
  snakemake-run-cellranger generate-test-data ATAC --output-dir 00_TEST_DATA
        """
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    # Subcommands
    subparsers = parser.add_subparsers(dest='subcommand', help='Available commands')
    
    # Run workflow subcommand
    run_parser = subparsers.add_parser(
        'run',
        help='Run the Snakemake workflow (auto-validates config)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Note: Config is automatically validated before running. No need to run validate-config first.

Additional Snakemake arguments can be passed in two ways:
1. Directly as arguments (e.g., --cores 8 --dry-run)
2. Via --additional-params for complex arguments (e.g., cluster config with quotes)

Common Snakemake arguments (--cores is REQUIRED):
  --cores N, -c N            Number of cores to use (REQUIRED - e.g., --cores 8)
  --cores all                Use all available cores
  --dry-run, -n              Perform a dry run (still requires --cores)
  --dag                      Generate DAG visualization
  --rulegraph                Generate rule graph
  --forceall, -F             Force re-run of all rules
  --unlock                   Unlock working directory
  --cluster "cmd"            Submit jobs to cluster
  --cluster-config FILE      Cluster configuration file
  --jobs N                   Use at most N CPU cluster/cloud jobs in parallel

Examples:
  # Dry run to check what will be executed
  snakemake-run-cellranger run --config-file config.yaml --cores 1 --dry-run
  
  # Run with 8 cores
  snakemake-run-cellranger run --config-file config.yaml --cores 8
  
  # Use all available cores
  snakemake-run-cellranger run --config-file config.yaml --cores all
  
  # Run with cluster execution
  snakemake-run-cellranger run --config-file config.yaml --cores 1 \\
    --additional-params "--cluster 'sbatch -J {rule} --ntasks=1 --cpus-per-task=12 --mem=40G' --jobs 10"
  
  # Unlock working directory
  snakemake-run-cellranger run --config-file config.yaml --cores 1 --additional-params "--unlock"
        """
    )
    run_parser.add_argument(
        '--config-file',
        required=True,
        help='Path to pipeline configuration file'
    )
    run_parser.add_argument(
        '--additional-params',
        default='',
        help='Additional parameters to pass to Snakemake (e.g., "--unlock" or cluster config)'
    )
    
    # Init config subcommand
    init_parser = subparsers.add_parser(
        'init-config',
        help='Generate configuration file interactively'
    )
    init_parser.add_argument(
        '--output',
        default='pipeline_config.yaml',
        help='Output path for configuration file'
    )
    
    # Validate config subcommand
    validate_parser = subparsers.add_parser(
        'validate-config',
        help='Validate a configuration file (quick check without running workflow)'
    )
    validate_parser.add_argument(
        '--config-file',
        required=True,
        help='Path to configuration file to validate'
    )
    
    # Show params subcommand
    params_parser = subparsers.add_parser(
        'show-params',
        help='Show available parameters for a method'
    )
    params_parser.add_argument(
        '--step',
        required=True,
        help='Pipeline step (cellranger, demultiplexing, doublet_detection, celltype_annotation)'
    )
    params_parser.add_argument(
        '--method',
        required=True,
        help='Method name'
    )
    
    # List methods subcommand
    subparsers.add_parser(
        'list-methods',
        help='List all available methods for each pipeline step'
    )
    
    # Check versions subcommand
    versions_parser = subparsers.add_parser(
        'check-versions',
        help='Check installed Cell Ranger versions'
    )
    versions_parser.add_argument(
        '--workflow',
        choices=['GEX', 'ATAC', 'ARC'],
        help='Check versions for specific workflow (default: check all)'
    )
    
    # Generate test data subcommand
    testdata_parser = subparsers.add_parser(
        'generate-test-data',
        help='Generate test data for a specific workflow'
    )
    testdata_parser.add_argument(
        'workflow',
        choices=['GEX', 'ATAC', 'ARC'],
        help='Workflow type (GEX, ATAC, or ARC)'
    )
    testdata_parser.add_argument(
        '--output-dir',
        default='00_TEST_DATA',
        help='Directory where test data will be generated (default: 00_TEST_DATA)'
    )

    # Parse known args to allow passing unknown args to snakemake
    args, unknown_args = parser.parse_known_args()

    # Handle subcommands
    if args.subcommand == 'run':
        # Get the path to the main.smk file
        package_dir = Path(__file__).parent
        snakefile = os.path.join(package_dir, "workflows", "main.smk")
        
        if not os.path.exists(snakefile):
            custom_logger.error(f"Snakefile not found: {snakefile}")
            sys.exit(1)
        
        # Validate config file exists
        if not os.path.exists(args.config_file):
            custom_logger.error(f"Config file not found: {args.config_file}")
            sys.exit(1)
        
        # Validate config before running
        try:
            config = ConfigValidator.validate_config_file(args.config_file)
            custom_logger.info(f"Config validated. Enabled steps: {', '.join(config.get_enabled_steps())}")
        except Exception:
            custom_logger.error("Config validation failed. Fix errors before running.")
            sys.exit(1)
        
        # Build snakemake command
        snakemake_cmd = [
            "snakemake",
            "--snakefile", str(snakefile),
            "--configfile", args.config_file
        ]
        
        # Add additional parameters if provided
        if args.additional_params:
            snakemake_cmd.extend(shlex.split(args.additional_params))
        
        # Add any additional snakemake arguments from command line
        snakemake_cmd.extend(unknown_args)
        
        # Run snakemake
        try:
            result = subprocess.run(snakemake_cmd)
            sys.exit(result.returncode)
        except FileNotFoundError:
            custom_logger.error("Snakemake not found. Please install snakemake: pip install snakemake")
            sys.exit(1)
        except KeyboardInterrupt:
            custom_logger.warning("Pipeline execution cancelled.")
            sys.exit(130)
    
    elif args.subcommand == 'init-config':
        generator = ConfigGenerator()
        try:
            config = generator.generate()
            generator.save_config(config, args.output)
            custom_logger.info(f"Enabled steps: {', '.join(config.get_enabled_steps())}")
        except KeyboardInterrupt:
            custom_logger.warning("Configuration generation cancelled.")
            sys.exit(0)
        except Exception as e:
            custom_logger.error(f"Failed to generate config: {e}")
            sys.exit(1)
    
    elif args.subcommand == 'validate-config':
        try:
            config = ConfigValidator.validate_config_file(args.config_file)
            custom_logger.info("Configuration is valid!")
            custom_logger.info(f"Enabled steps: {', '.join(config.get_enabled_steps())}")
        except Exception:
            sys.exit(1)
    
    elif args.subcommand == 'show-params':
        ConfigValidator.show_method_params(args.step, args.method)
    
    elif args.subcommand == 'list-methods':
        ConfigValidator.list_all_methods()
    
    elif args.subcommand == 'generate-test-data':
        workflow_type = args.workflow  # 'GEX', 'ATAC', or 'ARC'
        output_dir = Path(args.output_dir)
        
        custom_logger.info(f"Generating {workflow_type} test data in {output_dir}")
        
        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate test data files (libraries TSV, reference file)
        test_data_paths = TestDataGenerator.generate_test_data(workflow_type, output_dir)
        
        if not test_data_paths:
            custom_logger.error("Failed to generate test data")
            sys.exit(1)
        
        # Generate test configuration file with actual paths
        config_path = output_dir / f"test_config_{workflow_type.lower()}.yaml"
        ConfigValidator.generate_test_data_config(
            workflow=workflow_type,
            test_data_dir=str(output_dir),
            reference_path=test_data_paths['reference_path'],
            output_path=str(config_path)
        )
        
        custom_logger.info(f"\n{'='*60}")
        custom_logger.info(f"Test data generation complete!")
        custom_logger.info(f"{'='*60}")
        custom_logger.info(f"Libraries file: {test_data_paths['libraries_file']}")
        custom_logger.info(f"Reference path: {test_data_paths['reference_path']}")
        custom_logger.info(f"Config file: {config_path}")
        custom_logger.info(f"\nTo run the {workflow_type} test pipeline:")
        custom_logger.info(f"  snakemake-run-cellranger run --config-file {config_path} --cores 8")
    
    elif args.subcommand == 'check-versions':
        checker = CellRangerVersionChecker()
        
        if args.workflow:
            # Check specific workflow
            is_valid, messages = checker.validate_for_workflow(args.workflow)
            
            for msg in messages:
                if "not found" in msg.lower() or "below minimum" in msg.lower():
                    custom_logger.error(msg)
                else:
                    custom_logger.info(msg)
            
            if is_valid:
                custom_logger.info(f"✓ All required tools for {args.workflow} workflow are installed and meet minimum version requirements")
                sys.exit(0)
            else:
                custom_logger.error(f"✗ Version check failed for {args.workflow} workflow")
                sys.exit(1)
        else:
            # Check all tools
            checker.print_all_versions()
            versions = checker.get_all_versions()
            
            if all(v is not None for v in versions.values()):
                custom_logger.info("\n✓ All Cell Ranger tools are installed")
                sys.exit(0)
            else:
                custom_logger.warning("\n⚠ Some Cell Ranger tools are not installed")
                sys.exit(0)  # Don't exit with error if just checking all
    
    else:
        parser.print_help()
        sys.exit(0)


if __name__ == "__main__":
    main()