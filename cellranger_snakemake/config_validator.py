"""Configuration validation and parameter documentation utilities."""

import sys
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, Type
from pydantic import ValidationError

from cellranger_snakemake.schemas.config import PipelineConfig
from cellranger_snakemake.schemas.cellranger import (
    CellRangerGEXConfig, CellRangerATACConfig, CellRangerARCConfig
)
from cellranger_snakemake.schemas.demultiplexing import (
    DemuxletConfig, SouporcellConfig, FreemuxletConfig, ViralConfig
)
from cellranger_snakemake.schemas.doublet_detection import (
    ScrubletConfig, DoubletFinderConfig, ScdsConfig, ScDblFinderConfig
)
from cellranger_snakemake.schemas.annotation import (
    CelltypistConfig, AzimuthConfig, SingleRConfig, ScTypeConfig, CelltypistCustomConfig
)


class ConfigValidator:
    """Configuration validation and help utilities."""
    
    # Method-specific parameter schemas
    METHOD_SCHEMAS = {
        "cellranger": {
            "gex": CellRangerGEXConfig,
            "atac": CellRangerATACConfig,
            "arc": CellRangerARCConfig,
        },
        "demultiplexing": {
            "demuxlet": DemuxletConfig,
            "souporcell": SouporcellConfig,
            "freemuxlet": FreemuxletConfig,
            "vireo": ViralConfig,
        },
        "doublet_detection": {
            "scrublet": ScrubletConfig,
            "doubletfinder": DoubletFinderConfig,
            "scds": ScdsConfig,
            "scdblfinder": ScDblFinderConfig,
        },
        "celltype_annotation": {
            "celltypist": CelltypistConfig,
            "azimuth": AzimuthConfig,
            "singler": SingleRConfig,
            "sctype": ScTypeConfig,
            "celltypist_custom": CelltypistCustomConfig,
        },
    }
    
    @classmethod
    def validate_config_file(cls, config_path: str) -> PipelineConfig:
        """
        Validate a configuration file.
        
        Note: The 'run' command automatically validates configs before execution.
        This method is useful for:
        - Quick validation during config development
        - CI/CD pipelines that need standalone validation
        - Checking configs without triggering Snakemake setup
        
        Args:
            config_path: Path to YAML configuration file
            
        Returns:
            Validated PipelineConfig object
            
        Raises:
            ValidationError: If configuration is invalid
        """
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        try:
            config = PipelineConfig(**config_dict)
            return config
        except ValidationError as e:
            print(f"\n❌ Configuration validation failed for: {config_path}\n")
            cls._print_validation_errors(e)
            raise
    
    @classmethod
    def _print_validation_errors(cls, error: ValidationError):
        """Print detailed validation errors with helpful messages."""
        for err in error.errors():
            location = " -> ".join(str(loc) for loc in err['loc'])
            msg = err['msg']
            err_type = err['type']
            
            print(f"\n  Location: {location}")
            print(f"  Error: {msg}")
            
            # Provide helpful context for common errors
            if 'required' in err_type.lower():
                print(f"  💡 This field is required. Please provide a value.")
            elif 'literal' in err_type.lower():
                # Extract allowed values from error context
                if 'expected' in err:
                    print(f"  💡 Allowed values: {err['expected']}")
            elif 'extra' in err_type.lower():
                print(f"  💡 This parameter is not recognized. Check for typos or see available parameters.")
    
    @classmethod
    def show_method_params(cls, step: str, method: str):
        """
        Show available parameters for a specific method.
        
        Args:
            step: Pipeline step (e.g., 'cellranger', 'demultiplexing', 'doublet_detection', 'celltype_annotation')
            method: Method name (e.g., 'gex', 'atac', 'arc', 'scrublet', 'demuxlet', etc.)
        """
        if step not in cls.METHOD_SCHEMAS:
            print(f"❌ Unknown step: {step}")
            print(f"Available steps: {', '.join(cls.METHOD_SCHEMAS.keys())}")
            return
        
        if method not in cls.METHOD_SCHEMAS[step]:
            print(f"❌ Unknown method '{method}' for step '{step}'")
            print(f"Available methods: {', '.join(cls.METHOD_SCHEMAS[step].keys())}")
            return
        
        schema_class = cls.METHOD_SCHEMAS[step][method]
        cls._print_schema_params(schema_class, f"{step}.{method}")
    
    @classmethod
    def _print_schema_params(cls, schema_class: Type, title: str):
        """Print parameters from a pydantic schema."""
        print(f"\n{'='*60}")
        print(f"Parameters for: {title}")
        print(f"{'='*60}\n")
        
        schema = schema_class.model_json_schema()
        properties = schema.get('properties', {})
        required = schema.get('required', [])
        
        for param_name, param_info in properties.items():
            is_required = param_name in required
            param_type = param_info.get('type', 'unknown')
            default = param_info.get('default', 'N/A')
            description = param_info.get('description', '')
            
            # Handle enum/literal types
            if 'enum' in param_info:
                param_type = f"Literal[{', '.join(repr(v) for v in param_info['enum'])}]"
            elif 'anyOf' in param_info:
                types = [t.get('type', 'unknown') for t in param_info['anyOf']]
                param_type = ' | '.join(types)
            
            # Format parameter info
            required_str = "REQUIRED" if is_required else f"optional (default: {default})"
            
            print(f"  {param_name}")
            print(f"    Type: {param_type}")
            print(f"    Status: {required_str}")
            if description:
                print(f"    Description: {description}")
            
            # Show constraints if any
            if 'minimum' in param_info or 'maximum' in param_info:
                constraints = []
                if 'minimum' in param_info:
                    constraints.append(f">= {param_info['minimum']}")
                if 'maximum' in param_info:
                    constraints.append(f"<= {param_info['maximum']}")
                print(f"    Constraints: {', '.join(constraints)}")
            
            print()
    
    @classmethod
    def list_all_methods(cls):
        """List all available methods for each step."""
        print("\n" + "="*60)
        print("Available Methods by Pipeline Step")
        print("="*60 + "\n")
        
        for step, methods in cls.METHOD_SCHEMAS.items():
            print(f"\n{step.replace('_', ' ').title()}:")
            for method in methods.keys():
                print(f"  • {method}")
    
    @classmethod
    def generate_example_config(cls, output_path: Optional[str] = None) -> str:
        """
        Generate an example configuration file with all options.
        
        Args:
            output_path: Optional path to save the example config
            
        Returns:
            Path to generated config file
        """
        example_config = {
            "project_name": "my_scrna_project",
            "output_dir": "output",
            "samples": {
                "sample1": {
                    "sample_id": "sample1",
                    "batch": "batch1",
                    "gex_fastqs": "/path/to/gex/fastqs",
                    "expected_cells": 5000,
                }
            },
            "hpc": {
                "mode": "local",
                "mempercore": None
            },
            "resources": {
                "mem_gb": 32,
                "tmpdir": "/tmp",
            },
            "directories_suffix": "none",
            "cellranger_gex": {
                "enabled": True,
                "reference": "/path/to/reference",
                "libraries": "/path/to/libraries.tsv",
                "chemistry": "auto",
                "normalize": "none",
            },
            "doublet_detection": {
                "enabled": True,
                "method": "scrublet",
                "scrublet": {
                    "expected_doublet_rate": 0.06,
                    "min_counts": 2,
                    "min_cells": 3,
                },
            },
            "celltype_annotation": {
                "enabled": True,
                "method": "celltypist",
                "celltypist": {
                    "model": "Immune_All_Low.pkl",
                    "majority_voting": False,
                },
            },
        }
        
        if output_path is None:
            output_path = "example_pipeline_config.yaml"
        
        with open(output_path, 'w') as f:
            yaml.dump(example_config, f, default_flow_style=False, sort_keys=False, indent=2)
        
        print(f"\n✓ Example configuration saved to: {output_path}")
        return output_path


def main():
    """CLI for config validation and help."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Configuration validation and parameter documentation"
    )
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate a configuration file')
    validate_parser.add_argument('config', help='Path to configuration file')
    
    # Show params command
    params_parser = subparsers.add_parser('show-params', help='Show parameters for a method')
    params_parser.add_argument('--step', required=True, 
                               help='Pipeline step (cellranger, demultiplexing, doublet_detection, celltype_annotation)')
    params_parser.add_argument('--method', required=True, help='Method name')
    
    # List methods command
    subparsers.add_parser('list-methods', help='List all available methods')
    
    # Generate example command
    example_parser = subparsers.add_parser('generate-example', help='Generate example configuration')
    example_parser.add_argument('--output', default='example_pipeline_config.yaml',
                                help='Output path for example config')
    
    args = parser.parse_args()
    
    if args.command == 'validate':
        try:
            config = ConfigValidator.validate_config_file(args.config)
            print(f"\n✓ Configuration is valid!")
            print(f"\nEnabled steps: {', '.join(config.get_enabled_steps())}")
        except Exception:
            sys.exit(1)
    
    elif args.command == 'show-params':
        ConfigValidator.show_method_params(args.step, args.method)
    
    elif args.command == 'list-methods':
        ConfigValidator.list_all_methods()
    
    elif args.command == 'generate-example':
        ConfigValidator.generate_example_config(args.output)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
