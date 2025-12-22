"""Test data generation utilities for cellranger workflows."""

import os
import subprocess
from pathlib import Path
from typing import Optional, Tuple, Dict
from cellranger_snakemake.utils.custom_logger import custom_logger
from cellranger_snakemake.utils.version_check import CellRangerVersionChecker


class TestDataGenerator:
    """Generate test data files and configurations for cellranger workflows."""
    
    # Map workflow types to cellranger programs
    WORKFLOW_PROGRAMS = {
        'GEX': 'cellranger',
        'ATAC': 'cellranger-atac',
        'ARC': 'cellranger-arc'
    }
    
    @staticmethod
    def find_cellranger_path(program: str) -> Optional[Path]:
        """
        Find the installation path of a cellranger program.
        
        Args:
            program: Name of cellranger program (e.g., 'cellranger', 'cellranger-atac')
            
        Returns:
            Path to cellranger installation directory, or None if not found
        """
        try:
            # Use 'which' to find the program
            result = subprocess.run(
                ['which', program],
                capture_output=True,
                text=True,
                check=True
            )
            program_path = result.stdout.strip()
            
            # Resolve symlinks
            real_path = subprocess.run(
                ['readlink', '-f', program_path],
                capture_output=True,
                text=True,
                check=True
            )
            real_path_str = real_path.stdout.strip()
            
            # Get cellranger directory (two levels up from binary)
            cellranger_dir = Path(real_path_str).parent.parent
            return cellranger_dir
            
        except subprocess.CalledProcessError:
            custom_logger.warning(f"{program} not found in PATH")
            return None
    
    @classmethod
    def find_tiny_test_data(cls, workflow: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Find tiny test data paths for a specific workflow.
        
        Args:
            workflow: Workflow type ('GEX', 'ATAC', 'ARC')
            
        Returns:
            Tuple of (fastq_path, reference_path), or (None, None) if not found
        """
        program = cls.WORKFLOW_PROGRAMS.get(workflow)
        if not program:
            custom_logger.error(f"Unknown workflow type: {workflow}")
            return None, None
        
        cellranger_dir = cls.find_cellranger_path(program)
        if not cellranger_dir:
            return None, None
        
        # Get version-specific path patterns
        path_config = CellRangerVersionChecker.TEST_DATA_PATHS.get(workflow)
        if not path_config:
            custom_logger.error(f"No test data path configuration for workflow: {workflow}")
            return None, None
        
        fastq_patterns = path_config['fastq_patterns']
        ref_patterns = path_config['ref_patterns']
        
        # Try to find FASTQ directory
        fastq_path = None
        for pattern in fastq_patterns:
            candidate = cellranger_dir / pattern
            if candidate.exists():
                fastq_path = str(candidate)
                custom_logger.info(f"Found FASTQ path: {fastq_path}")
                break
        
        if not fastq_path:
            custom_logger.warning(f"FASTQ path not found. Tried patterns: {fastq_patterns}")
            return None, None
        
        # Try to find reference directory
        ref_path = None
        for pattern in ref_patterns:
            candidate = cellranger_dir / pattern
            if candidate.exists():
                ref_path = str(candidate)
                custom_logger.info(f"Found reference path: {ref_path}")
                break
        
        if not ref_path:
            custom_logger.warning(f"Reference path not found. Tried patterns: {ref_patterns}")
            return fastq_path, None
        
        return fastq_path, ref_path
    
    @classmethod
    def generate_libraries_tsv(cls, workflow: str, fastq_path: str, output_path: Path) -> Path:
        """
        Generate a libraries TSV file for test data.
        
        Args:
            workflow: Workflow type ('GEX', 'ATAC', 'ARC')
            fastq_path: Path to fastq directory
            output_path: Output file path
            
        Returns:
            Path to created TSV file
        """
        with open(output_path, 'w') as f:
            # Write header
            if workflow == 'ARC':
                f.write("batch\tcapture\tsample\tfastqs\tlibrary_type\n")
                # ARC requires both GEX and ATAC libraries
                f.write(f"1\tL001\ttiny{workflow.lower()}\t{fastq_path}\tGene Expression\n")
                f.write(f"1\tL001\ttiny{workflow.lower()}\t{fastq_path}\tChromatin Accessibility\n")
            else:
                f.write("batch\tcapture\tsample\tfastqs\n")
                f.write(f"1\tL001\ttiny{workflow.lower()}\t{fastq_path}\n")
                f.write(f"1\tL002\ttiny{workflow.lower()}\t{fastq_path}\n")
        
        custom_logger.info(f"Created libraries TSV: {output_path}")
        return output_path
    
    @classmethod
    def generate_test_data(cls, workflow: str, output_dir: Path) -> Dict[str, str]:
        """
        Generate complete test data setup for a workflow.
        
        Args:
            workflow: Workflow type ('GEX', 'ATAC', 'ARC')
            output_dir: Directory where test data files will be created
            
        Returns:
            Dictionary with paths to generated files
        """
        # Find test data paths
        fastq_path, ref_path = cls.find_tiny_test_data(workflow)
        
        if not fastq_path:
            custom_logger.error(f"Could not find tiny test data for {workflow}")
            custom_logger.info(f"Make sure {cls.WORKFLOW_PROGRAMS[workflow]} is installed and in your PATH")
            return {}
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate libraries TSV
        libraries_file = output_dir / f"libraries_list_{workflow.lower()}.tsv"
        cls.generate_libraries_tsv(workflow, fastq_path, libraries_file)
        
        # Create reference file (optional, for documentation)
        ref_file = output_dir / f"reference_{workflow.lower()}.txt"
        if ref_path:
            with open(ref_file, 'w') as f:
                f.write(ref_path)
            custom_logger.info(f"Created reference file: {ref_file}")
        else:
            custom_logger.warning(f"Reference path not found, skipping reference file creation")
            ref_path = f"/path/to/{workflow.lower()}/reference"
        
        return {
            'libraries_file': str(libraries_file),
            'reference_path': ref_path,
            'fastq_path': fastq_path
        }
