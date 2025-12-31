"""Cell Ranger version detection and validation utilities."""

import subprocess
import re
from typing import Optional, Dict, Tuple
from packaging import version
from cellranger_snakemake.utils.custom_logger import custom_logger


class CellRangerVersionChecker:
    """Detect and validate Cell Ranger tool versions."""
    
    # Minimum supported versions
    MIN_VERSIONS = {
        "cellranger": "9.0.0",
        "cellranger-atac": "2.0.0",
        "cellranger-arc": "2.0.0",
    }
    
    # Version-specific test data paths (relative to cellranger installation directory)
    # These patterns are tried in order until a match is found
    TEST_DATA_PATHS = {
        'GEX': {
            'fastq_patterns': [
                'external/cellranger_tiny_fastq',
            ],
            'ref_patterns': [
                'external/cellranger_tiny_ref',
            ]
        },
        'ATAC': {
            'fastq_patterns': [
                'external/cellranger_atac_tiny_fastq/1.0.0',
                'external/cellranger_atac_tiny_fastq',
            ],
            'ref_patterns': [
                'external/arc_testrun_files/reference/fasta',
                'external/arc_testrun_files/reference',
                'external/cellranger_atac_tiny_ref',
            ]
        },
        'ARC': {
            'fastq_patterns': [
                'external/arc_testrun_files/fastqs',
                'external/cellranger_arc_tiny_fastq',
            ],
            'ref_patterns': [
                'external/arc_testrun_files/reference',
                'external/cellranger_arc_tiny_ref',
            ]
        }
    }
    
    @staticmethod
    def _run_version_command(command: str) -> Optional[str]:
        """
        Run a command with --version flag and capture output.
        
        Args:
            command: Command to run (e.g., "cellranger", "cellranger-atac")
            
        Returns:
            Version string if successful, None if command not found
        """
        try:
            result = subprocess.run(
                [command, "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                # Cell Ranger outputs like: "cellranger cellranger-7.0.1"
                # or "cellranger-atac cellranger-atac-2.1.0"
                output = result.stdout.strip() or result.stderr.strip()
                
                # Extract version number (e.g., "7.0.1" or "2.1.0")
                # Pattern: tool-name followed by semantic version
                pattern = r'(\d+\.\d+\.\d+)'
                match = re.search(pattern, output)
                
                if match:
                    return match.group(1)
                else:
                    custom_logger.warning(
                        f"Could not parse version from '{command} --version' output: {output}"
                    )
                    return None
            else:
                return None
                
        except FileNotFoundError:
            # Command not found
            return None
        except subprocess.TimeoutExpired:
            custom_logger.warning(f"Timeout running '{command} --version'")
            return None
        except Exception as e:
            custom_logger.warning(f"Error checking {command} version: {e}")
            return None
    
    @classmethod
    def get_version(cls, tool: str) -> Optional[str]:
        """
        Get version of a Cell Ranger tool.
        
        Args:
            tool: Tool name ("cellranger", "cellranger-atac", or "cellranger-arc")
            
        Returns:
            Version string (e.g., "7.0.1") or None if not found
        """
        return cls._run_version_command(tool)
    
    @classmethod
    def get_all_versions(cls) -> Dict[str, Optional[str]]:
        """
        Get versions of all Cell Ranger tools.
        
        Returns:
            Dictionary mapping tool names to version strings (or None if not found)
        """
        return {
            "cellranger": cls.get_version("cellranger"),
            "cellranger-atac": cls.get_version("cellranger-atac"),
            "cellranger-arc": cls.get_version("cellranger-arc"),
        }
    
    @classmethod
    def check_version(cls, tool: str, required_version: Optional[str] = None) -> Tuple[bool, str]:
        """
        Check if a tool meets minimum version requirement.
        
        Args:
            tool: Tool name
            required_version: Minimum version required (defaults to MIN_VERSIONS)
            
        Returns:
            Tuple of (is_valid, message)
        """
        detected_version = cls.get_version(tool)
        
        if detected_version is None:
            return False, f"'{tool}' not found in PATH. Please install it."
        
        min_version = required_version or cls.MIN_VERSIONS.get(tool)
        
        if min_version is None:
            # No minimum version requirement
            return True, f"'{tool}' version {detected_version} detected"
        
        try:
            if version.parse(detected_version) >= version.parse(min_version):
                return True, f"'{tool}' version {detected_version} meets minimum requirement ({min_version})"
            else:
                return False, f"'{tool}' version {detected_version} is below minimum requirement ({min_version})"
        except Exception as e:
            return False, f"Error comparing versions: {e}"
    
    @classmethod
    def validate_for_workflow(cls, workflow: str) -> Tuple[bool, list]:
        """
        Validate that required Cell Ranger tools are installed for a workflow.
        
        Args:
            workflow: Workflow name ("GEX", "ATAC", or "ARC")
            
        Returns:
            Tuple of (all_valid, messages_list)
        """
        workflow_tools = {
            "GEX": ["cellranger"],
            "ATAC": ["cellranger-atac"],
            "ARC": ["cellranger-arc"],
        }
        
        required_tools = workflow_tools.get(workflow.upper(), [])
        
        if not required_tools:
            return False, [f"Unknown workflow: {workflow}"]
        
        messages = []
        all_valid = True
        
        for tool in required_tools:
            is_valid, message = cls.check_version(tool)
            messages.append(message)
            
            if not is_valid:
                all_valid = False
        
        return all_valid, messages
    
    @classmethod
    def print_all_versions(cls):
        """Print detected versions of all Cell Ranger tools."""
        versions = cls.get_all_versions()
        
        custom_logger.info("Detected Cell Ranger versions:")
        for tool, ver in versions.items():
            if ver:
                min_ver = cls.MIN_VERSIONS.get(tool, "N/A")
                status = "✓" if ver and version.parse(ver) >= version.parse(min_ver) else "✗"
                custom_logger.info(f"  {status} {tool}: {ver} (min: {min_ver})")
            else:
                custom_logger.warning(f"  ✗ {tool}: not found")


def validate_cellranger_versions(workflow: str = None, verbose: bool = True) -> bool:
    """
    Convenience function to validate Cell Ranger versions.
    
    Args:
        workflow: Workflow name ("GEX", "ATAC", "ARC"). If None, checks all tools.
        verbose: Whether to print validation results
        
    Returns:
        True if all required tools are valid, False otherwise
    """
    checker = CellRangerVersionChecker()
    
    if workflow:
        is_valid, messages = checker.validate_for_workflow(workflow)
        
        if verbose:
            for msg in messages:
                if "not found" in msg.lower() or "below minimum" in msg.lower():
                    custom_logger.error(msg)
                else:
                    custom_logger.info(msg)
        
        return is_valid
    else:
        # Check all tools
        if verbose:
            checker.print_all_versions()
        
        versions = checker.get_all_versions()
        return all(v is not None for v in versions.values())
