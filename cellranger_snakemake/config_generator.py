"""Interactive configuration generator for single-cell preprocessing pipeline."""

import sys
import yaml
from typing import Optional, Dict, Any
from pathlib import Path

try:
    from rich.console import Console
    from rich.prompt import Prompt, Confirm, IntPrompt, FloatPrompt
    from rich.panel import Panel
    from rich.markdown import Markdown
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False
    # Fallback to basic input
    class Console:
        def print(self, *args, **kwargs):
            print(*args)
    
    class Prompt:
        @staticmethod
        def ask(prompt, choices=None, default=None):
            if choices:
                prompt += f" ({'/'.join(choices)})"
            if default:
                prompt += f" [{default}]"
            response = input(f"{prompt}: ").strip()
            return response if response else default
    
    class Confirm:
        @staticmethod
        def ask(prompt, default=True):
            default_str = "Y/n" if default else "y/N"
            response = input(f"{prompt} ({default_str}): ").strip().lower()
            if not response:
                return default
            return response in ['y', 'yes']
    
    class IntPrompt:
        @staticmethod
        def ask(prompt, default=None):
            if default:
                prompt += f" [{default}]"
            while True:
                response = input(f"{prompt}: ").strip()
                if not response and default is not None:
                    return default
                try:
                    return int(response)
                except ValueError:
                    print("Please enter a valid integer.")
    
    class FloatPrompt:
        @staticmethod
        def ask(prompt, default=None):
            if default:
                prompt += f" [{default}]"
            while True:
                response = input(f"{prompt}: ").strip()
                if not response and default is not None:
                    return default
                try:
                    return float(response)
                except ValueError:
                    print("Please enter a valid number.")


from cellranger_snakemake.schemas.config import PipelineConfig, SampleMetadata
from cellranger_snakemake.schemas.base import ResourceConfig
from cellranger_snakemake.schemas.cellranger import (
    CellRangerGEXConfig, CellRangerATACConfig, CellRangerARCConfig
)
from cellranger_snakemake.schemas.demultiplexing import (
    DemultiplexingConfig, DemuxalotConfig,
    VireoConfig, CellSNPConfig
)
from cellranger_snakemake.schemas.doublet_detection import (
    DoubletDetectionConfig, ScrubletConfig
)


class ConfigGenerator:
    """Interactive configuration generator."""
    
    def __init__(self):
        self.console = Console()
    
    def generate_default_config(self) -> PipelineConfig:
        """Generate configuration with all parameters set to default values (no prompts)."""
        self.console.print("\n[bold blue]Generating default configuration with ALL parameters...[/bold blue]\n")
        self.console.print("[yellow]Note: All steps are disabled by default. Edit and enable the ones you need.[/yellow]\n")
        
        # Create config with all pipeline steps showing all available parameters
        # All steps are explicitly disabled so users must enable what they need
        config_dict = {
            "project_name": "my_project",
            "output_dir": "output",
            "resources": ResourceConfig(),
            "directories_suffix": "none",
            
            # Cell Ranger configs with placeholder paths - disabled by default
            "cellranger_gex": CellRangerGEXConfig(
                enabled=False,
                reference="/path/to/cellranger/reference",
                libraries="libraries_gex.tsv"
            ),
            "cellranger_atac": CellRangerATACConfig(
                enabled=False,
                reference="/path/to/cellranger-atac/reference",
                libraries="libraries_atac.tsv"
            ),
            "cellranger_arc": CellRangerARCConfig(
                enabled=False,
                reference="/path/to/cellranger-arc/reference",
                libraries="libraries_arc.tsv"
            ),
            
            # Demultiplexing config - disabled by default
            "demultiplexing": DemultiplexingConfig(
                enabled=False,
                method="demuxalot",
                demuxalot=DemuxalotConfig(
                    vcf="/path/to/genotypes.vcf",
                    genome_names="/path/to/genome_names.txt",
                    refine=False,
                ),
            ),
            
            # Doublet detection config - disabled by default
            "doublet_detection": DoubletDetectionConfig(
                enabled=False,
                method="scrublet",
                scrublet=ScrubletConfig()
            ),
            
        }
        
        return PipelineConfig(**config_dict)
    
    def generate(self) -> PipelineConfig:
        """Generate configuration interactively."""
        self.console.print("\n[bold blue]Single-Cell Preprocessing Pipeline Configuration Generator[/bold blue]\n")
        
        # Project metadata
        project_name = Prompt.ask("Project name")
        output_dir = Prompt.ask("Output directory", default="output")
        
        # Resource configuration
        self.console.print("\n[bold]Resource Configuration[/bold]")
        mem_gb = IntPrompt.ask("Total memory (GB)", default=32)
        tmpdir = Prompt.ask("Temporary directory", default="")
        resources = ResourceConfig(mem_gb=mem_gb, tmpdir=tmpdir)
        
        # Directory suffix
        directories_suffix = Prompt.ask("Directory suffix (or 'none')", default="none")
        
        # Pipeline steps
        config_dict = {
            "project_name": project_name,
            "output_dir": output_dir,
            "resources": resources,
            "directories_suffix": directories_suffix,
        }
        
        # Cell Ranger steps
        self.console.print("\n[bold]Cell Ranger Configuration[/bold]")
        
        if Confirm.ask("Enable Cell Ranger GEX?", default=False):
            config_dict["cellranger_gex"] = self._configure_cellranger_gex()
        
        if Confirm.ask("Enable Cell Ranger ATAC?", default=False):
            config_dict["cellranger_atac"] = self._configure_cellranger_atac()
        
        if Confirm.ask("Enable Cell Ranger ARC?", default=False):
            config_dict["cellranger_arc"] = self._configure_cellranger_arc()
        
        # Demultiplexing
        if Confirm.ask("\nEnable demultiplexing?", default=False):
            config_dict["demultiplexing"] = self._configure_demultiplexing()
        
        # Doublet detection
        if Confirm.ask("\nEnable doublet detection?", default=False):
            config_dict["doublet_detection"] = self._configure_doublet_detection()
        
        # Create and validate config
        return PipelineConfig(**config_dict)
    
    def _configure_cellranger_gex(self) -> CellRangerGEXConfig:
        """Configure Cell Ranger GEX."""
        reference = Prompt.ask("  Reference genome path")
        libraries = Prompt.ask("  Libraries TSV path")
        chemistry = Prompt.ask(
            "  Chemistry",
            choices=["auto", "threeprime", "fiveprime", "SC3Pv1", "SC3Pv2", "SC3Pv3", "SC5P-PE", "SC5P-R2"],
            default="auto"
        )
        normalize = Prompt.ask("  Normalization method", choices=["none", "mapped", "depth"], default="none")
        
        return CellRangerGEXConfig(
            reference=reference,
            libraries=libraries,
            chemistry=chemistry,
            normalize=normalize
        )
    
    def _configure_cellranger_atac(self) -> CellRangerATACConfig:
        """Configure Cell Ranger ATAC."""
        reference = Prompt.ask("  Reference genome path")
        libraries = Prompt.ask("  Libraries TSV path")
        chemistry = Prompt.ask("  Chemistry", choices=["auto", "ARC-v1"], default="auto")
        normalize = Prompt.ask("  Normalization method", choices=["none", "depth"], default="none")
        
        return CellRangerATACConfig(
            reference=reference,
            libraries=libraries,
            chemistry=chemistry,
            normalize=normalize
        )
    
    def _configure_cellranger_arc(self) -> CellRangerARCConfig:
        """Configure Cell Ranger ARC."""
        reference = Prompt.ask("  Reference genome path")
        libraries = Prompt.ask("  Libraries TSV path")
        normalize = Prompt.ask("  Normalization method", choices=["none", "depth"], default="none")
        
        return CellRangerARCConfig(
            reference=reference,
            libraries=libraries,
            normalize=normalize
        )
    
    def _configure_demultiplexing(self) -> DemultiplexingConfig:
        """Configure demultiplexing."""
        method = Prompt.ask(
            "  Demultiplexing method",
            choices=["demuxalot", "vireo"],
            default="demuxalot"
        )

        method_config = None
        if method == "demuxalot":
            vcf = Prompt.ask("    VCF file path")
            genome_names = Prompt.ask("    Genome names (comma-separated)")
            method_config = DemuxalotConfig(
                vcf=vcf,
                genome_names=genome_names,
            )
        
        elif method == "vireo":
            donors = IntPrompt.ask("    Number of donors")
            vcf = Prompt.ask("    Path to VCF file")
            threads = IntPrompt.ask("    Number of threads", default=4)
            min_maf = FloatPrompt.ask("    Minimum MAF", default=0.0)
            min_count = IntPrompt.ask("    Minimum UMI count", default=1)
            umi_tag = Prompt.ask("    UMI tag", default="Auto")
            cell_tag = Prompt.ask("    Cell barcode tag", default="CB")
            cellsnp = CellSNPConfig(
                vcf=vcf,
                threads=threads,
                min_maf=min_maf,
                min_count=min_count,
                umi_tag=umi_tag,
                cell_tag=cell_tag,
                gzip=True
            )
            method_config = VireoConfig(cellsnp=cellsnp, donors=donors)
        
        return DemultiplexingConfig(method=method, **{method: method_config})
    
    def _configure_doublet_detection(self) -> DoubletDetectionConfig:
        """Configure doublet detection."""
        method = Prompt.ask(
            "  Doublet detection method",
            choices=["scrublet"],
            default="scrublet"
        )

        method_config = None
        if method == "scrublet":
            rate = FloatPrompt.ask("    Expected doublet rate", default=0.06)
            filter_cells_min_genes = IntPrompt.ask("    Min genes per cell (filter_cells_min_genes)", default=100)
            filter_genes_min_cells = IntPrompt.ask("    Min cells per gene (filter_genes_min_cells)", default=3)
            min_gene_variability_pctl = FloatPrompt.ask("    Min gene variability percentile", default=85.0)
            n_prin_comps = IntPrompt.ask("    Number of principal components", default=30)
            method_config = ScrubletConfig(
                expected_doublet_rate=rate,
                filter_cells_min_genes=filter_cells_min_genes,
                filter_genes_min_cells=filter_genes_min_cells,
                min_gene_variability_pctl=min_gene_variability_pctl,
                n_prin_comps=n_prin_comps
            )

        return DoubletDetectionConfig(method=method, **{method: method_config})
    
    def save_config(self, config: PipelineConfig, output_path: str):
        """Save configuration to YAML file."""
        # Convert to dict for YAML serialization
        config_dict = config.model_dump(exclude_none=True, by_alias=True)
        
        # Check if file exists and prompt for new name if needed
        final_path = output_path
        if Path(output_path).exists():
            self.console.print(f"\n[yellow]⚠[/yellow]  File '{output_path}' already exists.")
            overwrite = Confirm.ask("Overwrite existing file?", default=False)
            
            if not overwrite:
                while True:
                    new_path = Prompt.ask("Enter a new filename", default=output_path)
                    if not Path(new_path).exists():
                        final_path = new_path
                        self.console.print(f"\n[bold cyan]→[/bold cyan] Using new filename: '{final_path}'")
                        break
                    else:
                        self.console.print(f"[yellow]File '{new_path}' also exists. Please choose another name.[/yellow]")
        
        self.console.print(f"\nSaving configuration to '{final_path}'...")
        
        with open(final_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False, indent=2)
        
        self.console.print(f"\n[bold green]✓[/bold green] Configuration saved to: {final_path}")


def main():
    """Main entry point for interactive config generation."""
    generator = ConfigGenerator()
    
    try:
        config = generator.generate()
        
        # Suggest filename based on enabled Cell Ranger modalities
        enabled_steps = config.get_enabled_steps()
        cr_modalities = []
        if "cellranger_gex" in enabled_steps:
            cr_modalities.append("GEX")
        if "cellranger_atac" in enabled_steps:
            cr_modalities.append("ATAC")
        if "cellranger_arc" in enabled_steps:
            cr_modalities.append("ARC")
        
        if cr_modalities:
            default_filename = f"pipeline_config_{'_'.join(cr_modalities)}.yaml"
        else:
            default_filename = "pipeline_config.yaml"
        
        output_file = Prompt.ask("\nSave configuration to", default=default_filename)
        generator.save_config(config, output_file)
        
        # Show enabled steps
        generator.console.print("\n[bold]Enabled pipeline steps:[/bold]")
        for step in enabled_steps:
            generator.console.print(f"  • {step}")
        
    except KeyboardInterrupt:
        generator.console.print("\n\n[yellow]Configuration generation cancelled.[/yellow]")
        sys.exit(0)
    except Exception as e:
        generator.console.print(f"\n[bold red]Error:[/bold red] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
