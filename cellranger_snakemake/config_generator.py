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
    DemultiplexingConfig, DemuxletConfig, SouporcellConfig, 
    FreemuxletConfig, ViralConfig
)
from cellranger_snakemake.schemas.doublet_detection import (
    DoubletDetectionConfig, ScrubletConfig, DoubletFinderConfig,
    ScdsConfig, ScDblFinderConfig
)
from cellranger_snakemake.schemas.annotation import (
    CelltypeAnnotationConfig, CelltypistConfig, AzimuthConfig,
    SingleRConfig, ScTypeConfig, CelltypistCustomConfig
)


class ConfigGenerator:
    """Interactive configuration generator."""
    
    def __init__(self):
        self.console = Console()
    
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
        
        # Cell type annotation
        if Confirm.ask("\nEnable cell type annotation?", default=False):
            config_dict["celltype_annotation"] = self._configure_celltype_annotation()
        
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
            choices=["demuxlet", "souporcell", "freemuxlet", "vireo"],
            default="demuxlet"
        )
        
        method_config = None
        if method == "demuxlet":
            vcf = Prompt.ask("    VCF file path")
            field = Prompt.ask("    VCF field", default="GT")
            alpha = FloatPrompt.ask("    Alpha parameter", default=0.5)
            method_config = DemuxletConfig(vcf=vcf, field=field, alpha=alpha)
        
        elif method == "souporcell":
            clusters = IntPrompt.ask("    Number of clusters")
            ploidy = IntPrompt.ask("    Ploidy", default=2)
            method_config = SouporcellConfig(clusters=clusters, ploidy=ploidy)
        
        elif method == "freemuxlet":
            nsample = IntPrompt.ask("    Number of samples")
            alpha = FloatPrompt.ask("    Alpha parameter", default=0.5)
            method_config = FreemuxletConfig(nsample=nsample, alpha=alpha)
        
        elif method == "vireo":
            donors = IntPrompt.ask("    Number of donors")
            method_config = ViralConfig(donors=donors)
        
        return DemultiplexingConfig(method=method, **{method: method_config})
    
    def _configure_doublet_detection(self) -> DoubletDetectionConfig:
        """Configure doublet detection."""
        method = Prompt.ask(
            "  Doublet detection method",
            choices=["scrublet", "doubletfinder", "scds", "scdblfinder"],
            default="scrublet"
        )
        
        method_config = None
        if method == "scrublet":
            rate = FloatPrompt.ask("    Expected doublet rate", default=0.06)
            min_counts = IntPrompt.ask("    Minimum counts", default=2)
            min_cells = IntPrompt.ask("    Minimum cells", default=3)
            method_config = ScrubletConfig(
                expected_doublet_rate=rate,
                min_counts=min_counts,
                min_cells=min_cells
            )
        
        elif method == "doubletfinder":
            pN = FloatPrompt.ask("    pN parameter", default=0.25)
            pK = FloatPrompt.ask("    pK parameter", default=0.09)
            method_config = DoubletFinderConfig(pN=pN, pK=pK)
        
        elif method == "scds":
            algorithm = Prompt.ask("    Algorithm", choices=["cxds", "bcds", "hybrid"], default="hybrid")
            method_config = ScdsConfig(algorithm=algorithm)
        
        elif method == "scdblfinder":
            method_config = ScDblFinderConfig()
        
        return DoubletDetectionConfig(method=method, **{method: method_config})
    
    def _configure_celltype_annotation(self) -> CelltypeAnnotationConfig:
        """Configure cell type annotation."""
        method = Prompt.ask(
            "  Annotation method",
            choices=["celltypist", "azimuth", "singler", "sctype", "celltypist_custom"],
            default="celltypist"
        )
        
        method_config = None
        if method == "celltypist":
            model = Prompt.ask("    Model name or path")
            majority_voting = Confirm.ask("    Use majority voting?", default=False)
            method_config = CelltypistConfig(model=model, majority_voting=majority_voting)
        
        elif method == "azimuth":
            reference = Prompt.ask(
                "    Reference dataset",
                choices=["pbmc", "bonemarrow", "lung", "kidney", "heart", "adipose", "liver", "pancreas", "motor_cortex", "fetal"],
                default="pbmc"
            )
            method_config = AzimuthConfig(reference=reference)
        
        elif method == "singler":
            reference = Prompt.ask("    Reference object path")
            labels = Prompt.ask("    Label column", default="label.main")
            method_config = SingleRConfig(reference=reference, labels=labels)
        
        elif method == "sctype":
            tissue = Prompt.ask("    Tissue type")
            method_config = ScTypeConfig(tissue=tissue)
        
        elif method == "celltypist_custom":
            training_data = Prompt.ask("    Training data path")
            label_column = Prompt.ask("    Label column")
            method_config = CelltypistCustomConfig(
                training_data=training_data,
                label_column=label_column
            )
        
        return CelltypeAnnotationConfig(method=method, **{method: method_config})
    
    def save_config(self, config: PipelineConfig, output_path: str):
        """Save configuration to YAML file."""
        # Convert to dict for YAML serialization
        config_dict = config.model_dump(exclude_none=True, mode='python')
        
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
