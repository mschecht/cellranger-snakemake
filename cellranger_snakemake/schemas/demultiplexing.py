"""Demultiplexing configuration schemas."""

from typing import Literal, Optional, Union
from pydantic import BaseModel, Field, model_validator
from .base import BaseStepConfig


class DemuxalotConfig(BaseModel):
    """demuxalot demultiplexing parameters."""
    
    vcf: str = Field(description="Path to VCF file with genotypes")
    genome_names: str = Field(description="Path to txt file with genome names")
    refine: bool = Field(description="Run genotype refinement step. Read more here: https://github.com/arogozhnikov/demuxalot?tab=readme-ov-file#running-complex-scenario")
    celltag: Literal["CB"] = Field(
        default="CB",
        description="Cell barcode tag in BAM file"
    )
    umitag: Literal["UB"] = Field(
        default="UB",
        description="UMI tag in BAM file"
    )
    
    class Config:
        extra = "forbid"


class CellSNPConfig(BaseModel):
    """cellsnp-lite parameters for SNP calling."""
    
    vcf: str = Field(description="Path to VCF reference file with known variants")
    threads: int = Field(default=4, ge=1, description="Number of threads for cellsnp-lite")
    min_maf: float = Field(default=0.0, ge=0.0, le=1.0, description="Minimum minor allele frequency")
    min_count: int = Field(default=1, ge=0, description="Minimum UMI count")
    umi_tag: str = Field(default="Auto", description="UMI tag in BAM file (e.g., 'Auto', 'UB')")
    cell_tag: str = Field(default="CB", description="Cell barcode tag in BAM file")
    gzip: bool = Field(default=True, description="Gzip output files")
    
    class Config:
        extra = "forbid"


class VireoConfig(BaseModel):
    """Vireo demultiplexing parameters (requires cellsnp-lite preprocessing)."""
    
    cellsnp: CellSNPConfig = Field(description="cellsnp-lite configuration for SNP calling")
    donors: int = Field(description="Number of donors to demultiplex")
    
    class Config:
        extra = "forbid"


class DemultiplexingConfig(BaseStepConfig):
    """Demultiplexing step configuration."""
    
    method: Literal["demuxalot", "vireo"] = Field(
        description="Demultiplexing method to use"
    )
    
    # Method-specific parameters (only one should be populated based on method)
    demuxalot: Optional[DemuxalotConfig] = None
    vireo: Optional[VireoConfig] = None
    
    @model_validator(mode='after')
    def validate_method_params(self):
        """Ensure the correct parameters are provided for the selected method."""
        method_configs = {
            "demuxalot": self.demuxalot,
            "vireo": self.vireo,
        }
        
        selected_config = method_configs.get(self.method)
        if selected_config is None:
            raise ValueError(
                f"Parameters for method '{self.method}' are required. "
                f"Please provide a '{self.method}' configuration block."
            )
        
        # Warn if other method configs are present (informational)
        other_configs = {k: v for k, v in method_configs.items() if k != self.method and v is not None}
        if other_configs:
            print(f"Warning: Ignoring unused demultiplexing configs: {list(other_configs.keys())}")
        
        return self
