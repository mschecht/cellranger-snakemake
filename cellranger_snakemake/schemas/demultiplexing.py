"""Demultiplexing configuration schemas."""

from typing import Literal, Optional, Union
from pydantic import BaseModel, Field, model_validator
from .base import BaseStepConfig


class DemuxletConfig(BaseModel):
    """Demuxlet demultiplexing parameters."""
    
    vcf: str = Field(description="Path to VCF file with genotypes")
    field: str = Field(default="GT", description="VCF field to use for genotypes")
    group_list: Optional[str] = Field(
        default=None,
        description="Path to file mapping sample IDs to groups"
    )
    alpha: float = Field(default=0.5, description="Grid parameter for demuxlet")
    
    class Config:
        extra = "forbid"


class SouporcellConfig(BaseModel):
    """Souporcell demultiplexing parameters."""
    
    clusters: int = Field(description="Number of genotypes to cluster")
    known_genotypes: Optional[str] = Field(
        default=None,
        description="Path to VCF file with known genotypes (optional)"
    )
    ploidy: int = Field(default=2, description="Ploidy of the samples")
    min_alt: int = Field(default=10, description="Minimum ALT UMIs for clustering")
    min_ref: int = Field(default=10, description="Minimum REF UMIs for clustering")
    
    class Config:
        extra = "forbid"


class FreemuxletConfig(BaseModel):
    """Freemuxlet demultiplexing parameters."""
    
    nsample: int = Field(description="Number of samples multiplexed together")
    group_list: Optional[str] = Field(
        default=None,
        description="Path to file mapping sample IDs to groups"
    )
    alpha: float = Field(default=0.5, description="Grid parameter")
    
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
    
    method: Literal["demuxlet", "souporcell", "freemuxlet", "vireo"] = Field(
        description="Demultiplexing method to use"
    )
    
    # Method-specific parameters (only one should be populated based on method)
    demuxlet: Optional[DemuxletConfig] = None
    souporcell: Optional[SouporcellConfig] = None
    freemuxlet: Optional[FreemuxletConfig] = None
    vireo: Optional[VireoConfig] = None
    
    @model_validator(mode='after')
    def validate_method_params(self):
        """Ensure the correct parameters are provided for the selected method."""
        method_configs = {
            "demuxlet": self.demuxlet,
            "souporcell": self.souporcell,
            "freemuxlet": self.freemuxlet,
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
