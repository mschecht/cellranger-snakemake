"""Doublet detection configuration schemas."""

from typing import Literal, Optional
from pydantic import BaseModel, Field, model_validator
from .base import BaseStepConfig


class ScrubletConfig(BaseModel):
    """Scrublet doublet detection parameters."""
    
    expected_doublet_rate: float = Field(
        default=0.06,
        ge=0.0,
        le=1.0,
        description="Expected doublet rate (typically 0.05-0.10)"
    )
    min_counts: int = Field(
        default=2,
        ge=1,
        description="Minimum UMI counts per cell"
    )
    min_cells: int = Field(
        default=3,
        ge=1,
        description="Minimum cells expressing a gene"
    )
    min_gene_variability_pctl: float = Field(
        default=85.0,
        ge=0.0,
        le=100.0,
        description="Percentile for gene variability filtering"
    )
    n_prin_comps: int = Field(
        default=30,
        ge=1,
        description="Number of principal components"
    )
    
    class Config:
        extra = "forbid"


class DoubletFinderConfig(BaseModel):
    """DoubletFinder doublet detection parameters."""
    
    pN: float = Field(
        default=0.25,
        ge=0.0,
        le=1.0,
        description="Number of generated artificial doublets"
    )
    pK: float = Field(
        default=0.09,
        ge=0.0,
        le=1.0,
        description="PC neighborhood size"
    )
    nExp: Optional[int] = Field(
        default=None,
        description="Number of expected doublets (if None, calculated from pN)"
    )
    sct: bool = Field(
        default=False,
        description="Use SCTransform normalized data"
    )
    
    class Config:
        extra = "forbid"


class ScdsConfig(BaseModel):
    """scds (cxds/bcds/hybrid) doublet detection parameters."""
    
    algorithm: Literal["cxds", "bcds", "hybrid"] = Field(
        default="hybrid",
        description="Which scds algorithm to use"
    )
    ntop: int = Field(
        default=500,
        ge=1,
        description="Number of top variable genes to use"
    )
    
    class Config:
        extra = "forbid"


class ScDblFinderConfig(BaseModel):
    """scDblFinder doublet detection parameters."""
    
    dbr: Optional[float] = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="Expected doublet rate (if None, estimated from data)"
    )
    clusters: Optional[str] = Field(
        default=None,
        description="Column name for cluster labels (optional)"
    )
    
    class Config:
        extra = "forbid"


class DoubletDetectionConfig(BaseStepConfig):
    """Doublet detection step configuration."""
    
    method: Literal["scrublet", "doubletfinder", "scds", "scdblfinder"] = Field(
        description="Doublet detection method to use"
    )
    
    # Method-specific parameters
    scrublet: Optional[ScrubletConfig] = None
    doubletfinder: Optional[DoubletFinderConfig] = None
    scds: Optional[ScdsConfig] = None
    scdblfinder: Optional[ScDblFinderConfig] = None
    
    @model_validator(mode='after')
    def validate_method_params(self):
        """Ensure the correct parameters are provided for the selected method."""
        method_configs = {
            "scrublet": self.scrublet,
            "doubletfinder": self.doubletfinder,
            "scds": self.scds,
            "scdblfinder": self.scdblfinder,
        }
        
        selected_config = method_configs.get(self.method)
        if selected_config is None:
            raise ValueError(
                f"Parameters for method '{self.method}' are required. "
                f"Please provide a '{self.method}' configuration block."
            )
        
        # Warn if other method configs are present
        other_configs = {k: v for k, v in method_configs.items() if k != self.method and v is not None}
        if other_configs:
            print(f"Warning: Ignoring unused doublet detection configs: {list(other_configs.keys())}")
        
        return self
