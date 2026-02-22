"""Doublet detection configuration schemas."""

from typing import ClassVar, Literal, Optional
from pydantic import BaseModel, Field, model_validator
from .base import BaseStepConfig, ToolMeta


class ScrubletConfig(BaseModel):
    """Scrublet doublet detection parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="scrublet",
        url="https://github.com/swolock/scrublet",
    )

    filter_cells_min_genes: int = Field(
        default=100,
        ge=1,
        description="sc.pp.filter_cells(min_genes=): remove cells with fewer than this many genes detected"
    )
    filter_genes_min_cells: int = Field(
        default=3,
        ge=1,
        description="sc.pp.filter_genes(min_cells=): remove genes detected in fewer than this many cells"
    )
    expected_doublet_rate: float = Field(
        default=0.06,
        ge=0.0,
        le=1.0,
        description="sc.pp.scrublet(expected_doublet_rate=): expected fraction of doublets (typically 0.05-0.10)"
    )
    min_gene_variability_pctl: float = Field(
        default=85.0,
        ge=0.0,
        le=100.0,
        description="sc.pp.scrublet(min_gene_variability_pctl=): percentile threshold for highly variable gene selection"
    )
    n_prin_comps: int = Field(
        default=30,
        ge=1,
        description="sc.pp.scrublet(n_prin_comps=): number of principal components for doublet detection"
    )
    sim_doublet_ratio: float = Field(
        default=2.0,
        gt=0.0,
        description="sc.pp.scrublet(sim_doublet_ratio=): ratio of simulated doublets to observed transcriptomes"
    )
    threshold: Optional[float] = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="sc.pp.scrublet(threshold=): doublet score cutoff; auto-determined if None"
    )
    n_neighbors: Optional[int] = Field(
        default=None,
        ge=1,
        description="sc.pp.scrublet(n_neighbors=): KNN neighbors; auto-set to 0.5*sqrt(n_obs) if None"
    )
    random_state: int = Field(
        default=0,
        ge=0,
        description="sc.pp.scrublet(random_state=): random seed for reproducibility"
    )

    class Config:
        extra = "forbid"


class SoloConfig(BaseModel):
    """SOLO (from scvi-tools) doublet detection parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="scvi-tools",
        url="https://docs.scvi-tools.org/en/stable/user_guide/models/solo.html",
    )

    n_hidden: int = Field(
        default=128,
        ge=1,
        description="Number of hidden units"
    )
    n_latent: int = Field(
        default=64,
        ge=1,
        description="Latent space dimensionality"
    )
    n_layers: int = Field(
        default=1,
        ge=1,
        description="Number of hidden layers"
    )
    learning_rate: float = Field(
        default=1e-3,
        gt=0.0,
        description="Learning rate"
    )
    max_epochs: int = Field(
        default=400,
        ge=1,
        description="Maximum training epochs"
    )

    class Config:
        extra = "forbid"


class DoubletDetectionConfig(BaseStepConfig):
    """Doublet detection step configuration (Python-only)."""

    method: Literal["scrublet", "solo"] = Field(
        description="Doublet detection method to use"
    )

    # Method-specific parameters
    scrublet: Optional[ScrubletConfig] = None
    solo: Optional[SoloConfig] = None

    @model_validator(mode='after')
    def validate_method_params(self):
        """Ensure the correct parameters are provided for the selected method."""
        method_configs = {
            "scrublet": self.scrublet,
            "solo": self.solo,
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
