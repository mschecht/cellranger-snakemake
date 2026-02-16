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
