"""Cell type annotation configuration schemas."""

from typing import ClassVar, Literal, Optional, List
from pydantic import BaseModel, Field, model_validator
from .base import BaseStepConfig, ToolMeta


class CelltypistConfig(BaseModel):
    """Celltypist annotation parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="celltypist",
        url="https://github.com/Teichlab/celltypist",
    )

    model: str = Field(
        description="Path to celltypist model file or model name from celltypist models"
    )
    majority_voting: bool = Field(
        default=False,
        description="Use majority voting for cell type assignment"
    )
    over_clustering: Optional[str] = Field(
        default=None,
        description="Column name for over-clustering results (optional)"
    )
    min_prop: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Minimum proportion for majority voting"
    )
    
    class Config:
        extra = "forbid"


class ScANVIConfig(BaseModel):
    """scANVI (scvi-tools) reference-based annotation parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="scvi-tools",
        url="https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html",
    )

    reference_path: str = Field(
        description="Path to annotated reference AnnData (.h5ad file)"
    )
    label_key: str = Field(
        description="Column name in reference.obs containing cell type labels"
    )
    n_hidden: int = Field(
        default=128,
        ge=1,
        description="Number of hidden units"
    )
    n_latent: int = Field(
        default=30,
        ge=1,
        description="Latent space dimensionality"
    )
    n_layers: int = Field(
        default=2,
        ge=1,
        description="Number of hidden layers"
    )
    max_epochs: int = Field(
        default=400,
        ge=1,
        description="Maximum training epochs"
    )

    class Config:
        extra = "forbid"


class DecouplerMarkerConfig(BaseModel):
    """Decoupler marker-based annotation parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="decoupler",
        url="https://decoupler-py.readthedocs.io/",
    )

    marker_database: str = Field(
        description="Path to marker gene set file (GMT/CSV format) or tissue type for built-in markers"
    )
    method: Literal["ulm", "wsum", "ora", "aucell"] = Field(
        default="ulm",
        description="Scoring method: ulm (univariate linear model), wsum (weighted sum), ora (over-representation), aucell"
    )
    min_score: float = Field(
        default=0.5,
        description="Minimum score threshold for cell type assignment"
    )

    class Config:
        extra = "forbid"


class CelltypistCustomConfig(BaseModel):
    """Custom celltypist model training parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="celltypist",
        url="https://github.com/Teichlab/celltypist",
    )

    training_data: str = Field(
        description="Path to training data (h5ad file with cell type labels)"
    )
    label_column: str = Field(
        description="Column name containing cell type labels"
    )
    feature_selection: bool = Field(
        default=True,
        description="Perform feature selection for model training"
    )
    n_jobs: int = Field(
        default=1,
        ge=1,
        description="Number of parallel jobs for training"
    )
    
    class Config:
        extra = "forbid"


class CelltypeAnnotationConfig(BaseStepConfig):
    """Cell type annotation step configuration (Python-only)."""

    method: Literal["celltypist", "scanvi", "decoupler_markers", "celltypist_custom"] = Field(
        description="Cell type annotation method to use"
    )

    # Method-specific parameters
    celltypist: Optional[CelltypistConfig] = None
    scanvi: Optional[ScANVIConfig] = None
    decoupler_markers: Optional[DecouplerMarkerConfig] = None
    celltypist_custom: Optional[CelltypistCustomConfig] = None

    @model_validator(mode='after')
    def validate_method_params(self):
        """Ensure the correct parameters are provided for the selected method."""
        method_configs = {
            "celltypist": self.celltypist,
            "scanvi": self.scanvi,
            "decoupler_markers": self.decoupler_markers,
            "celltypist_custom": self.celltypist_custom,
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
            print(f"Warning: Ignoring unused annotation configs: {list(other_configs.keys())}")

        return self
