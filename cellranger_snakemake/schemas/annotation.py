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


class AzimuthConfig(BaseModel):
    """Azimuth reference-based annotation parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="azimuth",
        url="https://github.com/satijalab/azimuth",
    )

    reference: Literal[
        "pbmc", "bonemarrow", "lung", "kidney", "heart", 
        "adipose", "liver", "pancreas", "motor_cortex", "fetal"
    ] = Field(
        description="Azimuth reference dataset to use"
    )
    
    class Config:
        extra = "forbid"


class SingleRConfig(BaseModel):
    """SingleR reference-based annotation parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="SingleR",
        url="https://github.com/dviraran/SingleR",
    )

    reference: str = Field(
        description="Path to SingleR reference object (RDS file) or built-in reference name"
    )
    labels: str = Field(
        default="label.main",
        description="Column name in reference for cell type labels"
    )
    de_method: Literal["classic", "wilcox", "t"] = Field(
        default="classic",
        description="Method for marker gene detection"
    )
    
    class Config:
        extra = "forbid"


class ScTypeConfig(BaseModel):
    """ScType marker-based annotation parameters."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="sctype",
        url="https://github.com/IanevskiAleksandr/sc-type",
    )

    tissue: str = Field(
        description="Tissue type (e.g., 'Immune system', 'Liver', 'Brain')"
    )
    marker_database: Optional[str] = Field(
        default=None,
        description="Path to custom marker gene database (uses built-in if None)"
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
    """Cell type annotation step configuration."""
    
    method: Literal["celltypist", "azimuth", "singler", "sctype", "celltypist_custom"] = Field(
        description="Cell type annotation method to use"
    )
    
    # Method-specific parameters
    celltypist: Optional[CelltypistConfig] = None
    azimuth: Optional[AzimuthConfig] = None
    singler: Optional[SingleRConfig] = None
    sctype: Optional[ScTypeConfig] = None
    celltypist_custom: Optional[CelltypistCustomConfig] = None
    
    @model_validator(mode='after')
    def validate_method_params(self):
        """Ensure the correct parameters are provided for the selected method."""
        method_configs = {
            "celltypist": self.celltypist,
            "azimuth": self.azimuth,
            "singler": self.singler,
            "sctype": self.sctype,
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
