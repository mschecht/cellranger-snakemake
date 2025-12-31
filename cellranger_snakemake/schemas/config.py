"""Main unified configuration schema for single-cell preprocessing pipeline."""

from typing import Optional, Dict, Any
from pydantic import BaseModel, Field, field_validator
from .base import ResourceConfig
from .cellranger import CellRangerGEXConfig, CellRangerATACConfig, CellRangerARCConfig
from .demultiplexing import DemultiplexingConfig
from .doublet_detection import DoubletDetectionConfig
from .annotation import CelltypeAnnotationConfig


class SampleMetadata(BaseModel):
    """Metadata for a single sample."""
    
    sample_id: str = Field(description="Unique sample identifier")
    batch: Optional[str] = Field(default=None, description="Batch identifier")
    description: Optional[str] = Field(default=None, description="Sample description")
    
    # Paths to data
    gex_fastqs: Optional[str] = Field(default=None, description="Path to GEX FASTQ files")
    atac_fastqs: Optional[str] = Field(default=None, description="Path to ATAC FASTQ files")
    arc_csv: Optional[str] = Field(default=None, description="Path to Cell Ranger ARC libraries CSV")
    
    # Sample-specific parameters
    expected_cells: Optional[int] = Field(default=None, description="Expected number of cells")
    
    class Config:
        extra = "allow"  # Allow additional metadata fields


class PipelineConfig(BaseModel):
    """Main unified configuration for single-cell preprocessing pipeline."""
    
    # Project metadata
    project_name: str = Field(description="Name of the project")
    output_dir: str = Field(default="output", description="Base output directory")
    
    # Sample information
    samples: Dict[str, SampleMetadata] = Field(
        default_factory=dict,
        description="Dictionary of sample_id -> sample metadata"
    )
    
    # Global settings
    resources: ResourceConfig = Field(
        default_factory=ResourceConfig,
        description="Computational resources"
    )
    directories_suffix: str = Field(
        default="none",
        description="Suffix to add to all output directories"
    )
    
    # Pipeline steps - all optional, enabled by default if configured
    cellranger_gex: Optional[CellRangerGEXConfig] = Field(
        default=None,
        description="Cell Ranger GEX configuration"
    )
    cellranger_atac: Optional[CellRangerATACConfig] = Field(
        default=None,
        description="Cell Ranger ATAC configuration"
    )
    cellranger_arc: Optional[CellRangerARCConfig] = Field(
        default=None,
        description="Cell Ranger ARC configuration"
    )
    demultiplexing: Optional[DemultiplexingConfig] = Field(
        default=None,
        description="Demultiplexing configuration"
    )
    doublet_detection: Optional[DoubletDetectionConfig] = Field(
        default=None,
        description="Doublet detection configuration"
    )
    celltype_annotation: Optional[CelltypeAnnotationConfig] = Field(
        default=None,
        description="Cell type annotation configuration"
    )
    
    class Config:
        extra = "forbid"
        validate_assignment = True
    
    @field_validator('output_dir')
    @classmethod
    def validate_output_dir(cls, v: str) -> str:
        """Ensure output directory path is valid."""
        if not v or v.isspace():
            raise ValueError("output_dir cannot be empty")
        return v
    
    def get_enabled_steps(self) -> list[str]:
        """Return list of enabled pipeline steps."""
        steps = []
        step_configs = {
            "cellranger_gex": self.cellranger_gex,
            "cellranger_atac": self.cellranger_atac,
            "cellranger_arc": self.cellranger_arc,
            "demultiplexing": self.demultiplexing,
            "doublet_detection": self.doublet_detection,
            "celltype_annotation": self.celltype_annotation,
        }
        
        for step_name, step_config in step_configs.items():
            if step_config is not None and step_config.enabled:
                steps.append(step_name)
        
        return steps
