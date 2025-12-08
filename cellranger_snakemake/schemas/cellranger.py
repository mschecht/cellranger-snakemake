"""Cell Ranger configuration schemas for all modalities."""

from typing import Literal, Optional
from pydantic import BaseModel, Field, field_validator
from .base import BaseStepConfig, DirectoryConfig


class CellRangerGEXConfig(BaseStepConfig):
    """Cell Ranger GEX (Gene Expression) configuration."""
    
    reference: str = Field(
        description="Path to Cell Ranger GEX reference genome"
    )
    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, sample, fastqs"
    )
    chemistry: Literal["auto", "threeprime", "fiveprime", "SC3Pv1", "SC3Pv2", "SC3Pv3", "SC5P-PE", "SC5P-R2"] = Field(
        default="auto",
        description="Assay chemistry configuration"
    )
    normalize: Literal["none", "mapped", "depth"] = Field(
        default="none",
        description="Normalization method for aggregation"
    )
    directories: DirectoryConfig = Field(
        default_factory=lambda: DirectoryConfig(
            LOGS_DIR="00_LOGS",
            CELLRANGERGEX_COUNT_DIR="01_CELLRANGERGEX_COUNT",
            CELLRANGERGEX_AGGR_DIR="02_CELLRANGERGEX_AGGR"
        )
    )
    
    class DirectoryConfig(DirectoryConfig):
        """GEX-specific directories."""
        CELLRANGERGEX_COUNT_DIR: str = Field(default="01_CELLRANGERGEX_COUNT")
        CELLRANGERGEX_AGGR_DIR: str = Field(default="02_CELLRANGERGEX_AGGR")


class CellRangerATACConfig(BaseStepConfig):
    """Cell Ranger ATAC configuration."""
    
    reference: str = Field(
        description="Path to Cell Ranger ATAC reference genome"
    )
    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, sample, fastqs"
    )
    chemistry: Literal["auto", "ARC-v1"] = Field(
        default="auto",
        description="Assay chemistry configuration"
    )
    normalize: Literal["none", "depth"] = Field(
        default="none",
        description="Normalization method for aggregation"
    )
    directories: DirectoryConfig = Field(
        default_factory=lambda: DirectoryConfig(
            LOGS_DIR="00_LOGS",
            CELLRANGERATAC_COUNT_DIR="01_CELLRANGERATAC_COUNT",
            CELLRANGERATAC_AGGR_DIR="02_CELLRANGERATAC_AGGR"
        )
    )
    
    class DirectoryConfig(DirectoryConfig):
        """ATAC-specific directories."""
        CELLRANGERATAC_COUNT_DIR: str = Field(default="01_CELLRANGERATAC_COUNT")
        CELLRANGERATAC_AGGR_DIR: str = Field(default="02_CELLRANGERATAC_AGGR")


class CellRangerARCConfig(BaseStepConfig):
    """Cell Ranger ARC (multiome: ATAC + RNA) configuration."""
    
    reference: str = Field(
        description="Path to Cell Ranger ARC reference genome"
    )
    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, CSV"
    )
    normalize: Literal["none", "depth"] = Field(
        default="none",
        description="Normalization method for aggregation"
    )
    directories: DirectoryConfig = Field(
        default_factory=lambda: DirectoryConfig(
            LOGS_DIR="00_LOGS",
            CELLRANGERARC_COUNT_DIR="01_CELLRANGERARC_COUNT",
            CELLRANGERARC_AGGR_DIR="02_CELLRANGERARC_AGGR"
        )
    )
    
    class DirectoryConfig(DirectoryConfig):
        """ARC-specific directories."""
        CELLRANGERARC_COUNT_DIR: str = Field(default="01_CELLRANGERARC_COUNT")
        CELLRANGERARC_AGGR_DIR: str = Field(default="02_CELLRANGERARC_AGGR")
