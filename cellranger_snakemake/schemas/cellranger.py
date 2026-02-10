"""Cell Ranger configuration schemas for all modalities."""

from typing import ClassVar, Literal, Optional
from pydantic import BaseModel, Field, field_validator
from .base import BaseStepConfig, DirectoryConfig, ToolMeta


class CellRangerGEXConfig(BaseStepConfig):
    """Cell Ranger GEX (Gene Expression) configuration."""

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellranger",
        url="https://www.10xgenomics.com/support/software/cell-ranger/latest",
        shell_version_cmd="cellranger --version",
    )

    reference: str = Field(
        description="Path to Cell Ranger GEX reference genome"
    )
    libraries: str = Field(
        description="Path to libraries TSV file with columns: batch, capture, sample, fastqs"
    )
    chemistry: Literal["auto", "ARC-v1", "threeprime", "fiveprime", "SC3Pv1", "SC3Pv2", "SC3Pv3", "SC5P-PE", "SC5P-R2"] = Field(
        default="auto",
        description="Assay chemistry configuration"
    )
    normalize: Literal["none", "mapped", "depth"] = Field(
        default="none",
        description="Normalization method for aggregation"
    )
    create_bam: bool = Field(
        default=False,
        description="Whether to create a BAM file",
        alias="create-bam"
    )
    threads: int = Field(
        default=10,
        ge=1,
        description="CPU threads for cellranger count"
    )
    mem_gb: int = Field(
        default=64,
        gt=0,
        description="Memory (GB) for cellranger count"
    )
    
    # Cell Ranger job submission options (--jobmode and --mempercore, --maxjobs): https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-job-submission-mode
    jobmode: Optional[Literal["local", "slurm", "sge"]] = Field(
        default=None,
        description="Cell Ranger's --jobmode flag (use Snakemake --profile instead for cluster execution)"
    )
    mempercore: Optional[int] = Field(
        default=None,
        description="GB of RAM per CPU core (Cell Ranger's --mempercore flag)"
    )
    maxjobs: Optional[int] = Field(
        default=None,
        description="Maximum number of jobs (Cell Ranger's --maxjobs flag)"
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

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellranger-atac",
        url="https://www.10xgenomics.com/support/software/cell-ranger-atac/latest",
        shell_version_cmd="cellranger-atac --version",
    )

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
    threads: int = Field(
        default=10,
        ge=1,
        description="CPU threads for cellranger-atac count"
    )
    mem_gb: int = Field(
        default=64,
        gt=0,
        description="Memory (GB) for cellranger-atac count"
    )
    
    # Cell Ranger job submission options (optional, for cellranger's built-in cluster submission)
    jobmode: Optional[Literal["local", "slurm", "sge"]] = Field(
        default=None,
        description="Cell Ranger's --jobmode flag (use Snakemake --profile instead for cluster execution)"
    )
    mempercore: Optional[int] = Field(
        default=None,
        description="GB of RAM per CPU core (Cell Ranger's --mempercore flag)"
    )
    maxjobs: Optional[int] = Field(
        default=None,
        description="Maximum number of jobs (Cell Ranger's --maxjobs flag)"
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

    tool_meta: ClassVar[ToolMeta] = ToolMeta(
        package="cellranger-arc",
        url="https://www.10xgenomics.com/support/software/cell-ranger-arc/latest",
        shell_version_cmd="cellranger-arc --version",
    )

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
    threads: int = Field(
        default=10,
        ge=1,
        description="CPU threads for cellranger-arc count"
    )
    mem_gb: int = Field(
        default=64,
        gt=0,
        description="Memory (GB) for cellranger-arc count"
    )
    
    # Cell Ranger job submission options (optional, for cellranger's built-in cluster submission)
    jobmode: Optional[Literal["local", "slurm", "sge"]] = Field(
        default=None,
        description="Cell Ranger's --jobmode flag (use Snakemake --profile instead for cluster execution)"
    )
    mempercore: Optional[int] = Field(
        default=None,
        description="GB of RAM per CPU core (Cell Ranger's --mempercore flag)"
    )
    maxjobs: Optional[int] = Field(
        default=None,
        description="Maximum number of jobs (Cell Ranger's --maxjobs flag)"
    )
    
    class DirectoryConfig(DirectoryConfig):
        """ARC-specific directories."""
        CELLRANGERARC_COUNT_DIR: str = Field(default="01_CELLRANGERARC_COUNT")
        CELLRANGERARC_AGGR_DIR: str = Field(default="02_CELLRANGERARC_AGGR")
