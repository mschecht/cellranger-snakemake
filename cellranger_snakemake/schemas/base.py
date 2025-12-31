"""Base configuration models and shared types."""

from typing import Optional, Literal
from pydantic import BaseModel, Field, field_validator


class ResourceConfig(BaseModel):
    """Computational resource configuration."""
    mem_gb: int = Field(default=32, description="Memory in GB")
    tmpdir: str = Field(default="", description="Temporary directory path")


class DirectoryConfig(BaseModel):
    """Output directory configuration."""
    LOGS_DIR: str = Field(default="00_LOGS", description="Logs directory")
    
    @field_validator('*')
    @classmethod
    def validate_no_spaces(cls, v: str) -> str:
        """Ensure directory names don't contain spaces."""
        if ' ' in v:
            raise ValueError(f"Directory names cannot contain spaces: {v}")
        return v


class BaseStepConfig(BaseModel):
    """Base configuration for pipeline steps."""
    enabled: bool = Field(default=True, description="Whether this step is enabled")
    
    class Config:
        extra = "forbid"  # Prevent unexpected parameters
        validate_assignment = True
