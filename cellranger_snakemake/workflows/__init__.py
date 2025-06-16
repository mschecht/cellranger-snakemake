#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
from cellranger_snakemake.utils.logger import logger

__all__ = [
    "logger",
    "write_default_config",
    "run_snakemake",
    "sanity_check",
    "__version__",
    "__description__",
]

__version__ = "1.0.0"
__description__ = "Snakemake wrapper for Cell Ranger workflows"

