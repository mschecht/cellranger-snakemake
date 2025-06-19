# config_templates.py

# Default configuration template for ARC pipeline
ARC_CONFIG = {
    "reference": "/path/to/reference-genome",
    "libraries": "/path/to/libraries_list.tsv",
    "HPC_mode": "",
    "mempercore": 8,
    "normalize": "none",
    "directories_prefix": "none",
    "directories": {
        "LOGS_DIR": "00_LOGS",
        "CELLRANGERARC_COUNT_DIR": "01_CELLRANGERARC_COUNT",
        "CELLRANGERARC_AGGR_DIR": "02_CELLRANGERARC_AGGR"
    },
}

ATAC_CONFIG = {
    "samples": {
        "sample1": {
            "fastqs": "/path/to/fastqs/sample1",
            "libraries_csv": "/path/to/sample1_libraries.csv"
        }
    },
    "reference": "/path/to/cellranger-arc-reference",
    "output_dir": "results",
    "cellranger_arc": {
        "extra_args": "--force-cells=5000",
        "localcores": 8,
        "localmem": 64
    },
    "resources": {
        "default_threads": 8,
        "default_mem_gb": 64
    }
}

GEX_CONFIG = {
    "samples": {
        "sample1": {
            "fastqs": "/path/to/fastqs/sample1",
            "libraries_csv": "/path/to/sample1_libraries.csv"
        }
    },
    "reference": "/path/to/cellranger-arc-reference",
    "output_dir": "results",
    "cellranger_arc": {
        "extra_args": "--force-cells=5000",
        "localcores": 8,
        "localmem": 64
    },
    "resources": {
        "default_threads": 8,
        "default_mem_gb": 64
    }
}