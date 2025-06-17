from setuptools import setup, find_packages

setup(
    name="cellranger-snakemake",
    version="1.0.0",
    author="Matthew S. Schechter",
    author_email="mschechter@uchicago.edu",
    description="Snakemake wrapper for Cell Ranger workflows",
    packages=find_packages(),
    include_package_data=True,  # <--- Add this line
    package_data={
        "cellranger_snakemake.workflows.ARC": ["*.smk"],
        "cellranger_snakemake.workflows.ATAC": ["*.smk"],
        "cellranger_snakemake.workflows.GEX": ["*.smk"],
    },
    install_requires=[
        # Add other dependencies here
    ],
    entry_points={
        "console_scripts": [
            "snakemake-run-cellranger=cellranger_snakemake.cli:main",
        ],
    },
    python_requires='>=3.7',
)
