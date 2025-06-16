from setuptools import setup, find_packages

setup(
    name="cellranger-snakemake",
    version="1.0.0",
    author="Matthew S. Schechter",
    author_email="mschechter@uchicago.edu",
    description="Snakemake wrapper for Cell Ranger workflows",
    packages=find_packages(),  # Automatically find all packages under the directory
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
