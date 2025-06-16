import os
import sys
import yaml
import subprocess

from pathlib import Path
from cellranger_snakemake.utils.logger import logger
from cellranger_snakemake.config_templates import ARC_CONFIG, ATAC_CONFIG, GEX_CONFIG

class Workflow:
    def __init__(self, name, get_default_config, config_file):
        self.name = name
        self.get_default_config_flag = bool(get_default_config)
        self.default_config_filename = None
        self.config_file = config_file

        if isinstance(get_default_config, str):
            self.default_config_filename = get_default_config
        elif get_default_config is True:
            self.default_config_filename = f"{self.name}_default_config.yaml"

        self.default_config = self._get_default_config()
        self.snakefile = self._resolve_snakefile()

        self.sanity_check()

    def _get_default_config(self):
        if self.name == "ARC":
            return ARC_CONFIG
        elif self.name == "ATAC":
            return ATAC_CONFIG
        elif self.name == "GEX":
            return GEX_CONFIG
        else:
            logger.error(f"Unknown workflow type: {self.name}")
            sys.exit(1)

    def _resolve_snakefile(self):
        snakefile_name = f"cellranger{self.name}.smk"
        path = Path(__file__).parent / self.name / snakefile_name
        if not path.exists():
            logger.error(f"Could not find Snakefile at expected path: {path}")
            sys.exit(1)
        return str(path)

    def sanity_check(self):
        """Validate default config filename, if present, and other preconditions."""
        if self.default_config_filename:
            if not isinstance(self.default_config_filename, str):
                logger.error("`--get-default-config` must be a string filename or True.")
                sys.exit(1)
            if not self.default_config_filename.endswith('.yaml'):
                logger.error("`--get-default-config` must be a filename ending in '.yaml'")
                sys.exit(1)
            if os.path.exists(self.default_config_filename):
                logger.error(f"File '{self.default_config_filename}' already exists and will not be overwritten.")
                sys.exit(1)


    def write_default_config(self):
        """Write default config only if requested."""
        if self.default_config_filename:
            logger.info(f"Writing '{self.name}' default config YAML to '{self.default_config_filename}'")
            with open(self.default_config_filename, 'w') as f:
                yaml.dump(self.default_config, f, indent=2)
    

    def run(self, config_file=None, snakefile=None, dry_run=False, cores=1, dag=False, additional_args=""):
        config_file = config_file or self.config_file
        snakefile = snakefile or self.snakefile

        if not os.path.exists(config_file):
            logger.error(f"Configuration file '{config_file}' not found.")
            sys.exit(1)

        if not os.path.exists(snakefile):
            logger.error(f"Snakefile '{snakefile}' not found.")
            sys.exit(1)

        if dag:
            # Export DAG to PDF instead of running workflow
            output_file = f"{self.name}_dag.pdf"
            logger.info(f"Exporting DAG to {output_file}")
            try:
                with open(output_file, "wb") as out_pdf:
                    dag_process = subprocess.Popen(
                        ["snakemake", "-s", snakefile, "--configfile", config_file, "--dag", "--forceall"],
                        stdout=subprocess.PIPE
                    )
                    subprocess.run(
                        ["dot", "-Tpdf"],
                        stdin=dag_process.stdout,
                        stdout=out_pdf,
                        check=True
                    )
                logger.info(f"DAG exported successfully to {output_file}")
            except FileNotFoundError:
                logger.error("Graphviz 'dot' not found. Install it with e.g. `conda install -c conda-forge graphviz`.")
                sys.exit(1)
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to export DAG: {e}")
                sys.exit(1)
            return  # Do not run workflow
        else:
            cmd = ["snakemake", "-s", snakefile, "--configfile", config_file]
            if dry_run:
                cmd.append("-n")
            cmd += ["-j", str(cores or 1)]
            if additional_args:
                cmd += additional_args.split()

            logger.info(f"Running command: {' '.join(cmd)}")

            try:
                subprocess.run(cmd, check=True)
                if not dry_run:
                    logger.info("Workflow completed successfully!")
            except subprocess.CalledProcessError as e:
                logger.error(f"Snakemake failed: {e}")
                sys.exit(1)
            except FileNotFoundError:
                logger.error("Snakemake not found. Please install it.")
                sys.exit(1)


