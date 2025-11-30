import os
import sys
import yaml
import shlex
import subprocess

from pathlib import Path
from cellranger_snakemake.utils.custom_logger import custom_logger
from cellranger_snakemake.config_templates import ARC_CONFIG, ARC_README_content, ATAC_CONFIG, ATAC_README_content, GEX_CONFIG, GEX_README_content

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
            custom_logger.error(f"Unknown workflow type: {self.name}")
            sys.exit(1)

    def _resolve_snakefile(self):
        snakefile_name = f"cellranger{self.name}.smk"
        path = Path(__file__).parent / self.name / snakefile_name
        if not path.exists():
            custom_logger.error(f"Could not find Snakefile at expected path: {path}")
            sys.exit(1)
        return str(path)

    def sanity_check(self):
        """Validate default config filename, if present, and other preconditions."""
        if self.default_config_filename:
            if not isinstance(self.default_config_filename, str):
                custom_logger.error("`--get-default-config` must be a string filename or True.")
                sys.exit(1)
            if not self.default_config_filename.endswith('.yaml'):
                custom_logger.error("`--get-default-config` must be a filename ending in '.yaml'")
                sys.exit(1)
            if os.path.exists(self.default_config_filename):
                custom_logger.error(f"File '{self.default_config_filename}' already exists and will not be overwritten.")
                sys.exit(1)


    def write_default_config(self):
        """Write default config only if requested."""
        if self.default_config_filename:
            custom_logger.info(f"Writing '{self.name}' default config YAML to '{self.default_config_filename}'")
            with open(self.default_config_filename, 'w') as f:
                yaml.dump(self.default_config, f, indent=2, sort_keys=False)
    

    def generate_config_readme(self):
        """Generate a README file with detailed configuration instructions."""
        if self.default_config_filename:
            if self.name == "ARC":
                readme_content = ARC_README_content.format(workflow_type=self.name, config_filename=self.default_config_filename)
            elif self.name == "ATAC":
                readme_content = ATAC_README_content.format(workflow_type=self.name, config_filename=self.default_config_filename)
            elif self.name == "GEX":
                readme_content = GEX_README_content.format(workflow_type=self.name, config_filename=self.default_config_filename)
            else:
                custom_logger.warning(f"Unknown workflow name '{self.name}'. Cannot generate README.")
                return None

            readme_filename = f"{self.name}_CONFIG_README.md"
            custom_logger.info(f"Generating configuration README for '{self.name}' workflow: {readme_filename}")
            with open(readme_filename, 'w') as f:
                f.write(readme_content)

            return readme_filename
        else:
            custom_logger.warning("Skipping README generation.")
            return None
    

    def run(self, config_file=None, snakefile=None, dry_run=False, cores=1, dag=False, additional_args=""):
        config_file = config_file or self.config_file
        snakefile = snakefile or self.snakefile

        if not os.path.exists(config_file):
            custom_logger.error(f"Configuration file '{config_file}' not found.")
            sys.exit(1)

        if not os.path.exists(snakefile):
            custom_logger.error(f"Snakefile '{snakefile}' not found.")
            sys.exit(1)

        if dag:
            # Export DAG to PDF instead of running workflow
            output_file = f"{self.name}_dag.pdf"
            custom_logger.info(f"Exporting DAG to {output_file}")
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
                custom_logger.info(f"DAG exported successfully to {output_file}")
            except FileNotFoundError:
                custom_logger.error("Graphviz 'dot' not found. Install it with e.g. `conda install -c conda-forge graphviz`.")
                sys.exit(1)
            except subprocess.CalledProcessError as e:
                custom_logger.error(f"Failed to export DAG: {e}")
                sys.exit(1)
            return  # Do not run workflow
        else:
            cmd = ["snakemake", "-s", snakefile, "--configfile", config_file]
            if dry_run:
                cmd.append("-n")
            cmd += ["-j", str(cores or 1)]
            if additional_args:
                custom_logger.info(f"Additional args: {additional_args.strip()}")
                cmd += shlex.split(additional_args.strip())

            custom_logger.info(f"Running command: {' '.join(cmd)}")

            try:
                subprocess.run(cmd, check=True)
                if not dry_run:
                    custom_logger.info("Workflow completed successfully!")
            except subprocess.CalledProcessError as e:
                custom_logger.error(f"Snakemake failed: {e}")
                sys.exit(1)
            except FileNotFoundError:
                custom_logger.error("Snakemake not found. Please install it.")
                sys.exit(1)


