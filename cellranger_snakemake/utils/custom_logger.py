# custom_logger.py
import logging

class ColorFormatter(logging.Formatter):
    COLORS = {
        'DEBUG': "\033[37m",
        'INFO': "\033[36m",
        'WARNING': "\033[33m",
        'ERROR': "\033[31m",
        'CRITICAL': "\033[41m",
    }
    RESET = "\033[0m"

    def format(self, record):
        color = self.COLORS.get(record.levelname, self.RESET)
        message = super().format(record)
        return f"{color}{message}{self.RESET}"

# Initialize custom_logger
custom_logger = logging.getLogger("cellranger_snakemake")
custom_logger.setLevel(logging.DEBUG)  # Adjust level as needed

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

formatter = ColorFormatter("[%(levelname)s] %(message)s")
ch.setFormatter(formatter)
custom_logger.addHandler(ch)
