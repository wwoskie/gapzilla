import logging
import yaml
import multiprocessing as mp


with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Constants
OUTPUT_FILE_NAME = config["output_file_name"]
PATH_TO_OUTPUT_FOLDER = config["path_to_output_folder"]
SUFFIX_FILE_NAME = config["suffix_file_name"]

OVERLAP_SIZE = config["overlap_size"]

MIN_GAP_LENGTH = config["min_gap_length"]
MAX_GAP_LENGTH = config["max_gap_length"]
MIN_FLANKS_LENGTH = config["min_flanks_length"]
MAX_FLANKS_LENGTH = config["max_flanks_length"]
MFE_THRESHOLD_HPT = config["mfe_threshold_hpt"]
MFE_THRESHOLD_HPA = config["mfe_threshold_hpa"]

BAR_FORMAT = config["bar_format"]

NUM_PROCESSES = (
    config["num_processes"] or mp.cpu_count()
)  # take maximum available if not specified

VERBOSITY = config["verbosity"]


def setup_logging(verbosity):
    log_format = "%(message)s"
    if verbosity == 0:
        level = logging.ERROR
    elif verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    else:
        level = logging.WARNING  # Default
    logging.basicConfig(level=level, format=log_format)
