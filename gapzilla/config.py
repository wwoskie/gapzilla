"""
This submodule handles parsing data from config file and setting up logging.
"""

import logging
import yaml
import multiprocessing as mp


with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)


# Constants
OUTPUT_FILE_NAME = config["output_file_name"]
"""Output file name. If not specified, will save input file name and add suffix"""
PATH_TO_OUTPUT_FOLDER = config["path_to_output_folder"]
"""Path to the gbk output folder"""
SUFFIX_FILE_NAME = config["suffix_file_name"]
"""Specify filename suffix to add"""

MIN_GAP_LENGTH = config["min_gap_length"]
"""Minimum gap length"""
MAX_GAP_LENGTH = config["max_gap_length"]
"""Maximum gap length"""

MIN_FLANKS_LENGTH = config["min_flanks_length"]
"""Minimum flanking CDS length"""
MAX_FLANKS_LENGTH = config["max_flanks_length"]
"""Maximum flanking CDS length"""

MFE_THRESHOLD_HPA = config["mfe_threshold_hpa"]
"""MFE (minimal free energy) threshold for RNA hairpins"""
MFE_THRESHOLD_HPT = config["mfe_threshold_hpt"]
"""MFE threshold to filter most stable RNA hairpins"""

BORDER_SHIFT = config["border_shift"]
"""Distance to the right and to the left from gap to search for hairpins"""

HAIRPIN_SIMILARITY_THRES = config["hairpin_similarity_thres"]
"""Threshold for filtering similar hairpins. Lower thres -> more hairpins dropped"""

AVOID_PLOTTING = None
"""Specify instances that will not be plotted. Example to drop hpa: -ap all_hairpins_f all_hairpins_r. Full list of plottable instances: gaps, top_hairpins_f, top_hairpins_r, all_hairpins_f, all_hairpins_r, insertion_sites"""

if config["avoid_plotting"]:
    AVOID_PLOTTING = config["avoid_plotting"].split(" ")
else:
    AVOID_PLOTTING = []


NUM_PROCESSES = config["num_processes"] or mp.cpu_count()
"""Maximum number of threads to use. Takes all available cores, if not specified"""

VERBOSITY = config["verbosity"]
"""
Specify verbosity of output (0 - silent, 1 - max)
"""

BAR_FORMAT = config["bar_format"]
"""
Format of the progress bar for tqdm module.
"""


def setup_logging(verbosity: int) -> None:
    """
    Configures the logging settings based on the provided verbosity level.

    Parameters
    ----------
    verbosity : int
        The verbosity level for logging. The levels are defined as:
        - 0: ERROR level, only error messages will be logged.
        - 1: INFO level, informational messages and above will be logged.
        - 2 or higher: DEBUG level, all messages including debugging details will be logged.

    Returns
    -------
    None

    Notes
    -----
    This function sets up the logging configuration with a simple message format and adjusts the logging
    level according to the verbosity parameter. If an invalid verbosity level is provided, it defaults to
    WARNING level.
    """

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
