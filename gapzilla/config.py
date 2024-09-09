"""
This submodule handles parsing data from config file and setting up logging.
"""

import logging
import yaml
import multiprocessing as mp


# Constants
OUTPUT_FILE_NAME = ""
"""Output file name. If not specified, will save input file name and add suffix"""
PATH_TO_OUTPUT_FOLDER = ""
"""Path to the gbk output folder.  If not specified, will use input folder"""
SUFFIX_FILE_NAME = "output"
"""Specify filename suffix to add"""

MIN_GAP_LENGTH = 0
"""Minimum gap length"""
MAX_GAP_LENGTH = 150
"""Maximum gap length"""

MIN_FLANKS_LENGTH = 700
"""Minimum flanking CDS length"""
MAX_FLANKS_LENGTH = 10000000
"""Maximum flanking CDS length"""

MFE_THRESHOLD_HPA = -7
"""MFE (minimal free energy) threshold for RNA hairpins"""
MFE_THRESHOLD_HPT = -15
"""MFE threshold to filter most stable RNA hairpins"""

BORDER_SHIFT = 75
"""Distance to the right and to the left from gap to search for hairpins"""

HAIRPIN_SIMILARITY_THRES = 0.9
"""Threshold for filtering similar hairpins. Lower thres -> more hairpins dropped"""

AVOID_PLOTTING = []
"""Specify instances that will not be plotted. Example to drop hpa: -ap all_hairpins_f all_hairpins_r. Full list of plottable instances: gaps, top_hairpins_f, top_hairpins_r, all_hairpins_f, all_hairpins_r, insertion_sites"""


NUM_PROCESSES = mp.cpu_count()
"""Maximum number of threads to use. Takes all available cores, if not specified"""

VERBOSITY = 1
"""
Specify verbosity of output (0 - silent, 1 - max)
"""

BAR_FORMAT = "{desc:<50.100} {percentage:3.0f}% |{bar:100}{r_bar}"
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
