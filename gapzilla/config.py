import logging
import yaml
import multiprocessing as mp


with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)


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
