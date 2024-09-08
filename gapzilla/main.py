"""
This module handles CLI behaviour of the module. It takes constants from config file and processes parameters from command line arguments or custom config file and then invokes `process_gbk` function from gbk_processing.py
"""

import argparse
import multiprocessing as mp

import yaml

from gapzilla.config import *
from gapzilla.gbk_processing import process_gbk
from gapzilla.utils import timeit, list_of_strings


@timeit
def main() -> None:
    """
    Parses command line input and invokes `process_gbk()` function.
    """

    parser = argparse.ArgumentParser(
        description="Find RNA hairpin-flanked gaps in GBK file for potential insert sites",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("path_to_gbk", help="Path to the gbk input file")

    parser.add_argument(
        "--output_file_name",
        "-f",
        default=OUTPUT_FILE_NAME,
        help="Output file name. If not specified, will save input file name and add suffix",
    )

    parser.add_argument(
        "--suffix_file_name",
        "-s",
        default=SUFFIX_FILE_NAME,
        help="Specify filename suffix to add",
    )

    parser.add_argument(
        "--path_to_output_folder",
        "-o",
        default=PATH_TO_OUTPUT_FOLDER,
        help="Path to the gbk output folder",
    )

    parser.add_argument(
        "--min_gap_length",
        "-min_gap",
        type=int,
        default=MIN_GAP_LENGTH,
        help="Minimum gap length",
    )

    parser.add_argument(
        "--max_gap_length",
        "-max_gap",
        type=int,
        default=MAX_GAP_LENGTH,
        help="Maximum gap length",
    )

    parser.add_argument(
        "--min_flanks_length",
        "-min_flanks",
        type=int,
        default=MIN_FLANKS_LENGTH,
        help="Minimum flanking CDS length",
    )
    parser.add_argument(
        "--max_flanks_length",
        "-max_flanks",
        type=int,
        default=MAX_FLANKS_LENGTH,
        help="Maximum flanking CDS length",
    )

    parser.add_argument(
        "--mfe_threshold_hpa",
        "-hpa",
        type=int,
        default=MFE_THRESHOLD_HPA,
        help="MFE (minimal free energy) threshold for RNA hairpins",
    )

    parser.add_argument(
        "--mfe_threshold_hpt",
        "-hpt",
        type=int,
        default=MFE_THRESHOLD_HPT,
        help="MFE threshold to filter most stable RNA hairpins",
    )

    parser.add_argument(
        "--hairpin_similarity_thres",
        "-similarity",
        type=float,
        default=HAIRPIN_SIMILARITY_THRES,
        help="Threshold for filtering similar hairpins. Lower thres -> more hairpins dropped",
    )

    parser.add_argument(
        "--border_shift",
        "-bs",
        type=int,
        default=BORDER_SHIFT,
        help="Distance to the right and to the left from gap to search for haiirpins",
    )

    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=NUM_PROCESSES,
        help="Maximum number of threads to use. Takes all available cores, if not specified",
    )

    parser.add_argument(
        "--verbosity",
        "-v",
        type=int,
        default=VERBOSITY,
        help="Specify verbosity of output (0 - silent, 1 - max)",
    )

    parser.add_argument(
        "--avoid_plotting",
        "-ap",
        default=AVOID_PLOTTING,
        metavar="N",
        type=str,
        nargs="*",
        help="Specify instances that will not be plotted. Example to drop hpa: -ap all_hairpins_f all_hairpins_r. Full list of plottable instances: gaps, top_hairpins_f, top_hairpins_r, all_hairpins_f, all_hairpins_r, insertion_sites",
    )

    parser.add_argument(
        "--custom_config",
        "-cc",
        default=None,
        metavar="N",
        type=str,
        help="Path to custom yaml config file that will override defaults and any other commandline args",
    )

    # Parse all args
    args = parser.parse_args()
    path_to_gbk = args.path_to_gbk
    custom_config = args.custom_config

    # Overrride all if custom config is provided

    if custom_config is not None:
        with open(custom_config, "r") as file:
            custom_config = yaml.safe_load(file)
        if custom_config["avoid_plotting"]:
            custom_config["avoid_plotting"] = custom_config["avoid_plotting"].split(" ")

    else:
        custom_config = {}

    # Set all necessary variable, if available, take from custom config, else use from args (defaults if other not provided)
    output_file_name = custom_config.get("output_file_name", args.output_file_name)
    path_to_output_folder = custom_config.get(
        "path_to_output_folder", args.path_to_output_folder
    )
    suffix_file_name = custom_config.get("suffix_file_name", args.suffix_file_name)

    min_gap_length = custom_config.get("min_gap_length", args.min_gap_length)
    max_gap_length = custom_config.get("max_gap_length", args.max_gap_length)

    min_flanks_length = custom_config.get("min_flanks_length", args.min_flanks_length)
    max_flanks_length = custom_config.get("max_flanks_length", args.max_flanks_length)

    mfe_threshold_hpa = custom_config.get("mfe_threshold_hpt", args.mfe_threshold_hpa)
    mfe_threshold_hpt = custom_config.get("mfe_threshold_hpa", args.mfe_threshold_hpt)

    hairpin_similarity_thres = custom_config.get(
        "hairpin_similarity_thres", args.hairpin_similarity_thres
    )
    border_shift = custom_config.get("border_shift", args.border_shift)

    num_processes = custom_config.get("num_processes", args.threads)

    verbosity = custom_config.get("verbosity", args.verbosity)

    avoid_plotting = custom_config.get("avoid_plotting", args.avoid_plotting)

    process_gbk(
        path_to_gbk,
        output_file_name,
        path_to_output_folder,
        suffix_file_name,
        min_gap_length,
        max_gap_length,
        min_flanks_length,
        max_flanks_length,
        mfe_threshold_hpa,
        mfe_threshold_hpt,
        hairpin_similarity_thres,
        border_shift,
        num_processes,
        verbosity,
        avoid_plotting,
    )


if __name__ == "__main__":
    main()
