import argparse
import logging
import shutil
import multiprocessing as mp

import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

from gapzilla.config import setup_logging, config
from gapzilla.file_handling import create_output_path, modify_first_line, uniquify_path
from gapzilla.feature_analysis import (
    create_hairpins_feature,
    create_insertion_sites_feature,
    create_uncovered_intervals_feature,
)
from gapzilla.hairpin_analysis import find_hairpins_in_subseqs
from gapzilla.models import IntervaledFeature
from gapzilla.sequence_processing import (
    merge_intervals,
    find_uncovered_intervals,
    filter_intervals_by_length,
    filter_intervals_by_flanking_legth,
    split_sequence,
)
from gapzilla.insertion_finder import (
    find_insertion_sites,
    find_overlapping_insertion_sites,
)
from gapzilla.utils import timeit, list_of_strings


# Constants
OUTPUT_FILE_NAME = config["output_file_name"]
PATH_TO_OUTPUT_FOLDER = config["path_to_output_folder"]
SUFFIX_FILE_NAME = config["suffix_file_name"]

MIN_GAP_LENGTH = config["min_gap_length"]
MAX_GAP_LENGTH = config["max_gap_length"]
MIN_FLANKS_LENGTH = config["min_flanks_length"]
MAX_FLANKS_LENGTH = config["max_flanks_length"]
MFE_THRESHOLD_HPT = config["mfe_threshold_hpt"]
MFE_THRESHOLD_HPA = config["mfe_threshold_hpa"]

HAIRPIN_SIMILARITY_THRES = config["hairpin_similaruty_thres"]

BAR_FORMAT = config["bar_format"]

NUM_PROCESSES = (
    config["num_processes"] or mp.cpu_count()
)  # take maximum available if not specified

VERBOSITY = config["verbosity"]

if config["avoid_plotting"]:
    AVOID_PLOTTING = config["avoid_plotting"].split(" ")
else:
    AVOID_PLOTTING = []


@timeit
def main():
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
        help="Minimum flanks length",
    )
    parser.add_argument(
        "--max_flanks_length",
        "-max_flanks",
        type=int,
        default=MAX_FLANKS_LENGTH,
        help="Maximum flanks length",
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
        type=int,
        default=HAIRPIN_SIMILARITY_THRES,
        help="Threshold for filtering similar hairpins. Lower thres -> more hairpins dropped",
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

    if custom_config is not None:
        with open(custom_config, "r") as file:
            config = yaml.safe_load(file)

        output_file_name = config["output_file_name"]
        path_to_output_folder = config["path_to_output_folder"]
        suffix_file_name = config["suffix_file_name"]

        min_gap_length = config["min_gap_length"]
        max_gap_length = config["max_gap_length"]

        min_flanks_length = config["min_flanks_length"]
        max_flanks_length = config["max_flanks_length"]

        mfe_threshold_hpt = config["mfe_threshold_hpt"]
        mfe_threshold_hpa = config["mfe_threshold_hpa"]

        hairpin_similarity_thres = config["hairpin_similarity_thres"]

        num_processes = config["num_processes"] or mp.cpu_count()

        verbosity = config["verbosity"]
        setup_logging(args.verbosity)

        if config["avoid_plotting"]:
            avoid_plotting = config["avoid_plotting"].split(" ")
        else:
            avoid_plotting = []
    else:
        output_file_name = args.output_file_name
        path_to_output_folder = args.path_to_output_folder
        suffix_file_name = args.suffix_file_name

        min_gap_length = args.min_gap_length
        max_gap_length = args.max_gap_length

        min_flanks_length = args.min_flanks_length
        max_flanks_length = args.max_flanks_length

        num_processes = args.threads

        mfe_threshold_hpa = args.mfe_threshold_hpa
        mfe_threshold_hpt = args.mfe_threshold_hpt

        hairpin_similarity_thres = args.hairpin_similarity_thres

        verbosity = args.verbosity
        setup_logging(args.verbosity)

        avoid_plotting = args.avoid_plotting

    path_to_output = uniquify_path(
        create_output_path(
            path_to_gbk, output_file_name, path_to_output_folder, suffix_file_name
        )
    )

    logging.info(f"Writing output to: {path_to_output}")

    feature_list = []
    intervaled_features_list = []

    # Try opening file as proper gbk
    logging.info("Reading gbk...")
    try:
        record = SeqIO.read(
            path_to_gbk,
            "genbank",
        )

        shutil.copy(path_to_gbk, path_to_output)

    # Try to fix bp LOCUS string issue
    except ValueError:
        logging.info("Trying to fix possible prokka LOCUS input error...")
        modify_first_line(path_to_gbk, path_to_output)
        record = SeqIO.read(
            path_to_output,
            "genbank",
        )

    for feature in record.features:
        feature_list.append(feature)

    for feature in tqdm(
        feature_list,
        desc="Processing annotations...",
        disable=logging.root.level > logging.INFO,
        bar_format=BAR_FORMAT,
    ):
        start, stop = int(feature.location.start), int(feature.location.end)
        feature.length = stop - start

        intervaled_features_list.append(IntervaledFeature([start, stop], feature))

    merged_intervals = merge_intervals(intervaled_features_list[1:])
    uncovered_intervals = find_uncovered_intervals(
        intervaled_features_list[0].interval, merged_intervals
    )

    filtered_intervals = filter_intervals_by_length(
        uncovered_intervals, min_gap_length, max_gap_length
    )

    filtered_intervals = filter_intervals_by_flanking_legth(
        filtered_intervals, min_flanks_length, max_flanks_length
    )

    sequence = Seq(record.seq)

    logging.info("Annotating RNA hairpins...")
    logging.info("This might take a while...")

    subsequences = split_sequence(str(sequence.transcribe()), filtered_intervals)

    logging.info(f"Total # of gaps: {len(subsequences)}")
    logging.info("Processing forward strand...")

    top_hairpins_f, all_hairpins_f = find_hairpins_in_subseqs(
        subsequences,
        mfe_threshold_hpt,
        mfe_threshold_hpa,
        num_processes,
        hairpin_similarity_thres,
    )

    logging.info("Processing reverse strand...")

    subsequences = split_sequence(str(sequence.complement_rna()), filtered_intervals)

    top_hairpins_r, all_hairpins_r = find_hairpins_in_subseqs(
        subsequences,
        mfe_threshold_hpt,
        mfe_threshold_hpa,
        num_processes,
        hairpin_similarity_thres,
    )

    logging.info("Finding insertion sites...")
    insertion_sites_f = find_insertion_sites(filtered_intervals, top_hairpins_f)
    insertion_sites_r = find_insertion_sites(filtered_intervals, top_hairpins_r)

    insertion_sites_overlapping = find_overlapping_insertion_sites(
        insertion_sites_f + insertion_sites_r
    )

    logging.info("Writing output...")

    print(avoid_plotting)

    plot_dict = {
        "gaps": create_uncovered_intervals_feature(filtered_intervals),
        "top_hairpins_f": create_hairpins_feature(top_hairpins_f, "hpt"),
        "all_hairpins_f": create_hairpins_feature(all_hairpins_f, "hpa"),
        "top_hairpins_r": create_hairpins_feature(
            top_hairpins_r, "hpt", complement=True
        ),
        "all_hairpins_r": create_hairpins_feature(
            all_hairpins_r, "hpa", complement=True
        ),
        "insertion_sites": create_insertion_sites_feature(insertion_sites_overlapping),
    }

    list_to_plot = []
    for key in plot_dict.keys():
        if key not in avoid_plotting:
            list_to_plot += plot_dict[key]

    record.features = record.features + list_to_plot

    with open(path_to_output, "w") as handle:
        SeqIO.write(record, handle, "genbank")

    logging.info("Finished!")


if __name__ == "__main__":
    main()
