import argparse
import logging
import shutil
import multiprocessing as mp

from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

from gapzilla.config import setup_logging, config
from gapzilla.file_handling import create_output_path, modify_first_line
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
from gapzilla.utils import timeit


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

    # Parse all args
    args = parser.parse_args()

    path_to_gbk = args.path_to_gbk

    output_file_name = args.output_file_name
    path_to_output_folder = args.path_to_output_folder
    suffix_file_name = args.suffix_file_name

    path_to_output = create_output_path(
        path_to_gbk, output_file_name, path_to_output_folder, suffix_file_name
    )

    min_gap_length = args.min_gap_length
    max_gap_length = args.max_gap_length

    min_flanks_length = args.min_flanks_length
    max_flanks_length = args.max_flanks_length

    num_processes = args.threads

    mfe_threshold_hpa = args.mfe_threshold_hpa
    mfe_threshold_hpt = args.mfe_threshold_hpt
    verbosity = args.verbosity

    setup_logging(args.verbosity)

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
        subsequences, mfe_threshold_hpt, mfe_threshold_hpa, num_processes
    )

    logging.info("Processing reverse strand...")

    subsequences = split_sequence(str(sequence.complement_rna()), filtered_intervals)

    top_hairpins_r, all_hairpins_r = find_hairpins_in_subseqs(
        subsequences, mfe_threshold_hpt, mfe_threshold_hpa, num_processes
    )

    logging.info("Finding insertion sites...")
    insertion_sites_f = find_insertion_sites(filtered_intervals, top_hairpins_f)
    insertion_sites_r = find_insertion_sites(filtered_intervals, top_hairpins_r)

    insertion_sites_overlapping = find_overlapping_insertion_sites(
        insertion_sites_f + insertion_sites_r
    )

    logging.info("Writing output...")

    record.features = (
        record.features
        + create_uncovered_intervals_feature(filtered_intervals)
        + create_hairpins_feature(top_hairpins_f, "hpt")
        + create_hairpins_feature(all_hairpins_f, "hpa")
        + create_hairpins_feature(top_hairpins_r, "hpt", complement=True)
        + create_hairpins_feature(all_hairpins_r, "hpa", complement=True)
        + create_insertion_sites_feature(insertion_sites_overlapping)
    )

    with open(path_to_output, "w") as handle:
        SeqIO.write(record, handle, "genbank")

    logging.info("Finished!")


if __name__ == "__main__":
    main()
