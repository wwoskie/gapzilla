"""
Module to process GenBank (gbk) file to identify and annotate RNA hairpins and insertion sites.
"""

import logging
import os

from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

from gapzilla.config import setup_logging, BAR_FORMAT
from gapzilla.file_processing import (
    create_output_path,
    modify_first_line,
    uniquify_path,
)
from gapzilla.feature_processing import (
    create_hairpins_feature,
    create_insertion_sites_feature,
    create_uncovered_intervals_feature,
)
from gapzilla.hairpin_processing import find_hairpins_in_subseqs
from gapzilla.models import IntervaledFeature
from gapzilla.interval_processing import (
    merge_intervals,
    find_uncovered_intervals,
    filter_intervals_by_length,
    filter_intervals_by_flanking_legth,
    filter_intervals_by_strand_direction,
    split_sequence,
)
from gapzilla.insertion_processing import (
    find_insertion_sites,
    find_overlapping_insertion_sites,
    filter_insertion_sites_by_max_score,
    filter_insertion_sites_by_hairpins,
)


def process_gbk(
    path_to_gbk: str | os.PathLike,
    output_file_name: str = "",
    path_to_output_folder: str | os.PathLike = "",
    suffix_file_name: str = "output",
    min_gap_length: int = 0,
    max_gap_length: int = 150,
    min_flanks_length: int = 700,
    max_flanks_length: int = 1e6,
    mfe_threshold_hpa: float = -7,
    mfe_threshold_hpt: float = -15,
    hairpin_similarity_thres: float = 0.9,
    border_shift: int = 75,
    num_processes: int = 2,
    verbosity: int = 0,
    avoid_plotting: list[str] | None = [],
) -> None:
    """
    Process a GenBank (gbk) file to identify and annotate RNA hairpins and insertion sites.

    This function reads a GenBank file with multiple records, identifies gaps between annotated features,
    filters these gaps based on provided criteria, and then searches for RNA hairpins
    within these gaps and in nearby regions according to the provided border shift. The identified hairpins and insertion sites are then annotated
    back into the GenBank file.

    Parameters
    ----------
    path_to_gbk : str or os.PathLike
        Path to the input GenBank file.
    output_file_name : str, optional
        Name of the output file. If not provided, a default name will be generated.
    path_to_output_folder : str or os.PathLike, optional
        Path to the folder where the output file will be saved. If not provided, the
        output file will be saved in the same directory as the input file.
    suffix_file_name : str, optional
        Suffix to append to the output file name. Default is "output".
    min_gap_length : int, optional
        Minimum length of gaps to consider. Default is 0.
    max_gap_length : int, optional
        Maximum length of gaps to consider. Default is 150.
    min_flanks_length : int, optional
        Minimum length of flanking regions around gaps. Default is 700.
    max_flanks_length : int, optional
        Maximum length of flanking regions around gaps. Default is 1e6.
    mfe_threshold_hpa : float, optional
        Minimum free energy threshold for hairpin acceptance. Default is -7.
    mfe_threshold_hpt : float, optional
        Minimum free energy threshold for hairpin truncation. Default is -15.
    hairpin_similarity_thres : float, optional
        Similarity threshold for hairpin filtering. Default is 0.9.
    border_shift : int, optional
        Number of nucleotides to shift borders when splitting sequences. Default is 75.
    num_processes : int, optional
        Number of processes to use for parallel computation. Default is 2.
    verbosity : int, optional
        Verbosity level for logging. Default is 0.
    avoid_plotting : list of str or None, optional
        List of feature types to avoid plotting. Default is an empty list.

    Returns
    -------
    None
        The function writes the annotated GenBank file to the specified output path.
    """

    setup_logging(verbosity)

    path_to_output = uniquify_path(
        create_output_path(
            path_to_gbk, output_file_name, path_to_output_folder, suffix_file_name
        )
    )

    logging.info(f"Writing output to: {path_to_output}")

    # List to read and append records to
    records_list = []
    records_output = []

    # Try opening file as proper gbk
    logging.info("Reading gbk...")
    try:
        with open(path_to_gbk) as handle:
            for record in SeqIO.parse(handle, format="gb"):
                records_list.append(record)

    # Try to fix bp LOCUS string issue
    except ValueError:
        logging.info("Trying to fix possible prokka LOCUS input error...")
        modify_first_line(path_to_gbk, path_to_output)

        with open(path_to_output) as handle:
            for record in SeqIO.parse(handle, format="gb"):
                records_list.append(record)

    number_of_records = len(records_list)

    for indx, record in enumerate(records_list):
        if indx == 0 and number_of_records > 1:
            logging.info(f"Total number of records in file: {number_of_records}")
        if number_of_records > 1:
            logging.info(f"Processing record # {indx + 1}: {record.name}")

        # Parsing features
        intervaled_features_list = []

        for feature in tqdm(
            record.features,
            desc="Processing annotations...",
            disable=logging.root.level > logging.INFO,
            bar_format=BAR_FORMAT,
        ):
            start, stop = int(feature.location.start), int(feature.location.end)
            feature.length = stop - start

            intervaled_features_list.append(IntervaledFeature(start, stop, feature))

        merged_intervals = merge_intervals(
            intervaled_features_list[1:]
        )  # Merging overlapping features
        uncovered_intervals = find_uncovered_intervals(
            intervaled_features_list[0].start,
            intervaled_features_list[0].end,
            merged_intervals,
        )  # Finding uncovered intervals (gaps)

        # Filtering gaps
        filtered_intervals = filter_intervals_by_length(
            uncovered_intervals, min_gap_length, max_gap_length
        )

        filtered_intervals = filter_intervals_by_flanking_legth(
            filtered_intervals, min_flanks_length, max_flanks_length
        )

        filtered_intervals = filter_intervals_by_strand_direction(filtered_intervals)

        sequence = Seq(record.seq)

        logging.info("Annotating RNA hairpins...")
        subsequences = split_sequence(
            str(sequence.transcribe()), filtered_intervals, border_shift
        )
        # Finding hairpins
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
        subsequences = split_sequence(
            str(sequence.complement_rna()), filtered_intervals, border_shift
        )
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
        insertion_sites_overlapping = filter_insertion_sites_by_max_score(
            insertion_sites_overlapping
        )
        insertion_sites_overlapping = filter_insertion_sites_by_hairpins(
            insertion_sites_overlapping, all_hairpins_f + all_hairpins_r
        )

        # Handle what to plot
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
            "insertion_sites": create_insertion_sites_feature(
                insertion_sites_overlapping
            ),
        }

        # Feature will be plotted only if it's not in avoid_plotting list
        list_to_plot = [
            feature
            for key, features in plot_dict.items()
            if key not in avoid_plotting
            for feature in features
        ]

        record.features = record.features + list_to_plot

        records_output.append(record)

    logging.info("Writing output...")
    with open(path_to_output, "w") as handle:
        SeqIO.write(records_list, handle, "genbank")

    logging.info("Finished!")
