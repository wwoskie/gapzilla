import argparse
import datetime
import os
import re
import shutil
import time
import multiprocessing as mp

import RNA
import yaml
from Bio import GenBank
from Bio.Seq import Seq
from functools import wraps, reduce
from tqdm import tqdm

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Constants
OUTPUT_FILE_NAME = config["output_file_name"]
PATH_TO_OUTPUT_FOLDER = config["path_to_output_folder"]
SUFFIX_FILE_NAME = config["suffix_file_name"]

BATCH_SIZE = config["batch_size"]
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


class IntervaledFeature:
    def __init__(self, interval, feature):
        self.interval = interval
        self.feature_list = [feature]
        self.feature_num = len(self.feature_list)
        self.feature_lengths = [feature.length]

    def __repr__(self):
        return f"IntervaledFeature(interval={self.interval}, feature_list={self.feature_list}, feature_num={self.feature_num}, feature_lengths={self.feature_lengths})"


class IntervaledGap:
    def __init__(self, interval, features_left, features_right):
        self.interval = interval
        self.features_left = features_left
        self.features_right = features_right

    def __repr__(self):
        return f"IntervaledGap(interval={self.interval}, features_left={self.features_left}, features_right={self.features_right})"


class InsertionSite:
    def __init__(self, start, end, score):
        self.interval = [start, end]
        self.score = score

    def __repr__(self):
        return f"InsertionSite(interval={self.interval}, score={self.score})"


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = datetime.datetime.now()
        result = func(*args, **kwargs)
        end_time = datetime.datetime.now()
        total_time = str(end_time - start_time).split(".")[0]
        print(f"Total execution time {total_time}")
        return result

    return timeit_wrapper


def create_dirs_to_output(path_to_dirs):
    if not os.path.isdir(path_to_dirs):
        os.makedirs(path_to_dirs)


def create_output_path(
    input_path, output_file_name, path_to_output_folder, suffix_file_name="output"
):
    if path_to_output_folder:
        create_dirs_to_output(path_to_output_folder)
    if not output_file_name and not path_to_output_folder:
        directory, base_name = os.path.split(input_path)
        file_name, file_extension = os.path.splitext(base_name)
        new_file_name = f"{file_name}_{suffix_file_name}{file_extension}"

    elif not output_file_name and path_to_output_folder:
        directory, base_name = os.path.split(input_path)
        directory = os.path.join(directory, path_to_output_folder)
        file_name, file_extension = os.path.splitext(base_name)
        new_file_name = f"{file_name}_{suffix_file_name}{file_extension}"

    elif output_file_name and not path_to_output_folder:
        directory = os.path.split(input_path)[0]
        new_file_name = output_file_name

    elif output_file_name and path_to_output_folder:
        directory = path_to_output_folder
        new_file_name = output_file_name

    output_path = os.path.join(directory, new_file_name)

    return output_path


def modify_first_line(input_path, output_path):
    with open(input_path, "r") as file:
        lines = file.readlines()

    if lines and "LOCUS" in lines[0]:
        first_line = lines[0]
        bp_index = first_line.find("bp")
        if bp_index != -1:
            modified_first_line = first_line[:bp_index] + "0 " + first_line[bp_index:]
            lines[0] = modified_first_line

    with open(output_path, "w") as file:
        file.writelines(lines)


def insert_before_origin(file_path, insert_strings):
    with open(file_path, "r") as file:
        content = file.read()

    origin_index = content.find("ORIGIN")
    if origin_index == -1:
        return

    insert_content = "\n".join(insert_strings) + "\n"

    modified_content = content[:origin_index] + insert_content + content[origin_index:]

    with open(file_path, "w") as file:
        file.write(modified_content)


def split_sequence(sequence, window_size=10e6, overlap=10e3):
    step_size = int(window_size - overlap)
    subsequences = []
    for i in range(0, len(sequence), step_size):
        end_index = i + window_size
        if end_index > len(sequence):
            end_index = len(sequence)
        subsequences.append(sequence[i:end_index])
        if end_index == len(sequence):
            break
    return subsequences


def merge_intervals(intervaled_features):
    sorted_intervaled_features = sorted(
        intervaled_features, key=lambda x: x.interval[0]
    )

    merged_intervaled_features = []

    for intervaled_feature in tqdm(
        sorted_intervaled_features, desc="Merging intervals...", bar_format=BAR_FORMAT
    ):
        if (
            not merged_intervaled_features
            or merged_intervaled_features[-1].interval[1]
            < intervaled_feature.interval[0]
        ):
            merged_intervaled_features.append(intervaled_feature)
        else:
            merged_intervaled_features[-1].interval[1] = max(
                merged_intervaled_features[-1].interval[1],
                intervaled_feature.interval[1],
            )
            if (
                intervaled_feature.feature_list[0]
                not in merged_intervaled_features[-1].feature_list
            ):
                merged_intervaled_features[-1].feature_list.append(
                    intervaled_feature.feature_list[0]
                )

                merged_intervaled_features[-1].feature_num += 1
                merged_intervaled_features[-1].feature_lengths.append(
                    intervaled_feature.feature_lengths[0]
                )

    return merged_intervaled_features


def find_uncovered_intervals(interval, intervaled_features):

    start, end = interval
    intervaled_features.sort(key=lambda x: x.interval)  # Sort by intervals

    uncovered_intervals = []
    current_start = start
    prev_feature = None

    for feature in tqdm(
        intervaled_features,
        desc="Finding uncovered intervals...",
        bar_format=BAR_FORMAT,
    ):
        i_start, i_end = feature.interval
        if i_end <= start:
            prev_feature = feature
            continue
        if i_start >= end:
            break

        if i_start > current_start and (i_start - 1) != current_start:
            uncovered_interval = (current_start + 1, i_start - 1)
            features_left = (
                prev_feature
                if prev_feature and prev_feature.interval[1] <= uncovered_interval[0]
                else None
            )
            features_right = (
                feature if feature.interval[0] >= uncovered_interval[1] else None
            )
            uncovered_intervals.append(
                IntervaledGap(uncovered_interval, features_left, features_right)
            )

        current_start = max(current_start, i_end)
        prev_feature = feature

    if current_start < end:
        uncovered_interval = (current_start + 1, end - 1)
        features_left = (
            prev_feature
            if prev_feature and prev_feature.interval[1] <= uncovered_interval[0]
            else None
        )
        features_right = None
        uncovered_intervals.append(
            IntervaledGap(uncovered_interval, features_left, features_right)
        )

    return uncovered_intervals


def filter_intervals_by_length(intervals, min_length, max_length):
    filtered_intervals = []
    for intervaled_gap in tqdm(
        intervals, desc="Filtering intervals by length...", bar_format=BAR_FORMAT
    ):

        if intervaled_gap.interval[0] and intervaled_gap.interval[1]:
            start, end = intervaled_gap.interval[0], intervaled_gap.interval[1]
            length = end - start

            if min_length <= length + 1 <= max_length:
                filtered_intervals.append(intervaled_gap)

    return filtered_intervals


def filter_intervals_by_flanking_legth(intervals, min_flanks_length, max_flanks_length):
    filtered_intervals = []
    for intervaled_gap in tqdm(
        intervals,
        desc="Filtering intervals by flanking genes...",
        bar_format=BAR_FORMAT,
    ):
        if intervaled_gap:
            try:
                length_left = intervaled_gap.features_left.feature_lengths[-1]
                length_right = intervaled_gap.features_right.feature_lengths[0]

                if (
                    min_flanks_length <= length_left + 1 <= max_flanks_length
                    and min_flanks_length <= length_right + 1 <= max_flanks_length
                ):
                    filtered_intervals.append(intervaled_gap)
            except AttributeError:
                continue

    return filtered_intervals


def create_uncovered_intervals_str(intervals):

    return [
        f"     gap             {interval.interval[0]}..{interval.interval[1]}"
        for interval in intervals
    ]


def create_hairpins_str(hairpins, type_of_hairpin, complement=False):
    if complement:

        output_str_list = [
            f"     {type_of_hairpin}             complement({hairpin[0]}..{hairpin[1]})\n"
            + f"                     /sequence={hairpin[2]}\n"
            + f"                     /structure={hairpin[3]}\n"
            + f"                     /mfe={hairpin[4]}"
            for hairpin in hairpins
        ]

    else:
        output_str_list = [
            f"     {type_of_hairpin}             {hairpin[0]}..{hairpin[1]}\n"
            + f"                     /sequence={hairpin[2]}\n"
            + f"                     /structure={hairpin[3]}\n"
            + f"                     /mfe={hairpin[4]}"
            for hairpin in hairpins
        ]

    return output_str_list


def find_insertion_sites(gaps, coverages):
    uncovered_intervals = []
    coverages = [(c[0], c[1]) for c in sorted(coverages, key=lambda x: x[0])]

    for gap in gaps:
        g_start, g_end = gap.interval
        current_pos = g_start

        for c_start, c_end in coverages:
            if c_start > g_end:
                break
            if c_end < current_pos:
                continue
            if c_start > current_pos and c_start < g_end:
                uncovered_intervals.append((current_pos, min(c_start, g_end)))
            current_pos = max(current_pos, c_end)

        if current_pos < g_end:
            uncovered_intervals.append((current_pos, g_end))

    insertion_sites = []
    for u_start, u_end in uncovered_intervals:
        left_neighbour = any(c[1] == u_start for c in coverages)
        right_neighbour = any(c[0] == u_end for c in coverages)

        if left_neighbour:
            u_start += 1

        if right_neighbour:
            u_end -= 1

        if left_neighbour and right_neighbour:
            score = 2
        elif left_neighbour or right_neighbour:
            score = 1
        else:
            continue

        insertion_sites.append(InsertionSite(u_start, u_end, score))

    return insertion_sites


def find_overlapping_insertion_sites(insertion_sites):
    insertion_sites = sorted(insertion_sites, key=lambda x: x.interval[0])
    result = []
    n = len(insertion_sites)

    for i in range(n):
        current_site = insertion_sites[i]
        current_start, current_end = current_site.interval
        current_score = current_site.score

        # Always add the original interval first
        result.append(current_site)

        for j in range(i + 1, n):
            next_site = insertion_sites[j]
            next_start, next_end = next_site.interval

            if next_start > current_end:
                break  # No more overlaps possible

            if current_start <= next_end and current_end >= next_start:
                overlap_start = max(current_start, next_start)
                overlap_end = min(current_end, next_end)
                overlap_score = current_score + next_site.score
                result.append(InsertionSite(overlap_start, overlap_end, overlap_score))

    return result


def create_insertion_sites_str(insertion_sites):

    output_str_list = [
        f"     ins{insertion_site.score}            {insertion_site.interval[0]}..{insertion_site.interval[1]}"
        for insertion_site in insertion_sites
    ]

    return output_str_list


with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

BAR_FORMAT = config["bar_format"]


def find_hairpins(sequence, structure):
    hairpins = []
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")":
            start = stack.pop()
            end = i
            if end - start > 3:  # Minimum length for a hairpin loop
                hairpin_seq = sequence[start : end + 1]
                hairpin_struct = structure[start : end + 1]
                hairpins.append((start + 1, end + 1, hairpin_seq, hairpin_struct))
    return hairpins


def evaluate_hairpins(hairpins):
    evaluated_hairpins = []
    for start, end, seq, struct in hairpins:
        fc = RNA.fold_compound(seq)
        mfe = fc.eval_structure(struct)
        evaluated_hairpins.append((start, end, seq, struct, mfe))
    return evaluated_hairpins


def merge_similar_hairpins(
    hairpins,
    overlap_threshold=0.9,
    backward_search_thres=20,
    message="Merging similar hairpins...",
):
    merged_hairpins = []
    hairpins.sort(key=lambda x: (x[0], x[1]))  # Sort by start and end

    for hairpin in tqdm(hairpins, desc=message, bar_format=BAR_FORMAT):
        start, end, seq, struct, mfe = hairpin
        overlap_found = False
        for merged in merged_hairpins[-backward_search_thres:]:
            m_start, m_end, m_seq, m_struct, m_mfe = merged
            overlap = max(0, min(end, m_end) - max(start, m_start))
            if overlap / min(end - start, m_end - m_start) > overlap_threshold:
                overlap_found = True
                if mfe < m_mfe:
                    merged_hairpins.remove(merged)
                    merged_hairpins.append(hairpin)
                break
        if not overlap_found:
            merged_hairpins.append(hairpin)

    return merged_hairpins


def process_window(args):
    sequence, window_size, start = args
    window_seq = sequence[start : start + window_size]
    structure, mfe = RNA.fold(window_seq)
    hairpins = find_hairpins(window_seq, structure)
    evaluated_hairpins = evaluate_hairpins(hairpins)
    adjusted_hairpins = [
        (start + hp[0], start + hp[1], hp[2], hp[3], hp[4]) for hp in evaluated_hairpins
    ]
    return adjusted_hairpins


def find_top_hairpins(
    sequence,
    window_size=50,
    step_size=25,
    mfe_threshold_hpt=-15.0,
    mfe_threshold_hpa=-1,
    batch_num=None,
    num_processes=mp.cpu_count(),
):

    pool = mp.Pool(processes=num_processes)
    args = [
        (sequence, window_size, i)
        for i in range(0, len(sequence) - window_size + 1, step_size)
    ]
    all_hairpins = []

    with tqdm(
        total=len(args),
        desc=f"Calculating hairpins in batch #{batch_num}...",
        bar_format=BAR_FORMAT,
    ) as pbar:
        for result in pool.imap_unordered(process_window, args):
            all_hairpins.extend(result)
            pbar.update()

    pool.close()
    pool.join()

    # all_hairpins = merge_similar_hairpins(all_hairpins)
    all_hairpins = [hp for hp in all_hairpins if hp[4] <= mfe_threshold_hpa]
    filtered_hairpins = [hp for hp in all_hairpins if hp[4] <= mfe_threshold_hpt]
    all_hairpins = [
        hairpin for hairpin in all_hairpins if hairpin not in filtered_hairpins
    ]

    return filtered_hairpins, all_hairpins


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
        "--batch_size",
        "-b",
        type=int,
        default=BATCH_SIZE,
        help="Size of genome batches for hairpin calculations",
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

    # Parse all args
    args = parser.parse_args()

    path_to_gbk = args.path_to_gbk

    output_file_name = args.output_file_name
    path_to_output_folder = args.path_to_output_folder
    suffix_file_name = args.suffix_file_name

    path_to_output = create_output_path(
        path_to_gbk, output_file_name, path_to_output_folder, suffix_file_name
    )

    batch_size = args.batch_size

    min_gap_length = args.min_gap_length
    max_gap_length = args.max_gap_length

    min_flanks_length = args.min_flanks_length
    max_flanks_length = args.max_flanks_length

    num_processes = args.threads

    mfe_threshold_hpa = args.mfe_threshold_hpa
    mfe_threshold_hpt = args.mfe_threshold_hpt

    print(f"Writing output to: {path_to_output}")

    record_list = []
    feature_list = []
    intervaled_features_list = []

    # Try opening file as proper gbk
    print("Reading gbk...")
    try:
        with open(path_to_gbk) as handle:
            for record in GenBank.parse(handle):
                record_list.append(record)

        shutil.copy(path_to_gbk, path_to_output)

    # Try to fix bp LOCUS string issue
    except ValueError:
        print("Trying to fix possible prokka LOCUS input error...")
        modify_first_line(path_to_gbk, path_to_output)

        with open(path_to_output) as handle:
            for record in GenBank.parse(handle):
                record_list.append(record)

    for feature in record_list[0].features:
        feature_list.append(feature)

    for feature in tqdm(
        feature_list, desc="Processing annotations...", bar_format=BAR_FORMAT
    ):
        location = (
            feature.location.replace("complement(", "")
            .replace(">", "")
            .replace("<", "")
            .replace(")", "")
            .split("..")  # get length for normal and complement instances
        )
        start, stop = int(location[0]), int(location[1])
        feature.length = stop - start

        intervaled_features_list.append(IntervaledFeature([start, stop], feature))

    merged_intervals = merge_intervals(intervaled_features_list[1:])
    uncovered_intervals = find_uncovered_intervals(
        intervaled_features_list[0].interval, merged_intervals
    )

    # Filter by gaps length
    filtered_intervals = filter_intervals_by_length(
        uncovered_intervals, min_gap_length, max_gap_length
    )

    # Filter by flanks length
    filtered_intervals = filter_intervals_by_flanking_legth(
        filtered_intervals, min_flanks_length, max_flanks_length
    )

    output_intervals_strings = create_uncovered_intervals_str(filtered_intervals)

    sequence = Seq(record_list[0].sequence)

    print("Annotating RNA hairpins...")
    print("This might take a while...")

    print("Calculating batches...")
    subsequences = split_sequence(str(sequence.transcribe()), batch_size, OVERLAP_SIZE)
    num_batches = len(subsequences)
    print(f"Total # of batches: {num_batches}")

    print("Processing forward strand...")

    top_hairpins_f, all_hairpins_f = [], []

    for indx, sequence_ in enumerate(subsequences):
        top_hairpins_f_subs, all_hairpins_f_subs = find_top_hairpins(
            sequence_,
            mfe_threshold_hpt=mfe_threshold_hpt,
            mfe_threshold_hpa=mfe_threshold_hpa,
            batch_num=indx + 1,
            num_processes=num_processes,
        )

        top_hairpins_f_subs = [
            (
                x[0] + indx * batch_size - indx * OVERLAP_SIZE,
                x[1] + indx * batch_size - indx * OVERLAP_SIZE,
                x[2],
                x[3],
                x[4],
            )
            for x in top_hairpins_f_subs
        ]
        all_hairpins_f_subs = [
            (
                x[0] + indx * batch_size - indx * OVERLAP_SIZE,
                x[1] + indx * batch_size - indx * OVERLAP_SIZE,
                x[2],
                x[3],
                x[4],
            )
            for x in all_hairpins_f_subs
        ]

        top_hairpins_f.append(top_hairpins_f_subs)
        all_hairpins_f.append(all_hairpins_f_subs)

    top_hairpins_f = list(reduce(lambda x, y: x + y, top_hairpins_f, []))
    all_hairpins_f = list(reduce(lambda x, y: x + y, all_hairpins_f, []))

    all_hairpins_f = merge_similar_hairpins(all_hairpins_f, overlap_threshold=0.95)
    top_hairpins_f = merge_similar_hairpins(
        top_hairpins_f,
        overlap_threshold=0.95,
        message="Merging similar top hairpins...",
    )

    print("Processing reverse strand...")

    subsequences = split_sequence(str(sequence.complement_rna()), batch_size, 1000)

    top_hairpins_r, all_hairpins_r = [], []

    for indx, sequence_ in enumerate(subsequences):
        top_hairpins_r_subs, all_hairpins_r_subs = find_top_hairpins(
            sequence_,
            mfe_threshold_hpt=mfe_threshold_hpt,
            mfe_threshold_hpa=mfe_threshold_hpa,
            batch_num=indx + 1,
            num_processes=num_processes,
        )

        top_hairpins_r_subs = [
            (
                x[0] + indx * batch_size - indx * OVERLAP_SIZE,
                x[1] + indx * batch_size - indx * OVERLAP_SIZE,
                x[2],
                x[3],
                x[4],
            )
            for x in top_hairpins_r_subs
        ]
        all_hairpins_r_subs = [
            (
                x[0] + indx * batch_size - indx * OVERLAP_SIZE,
                x[1] + indx * batch_size - indx * OVERLAP_SIZE,
                x[2],
                x[3],
                x[4],
            )
            for x in all_hairpins_r_subs
        ]

        top_hairpins_r.append(top_hairpins_r_subs)
        all_hairpins_r.append(all_hairpins_r_subs)

    top_hairpins_r = list(reduce(lambda x, y: x + y, top_hairpins_r, []))
    all_hairpins_r = list(reduce(lambda x, y: x + y, all_hairpins_r, []))

    all_hairpins_r = merge_similar_hairpins(all_hairpins_r, overlap_threshold=0.95)
    top_hairpins_r = merge_similar_hairpins(
        top_hairpins_r,
        overlap_threshold=0.95,
        message="Merging similar top hairpins...",
    )

    print("Finding insertion sites...")
    insertion_sites_f = find_insertion_sites(filtered_intervals, top_hairpins_f)
    insertion_sites_r = find_insertion_sites(filtered_intervals, top_hairpins_r)

    insertion_sites_overlapping = find_overlapping_insertion_sites(
        insertion_sites_f + insertion_sites_r
    )

    print("Writing output...")
    insert_before_origin(
        path_to_output,
        (
            output_intervals_strings
            + create_hairpins_str(top_hairpins_f, "hpt")
            + create_hairpins_str(all_hairpins_f, "hpa")
            + create_hairpins_str(top_hairpins_r, "hpt", complement=True)
            + create_hairpins_str(all_hairpins_r, "hpa", complement=True)
            + create_insertion_sites_str(insertion_sites_overlapping)
        ),
    )

    print("Finished!")


if __name__ == "__main__":
    main()
