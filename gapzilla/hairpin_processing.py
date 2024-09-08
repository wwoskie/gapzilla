import logging
import multiprocessing as mp
from functools import reduce

import RNA
from tqdm import tqdm

from gapzilla.config import config

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

    for hairpin in tqdm(
        hairpins,
        desc=message,
        disable=logging.root.level > logging.INFO,
        bar_format=BAR_FORMAT,
    ):
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
    num_processes=mp.cpu_count(),
):

    pool = mp.Pool(processes=num_processes)
    args = [
        (sequence, window_size, i)
        for i in range(0, len(sequence) - window_size + 1, step_size)
    ]
    all_hairpins = []

    for result in pool.imap_unordered(process_window, args):
        all_hairpins.extend(result)

    pool.close()
    pool.join()

    # all_hairpins = merge_similar_hairpins(all_hairpins)
    all_hairpins = [hp for hp in all_hairpins if hp[4] <= mfe_threshold_hpa]
    filtered_hairpins = [hp for hp in all_hairpins if hp[4] <= mfe_threshold_hpt]
    all_hairpins = [
        hairpin for hairpin in all_hairpins if hairpin not in filtered_hairpins
    ]

    return filtered_hairpins, all_hairpins


def find_hairpins_in_subseqs(
    subsequences,
    mfe_threshold_hpt,
    mfe_threshold_hpa,
    num_processes,
    overlap_threshold=0.9,
):

    top_hairpins, all_hairpins = [], []
    for subsequence in tqdm(
        subsequences,
        desc="Processing subsequences...",
        disable=logging.root.level > logging.INFO,
        bar_format=BAR_FORMAT,
    ):

        sequence_ = subsequences[subsequence]["seq"]

        top_hairpins_subs, all_hairpins_subs = find_top_hairpins(
            sequence_,
            mfe_threshold_hpt=mfe_threshold_hpt,
            mfe_threshold_hpa=mfe_threshold_hpa,
            num_processes=num_processes,
        )

        top_hairpins_subs = [
            (
                x[0] + subsequences[subsequence]["start"],
                x[1] + subsequences[subsequence]["start"],
                x[2],
                x[3],
                x[4],
            )
            for x in top_hairpins_subs
        ]

        all_hairpins_subs = [
            (
                x[0] + subsequences[subsequence]["start"],
                x[1] + subsequences[subsequence]["start"],
                x[2],
                x[3],
                x[4],
            )
            for x in all_hairpins_subs
        ]

        top_hairpins.append(top_hairpins_subs)
        all_hairpins.append(all_hairpins_subs)

    top_hairpins = list(reduce(lambda x, y: x + y, top_hairpins, []))
    all_hairpins = list(reduce(lambda x, y: x + y, all_hairpins, []))

    all_hairpins = merge_similar_hairpins(
        all_hairpins, overlap_threshold=overlap_threshold
    )
    top_hairpins = merge_similar_hairpins(
        top_hairpins,
        overlap_threshold=overlap_threshold,
        message="Merging similar top hairpins...",
    )

    return top_hairpins, all_hairpins
