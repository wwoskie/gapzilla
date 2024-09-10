"""
Funtions to find, filter and merge RNA-hairpins
"""

import logging
import multiprocessing as mp
from functools import reduce

import RNA
from Bio.Seq import Seq
from tqdm import tqdm

from gapzilla.config import BAR_FORMAT
from gapzilla.models import Hairpin


def find_hairpins(sequence: str | Seq, structure: str) -> list[Hairpin]:
    """
    Find hairpin structures in a given RNA secondary structure.

    Parameters
    ----------
    sequence : str or Seq
        The RNA sequence corresponding to the secondary structure.
    structure : str
        The RNA secondary structure in dot-bracket notation.

    Returns
    -------
    list of Hairpin
        A list of `Hairpin` objects representing the hairpin structures found in
        the given sequence and structure. Each `Hairpin` contains the start and end
        indices, the hairpin sequence, and the hairpin structure.

    Notes
    -----
    - The function scans through the `structure` string using a stack to identify
      matching parentheses, which indicate hairpin loops.
    - A hairpin is considered valid if the distance between the matched parentheses
      is greater than 3 nucleotides.
    """

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
                hairpins.append(Hairpin(start, end + 1, hairpin_seq, hairpin_struct))
    return hairpins


def evaluate_hairpins(hairpins: list[Hairpin]) -> list[Hairpin]:
    """
    Evaluate the minimum free energy (MFE) of a list of hairpin structures.

    Parameters
    ----------
    hairpins : list of Hairpin
        A list of `Hairpin` objects to be evaluated. Each `Hairpin` should have
        a `sequence` and `structure` attribute representing the hairpin's
        sequence and secondary structure, respectively.

    Returns
    -------
    list of Hairpin
        A list of `Hairpin` objects with the `mfe` attribute set to the calculated
        minimum free energy of the hairpin's secondary structure.
    """

    evaluated_hairpins = []
    for hairpin in hairpins:
        fc = RNA.fold_compound(hairpin.sequence)
        mfe = fc.eval_structure(hairpin.structure)
        hairpin.mfe = mfe
        evaluated_hairpins.append(hairpin)
    return evaluated_hairpins


def merge_similar_hairpins(
    hairpins: list[Hairpin],
    overlap_threshold: float = 0.9,
    backward_search_thres: int = 20,
    message: str = "Merging similar hairpins...",
) -> list[Hairpin]:
    """
    Merge similar hairpin structures based on their overlap.

    Parameters
    ----------
    hairpins : list of Hairpin
        A list of `Hairpin` objects to be merged. Each `Hairpin` should have
        attributes like `start`, `end`, `sequence`, `structure`, and `mfe`.
    overlap_threshold : float, optional
        The threshold for considering two hairpins as overlapping. Default is 0.9.
    backward_search_thres : int, optional
        The number of recent hairpins to search backward for overlapping. Default
        is 20.
    message : str, optional
        The message to display on the progress bar. Default is "Merging similar hairpins...".

    Returns
    -------
    list of Hairpin
        A list of `Hairpin` objects with similar and overlapping hairpins merged.
        In case of overlap, the hairpin with the lower minimum free energy (MFE)
        is retained.

    Notes
    -----
    - The function sorts the list of hairpins by their start and end positions.
    - It uses a sliding window approach to check for overlap within a specified
      number of recent hairpins.
    - When overlapping hairpins are found, it retains the one with the lower MFE
      and discards the other.
    """

    merged_hairpins = []
    hairpins.sort(key=lambda x: (x.start, x.end))  # Sort by start and end

    for hairpin in tqdm(
        hairpins,
        desc=message,
        disable=logging.root.level > logging.INFO,
        bar_format=BAR_FORMAT,
    ):
        start, end, seq, struct, mfe = (
            hairpin.start,
            hairpin.end,
            hairpin.sequence,
            hairpin.structure,
            hairpin.mfe,
        )
        overlap_found = False
        for merged in merged_hairpins[-backward_search_thres:]:
            m_start, m_end, m_seq, m_struct, m_mfe = (
                merged.start,
                merged.end,
                merged.sequence,
                merged.structure,
                merged.mfe,
            )
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


def process_window(args: tuple[str, int, int]) -> list[Hairpin]:
    """
    Process a sliding window on a sequence to find and evaluate hairpin structures.

    Parameters
    ----------
    args : tuple of (str, int, int)
        A tuple containing the sequence, window size, and start position.
        - sequence: The RNA sequence to be processed.
        - window_size: The size of the window to slide over the sequence.
        - start: The starting position of the window on the sequence.

    Returns
    -------
    list of Hairpin
        A list of `Hairpin` objects found and evaluated within the given window.
        Each hairpin's start and end positions are adjusted to the original
        sequence coordinates.

    Notes
    -----
    - The function extracts a subsequence of the given window size starting
      from the provided start position.
    - RNA folding is performed on the window sequence to predict its secondary
      structure and minimum free energy (MFE).
    - Hairpin structures are identified within the window sequence and evaluated
      for their MFE.
    - The start and end positions of each identified hairpin are adjusted to match
      the original sequence coordinates.
    """

    sequence, window_size, start = args
    window_seq = sequence[start : start + window_size]
    structure, mfe = RNA.fold(window_seq)
    hairpins = find_hairpins(window_seq, structure)
    evaluated_hairpins = evaluate_hairpins(hairpins)

    adjusted_hairpins = []
    for hairpin in evaluated_hairpins:
        hairpin.start = hairpin.start + start
        hairpin.end = hairpin.end + start
        adjusted_hairpins.append(hairpin)

    return adjusted_hairpins


def find_top_hairpins(
    sequence: str | Seq,
    window_size: int = 50,
    step_size: int = 25,
    mfe_threshold_hpt: float = -15.0,
    mfe_threshold_hpa: float = -1,
    num_processes: int = mp.cpu_count(),
) -> tuple[list[Hairpin], list[Hairpin]]:
    """
    Find and evaluate top hairpin structures in an RNA sequence using a sliding window approach.

    Parameters
    ----------
    sequence : str or Seq
        The RNA sequence to be analyzed for hairpin structures.
    window_size : int, optional
        The size of the window to slide over the sequence. Default is 50.
    step_size : int, optional
        The step size for the sliding window. Default is 25.
    mfe_threshold_hpt : float, optional
        The threshold for minimum free energy (MFE) to classify hairpins as
        'top' hairpins. Default is -15.0.
    mfe_threshold_hpa : float, optional
        The threshold for MFE to filter hairpins after merging similar ones.
        Default is -1.
    num_processes : int, optional
        The number of processes to use for parallel processing. Default is the
        number of CPU cores available.

    Returns
    -------
    tuple of (list of Hairpin, list of Hairpin)

        - filtered_hairpins: A list of 'top' Hairpin objects that meet the MFE
          threshold specified by `mfe_threshold_hpt`.
        - all_hairpins: A list of remaining Hairpin objects after filtering out
          the 'top' hairpins.

    Notes
    -----
    - The function uses a sliding window approach to scan the sequence for hairpin
      structures. The windows are processed in parallel using multiple processes.
    - RNA folding is performed within each window to predict secondary structures
      and evaluate their MFE.
    - Similar hairpins within the sequence are merged.
    - Hairpins are filtered based on their MFE, and those meeting the 'top' MFE
      threshold are separated from the remaining hairpins.
    """

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
    all_hairpins = [hp for hp in all_hairpins if hp.mfe <= mfe_threshold_hpa]
    filtered_hairpins = [hp for hp in all_hairpins if hp.mfe <= mfe_threshold_hpt]
    all_hairpins = [
        hairpin for hairpin in all_hairpins if hairpin not in filtered_hairpins
    ]

    return filtered_hairpins, all_hairpins


def find_hairpins_in_subseqs(
    subsequences: list[str | Seq],
    mfe_threshold_hpt: float,
    mfe_threshold_hpa: float,
    num_processes: int = mp.cpu_count(),
    overlap_threshold: float = 0.9,
) -> tuple[list[Hairpin], list[Hairpin]]:
    """
    Find and evaluate hairpin structures in multiple RNA subsequences.

    Parameters
    ----------
    subsequences : list of str or Seq
        A list of RNA subsequences, each represented as a string or Seq object.
    mfe_threshold_hpt : float
        The threshold for minimum free energy (MFE) to classify hairpins as 'top'
        hairpins.
    mfe_threshold_hpa : float
        The threshold for MFE to filter hairpins after merging similar ones.
    num_processes : int, optional
        The number of processes to use for parallel processing. Default is the
        number of CPU cores available.
    overlap_threshold : float, optional
        The threshold for considering two hairpins as overlapping. Default is 0.9.

    Returns
    -------
    tuple of (list of Hairpin, list of Hairpin)
        - top_hairpins: A list of 'top' Hairpin objects that meet the MFE threshold
          specified by `mfe_threshold_hpt`.
        - all_hairpins: A list of remaining Hairpin objects after filtering out the
          'top' hairpins and merging similar hairpins.

    Notes
    -----
    - The function processes each subsequence to find and evaluate hairpin structures
      using parallel processing.
    - RNA folding is performed within each subsequence, and hairpins are evaluated
      based on their MFE.
    - The start and end positions of each identified hairpin are adjusted to match the
      original sequence coordinates of the subsequences.
    - Similar hairpins within each subsequence are merged.
    - Finally, hairpins are filtered based on their MFE, and those meeting the 'top'
      MFE threshold are separated from the remaining hairpins.
    """

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

        top_hairpins_subs_corrected = []
        for hairpin in top_hairpins_subs:
            hairpin.start += subsequences[subsequence]["start"]
            hairpin.end += subsequences[subsequence]["start"]
            top_hairpins_subs_corrected.append(hairpin)

        all_hairpins_subs_corrected = []
        for hairpin in all_hairpins_subs:
            hairpin.start += subsequences[subsequence]["start"]
            hairpin.end += subsequences[subsequence]["start"]
            all_hairpins_subs_corrected.append(hairpin)

        top_hairpins.append(top_hairpins_subs_corrected)
        all_hairpins.append(all_hairpins_subs_corrected)

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


def merge_overlapping_hairpins(hairpins: list[Hairpin]) -> list[Hairpin]:
    """
    Merge overlapping hairpins into a single hairpin.

    This function takes a list of Hairpin objects and merges any that overlap
    (i.e., if the `end` of one hairpin is greater than or equal to the `start`
    of the next hairpin, they are considered overlapping). The merged hairpin
    will have the lowest start and the highest end among the overlapping ones.

    Parameters
    ----------
    hairpins : list[Hairpin]
        A list of Hairpin objects with `.start` and `.end` attributes.

    Returns
    -------
    list[Hairpin]
        A list of Hairpin objects with overlapping hairpins merged.
    """

    # Sort hairpins based on the start position
    sorted_hairpins = sorted(hairpins, key=lambda x: x.start)
    merged_hairpins = [sorted_hairpins[0]]

    for current in sorted_hairpins[1:]:
        last_merged = merged_hairpins[-1]

        if current.start <= last_merged.end:
            # Overlapping, merge intervals
            last_merged.end = max(last_merged.end, current.end)
        else:
            # No overlap, add to result
            merged_hairpins.append(current)

    return merged_hairpins
