"""
Functions to process intervals and gaps
"""

import logging

from Bio.Seq import Seq
from tqdm import tqdm

from gapzilla.config import BAR_FORMAT
from gapzilla.models import IntervaledGap, IntervaledFeature


def split_sequence(
    sequence: Seq | str, gaps: IntervaledGap, border_shift: int = 75
) -> list[Seq | str]:
    """
    Split the given sequence based on intervals and a specified border shift.

    Parameters
    ----------
    sequence : Seq or str
        The sequence to be split.
    gaps : IntervaledGap
        An iterable containing gap intervals with `start` and `end` attributes.
    border_shift : int, optional
        The number of positions to extend the borders of each gap interval. Default is 75.

    Returns
    -------
    list of Seq or str
        A list of subsequences extracted from the original sequence, each extended by the border shift value.
    """

    subsequences = {}
    for indx, gap in enumerate(gaps):
        subsequences[indx] = {}
        subsequences[indx]["seq"] = sequence[
            slice(gap.start - border_shift, gap.end + border_shift)
        ]
        subsequences[indx]["start"] = gap.start - border_shift

    return subsequences


def merge_intervals(
    intervaled_features: list[IntervaledFeature],
) -> list[IntervaledFeature]:
    """
    Merge overlapping intervals from a list of intervaled features.

    Parameters
    ----------
    intervaled_features : list of IntervaledFeature
        A list of intervaled features to be merged. Each feature must have `start` and `end` attributes,
        as well as `feature_list`, `feature_num`, and `feature_lengths` attributes.

    Returns
    -------
    list of IntervaledFeature
        A list of merged intervaled features with updated attributes including combined feature lists,
        counts, and lengths.
    """

    sorted_intervaled_features = sorted(intervaled_features, key=lambda x: x.start)

    merged_intervaled_features = []

    for intervaled_feature in tqdm(
        sorted_intervaled_features,
        desc="Merging intervals...",
        disable=logging.root.level > logging.INFO,
        bar_format=BAR_FORMAT,
    ):
        if (
            not merged_intervaled_features
            or merged_intervaled_features[-1].end < intervaled_feature.start
        ):
            merged_intervaled_features.append(intervaled_feature)
        else:
            merged_intervaled_features[-1].end = max(
                merged_intervaled_features[-1].end,
                intervaled_feature.end,
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


def find_uncovered_intervals(
    start: int, end: int, intervaled_features: list[IntervaledFeature]
) -> list[IntervaledGap]:
    """
    Identify intervals within a specified range that are not covered by given features.

    Parameters
    ----------
    start : int
        The start position of the interested range.
    end : int
        The end position of the interested range.
    intervaled_features : list of IntervaledFeature
        A list of intervaled features, each having `start` and `end` attributes.

    Returns
    -------
    list of IntervaledGap
        A list of uncovered intervals, each represented as an `IntervaledGap` object.
        Each gap object includes information about the interval and any adjacent features on the left and right.
    """

    intervaled_features.sort(key=lambda x: x.start)  # Sort by intervals

    uncovered_intervals = []
    current_start = start
    prev_feature = None

    for feature in tqdm(
        intervaled_features,
        disable=logging.root.level > logging.INFO,
        desc="Finding uncovered intervals...",
        bar_format=BAR_FORMAT,
    ):
        i_start, i_end = feature.start, feature.end
        if i_end <= start:
            prev_feature = feature
            continue
        if i_start >= end:
            break

        if i_start > current_start and (i_start - 1) != current_start:
            uncovered_interval = (current_start, i_start)
            features_left = (
                prev_feature
                if prev_feature and prev_feature.end <= uncovered_interval[0]
                else None
            )
            features_right = feature if feature.start >= uncovered_interval[1] else None
            uncovered_intervals.append(
                IntervaledGap(
                    uncovered_interval[0],
                    uncovered_interval[1],
                    features_left,
                    features_right,
                )
            )

        current_start = max(current_start, i_end)
        prev_feature = feature

    if current_start < end:
        uncovered_interval = (current_start, end)
        features_left = (
            prev_feature
            if prev_feature and prev_feature.end <= uncovered_interval[0]
            else None
        )
        features_right = None
        uncovered_intervals.append(
            IntervaledGap(
                uncovered_interval[0],
                uncovered_interval[1],
                features_left,
                features_right,
            )
        )

    return uncovered_intervals


def filter_intervals_by_length(
    intervals: list[IntervaledGap], min_length: int, max_length: int
) -> list[IntervaledGap]:
    """
    Filter a list of intervaled gaps based on specified minimum and maximum lengths.

    Parameters
    ----------
    intervals : list of IntervaledGap
        A list of intervaled gaps, each with `start` and `end` attributes.
    min_length : int
        The minimum length an interval must have to be included in the output.
    max_length : int
        The maximum length an interval can have to be included in the output.

    Returns
    -------
    list of IntervaledGap
        A list of intervaled gaps that satisfy the specified length constraints.
    """

    filtered_intervals = []
    for intervaled_gap in tqdm(
        intervals,
        desc="Filtering intervals by length...",
        disable=logging.root.level > logging.INFO,
        bar_format=BAR_FORMAT,
    ):

        if intervaled_gap.start and intervaled_gap.end:
            start, end = intervaled_gap.start, intervaled_gap.end
            length = end - start

            if min_length <= length + 1 <= max_length:
                filtered_intervals.append(intervaled_gap)

    return filtered_intervals


def filter_intervals_by_flanking_legth(
    intervals: list[IntervaledGap], min_flanks_length: int, max_flanks_length: int
):
    """
    Filter a list of intervaled gaps based on the lengths of their flanking features.

    Parameters
    ----------
    intervals : list of IntervaledGap
        A list of intervaled gaps, each potentially having `features_left` and `features_right` attributes.
    min_flanks_length : int
        The minimum length a flanking feature must have for the gap to be included.
    max_flanks_length : int
        The maximum length a flanking feature can have for the gap to be included.

    Returns
    -------
    list of IntervaledGap
        A list of intervaled gaps where both flanking features fall within the specified length constraints.
    """

    filtered_intervals = []
    for intervaled_gap in tqdm(
        intervals,
        disable=logging.root.level > logging.INFO,
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
