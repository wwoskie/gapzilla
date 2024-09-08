import logging

from Bio.Seq import Seq
from tqdm import tqdm

from gapzilla.config import BAR_FORMAT
from gapzilla.models import IntervaledGap, IntervaledFeature


def split_sequence(
    sequence: Seq | str, gaps: IntervaledGap, border_shift: int = 75
) -> list[str]:

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


def find_uncovered_intervals(start, end, intervaled_features) -> list[IntervaledGap]:

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
            uncovered_interval = (current_start + 1, i_start - 1)
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
        uncovered_interval = (current_start + 1, end - 1)
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


def filter_intervals_by_length(intervals, min_length, max_length):
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


def filter_intervals_by_flanking_legth(intervals, min_flanks_length, max_flanks_length):
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
