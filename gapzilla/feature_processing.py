from Bio.SeqFeature import SeqFeature, FeatureLocation

from gapzilla.models import IntervaledGap, InsertionSite


def create_uncovered_intervals_feature(
    intervals: list[IntervaledGap],
) -> list[SeqFeature]:
    """
    Create a list of SeqFeature objects representing uncovered intervals.

    Parameters
    ----------
    intervals : list of IntervaledGap
        A list of IntervaledGap objects, each representing an interval with a start and end position.

    Returns
    -------
    list of SeqFeature
        A list of SeqFeature objects, each corresponding to an interval from the input list.
        Each SeqFeature has a FeatureLocation with the start and end positions of the interval,
        a type of "gap", and a qualifier with a label "gap".
    """

    return [
        SeqFeature(
            FeatureLocation(interval.start, interval.end, strand=1),
            type="gap",
            qualifiers={
                "label": "gap",
            },
        )
        for interval in intervals
    ]


def create_hairpins_feature(hairpins, type_of_hairpin, complement=None):
    if complement:
        strand = -1
    else:
        strand = 1

    hairpins_feature_list = []

    for hairpin in hairpins:
        qual = {
            "label": type_of_hairpin,
            "sequence": hairpin[2],
            "structure": hairpin[3],
            "mfe": hairpin[4],
        }

        hairpins_feature_list.append(
            SeqFeature(
                FeatureLocation(hairpin[0], hairpin[1], strand=strand),
                qualifiers=qual,
                type="misc_feature",
            )
        )

    return hairpins_feature_list


def create_insertion_sites_feature(
    insertion_sites: list[InsertionSite],
) -> list[SeqFeature]:

    return [
        SeqFeature(
            FeatureLocation(insertion_site.start, insertion_site.end, strand=1),
            qualifiers={
                "score": insertion_site.score,
                "label": f"ins{insertion_site.score}",
            },
            type="misc_feature",
        )
        for insertion_site in insertion_sites
    ]
