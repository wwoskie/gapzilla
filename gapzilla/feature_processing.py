"""
Functions to transform data to biopython-acceptable `SeqRecord` format for .gbk files.
"""

from Bio.SeqFeature import SeqFeature, FeatureLocation

from gapzilla.models import IntervaledGap, InsertionSite, Hairpin


def create_hairpins_feature(
    hairpins: list[Hairpin], type_of_hairpin: str, complement: bool = None
) -> list[SeqFeature]:
    """
    Create a list of SeqFeature objects from a list of Hairpin objects.

    Parameters
    ----------
    hairpins : list of Hairpin
        A list of Hairpin objects, each containing information about a hairpin structure.
    type_of_hairpin : str
        A string representing the type of hairpin: `hpa` or `hpt`.
    complement : bool, optional
        If True, the strand will be set to -1 (complementary strand). If False or None, the strand will be set to 1 (forward strand).

    Returns
    -------
    list of SeqFeature
        A list of SeqFeature objects, each representing a hairpin with associated metadata.
    """

    if complement:
        strand = -1
    else:
        strand = 1

    hairpins_feature_list = []

    for hairpin in hairpins:
        qual = {
            "label": type_of_hairpin,
            # "sequence": hairpin.sequence,
            "structure": hairpin.structure,
            "mfe": hairpin.mfe,
        }

        hairpins_feature_list.append(
            SeqFeature(
                FeatureLocation(hairpin.start, hairpin.end, strand=strand),
                qualifiers=qual,
                type="misc_feature",
            )
        )

    return hairpins_feature_list


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


def create_insertion_sites_feature(
    insertion_sites: list[InsertionSite],
) -> list[SeqFeature]:
    """
    Create a list of sequence features from insertion sites.

    Parameters
    ----------
    insertion_sites : list of InsertionSite
        A list of `InsertionSite` objects representing the insertion sites.
        Each `InsertionSite` object should have `start`, `end`, and `score`
        attributes.

    Returns
    -------
    list of SeqFeature
        A list of `SeqFeature` objects created from the given insertion sites.
        Each `SeqFeature` object will have a `FeatureLocation` from `start` to
        `end` with a strand of 1, and qualifiers containing the `score` of the
        insertion site and a `label` in the format of "ins<score>".
    """

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
