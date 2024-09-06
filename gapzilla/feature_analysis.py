from Bio.SeqFeature import SeqFeature, FeatureLocation


def create_uncovered_intervals_feature(intervals):

    return [
        SeqFeature(
            FeatureLocation(interval.interval[0], interval.interval[1], strand=1),
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


def create_insertion_sites_feature(insertion_sites):

    return [
        SeqFeature(
            FeatureLocation(
                insertion_site.interval[0], insertion_site.interval[1], strand=1
            ),
            qualifiers={
                "score": insertion_site.score,
                "label": f"ins{insertion_site.score}",
            },
            type="misc_feature",
        )
        for insertion_site in insertion_sites
    ]
