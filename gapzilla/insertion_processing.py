"""
Functions to create and handle insertion sites
"""

from gapzilla.models import InsertionSite, IntervaledGap, Hairpin


def find_insertion_sites(
    gaps: list[IntervaledGap], coverages: list[Hairpin]
) -> list[InsertionSite]:
    """
    Identify potential insertion sites in gaps not covered by hairpin structures.

    Parameters
    ----------
    gaps : list of IntervaledGap
        A list of `IntervaledGap` objects representing gaps in the genome.
        Each `IntervaledGap` should have `start` and `end` attributes.
    coverages : list of Hairpin
        A list of `Hairpin` objects representing regions covered by hairpin structures.

    Returns
    -------
    list of InsertionSite
        A list of `InsertionSite` objects representing potential insertion sites.
        Each `InsertionSite` includes `start`, `end`, and `score` to indicate
        the insertion site's position and its score based on neighbouring `hpt` regions.
    """

    uncovered_intervals = []
    coverages = [(c.start, c.end) for c in sorted(coverages, key=lambda x: x.start)]

    for gap in gaps:
        g_start, g_end = gap.start, gap.end
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

        if left_neighbour and right_neighbour:
            score = 2
        elif left_neighbour or right_neighbour:
            score = 1
        else:
            continue

        insertion_sites.append(InsertionSite(u_start, u_end, score))

    return insertion_sites


def find_overlapping_insertion_sites(
    insertion_sites: list[InsertionSite],
) -> list[InsertionSite]:
    """
    Identify overlapping insertion sites and calculate their cumulative scores.

    Parameters
    ----------
    insertion_sites : list of InsertionSite
        A list of `InsertionSite` objects representing potential insertion sites.
        Each `InsertionSite` should have `start`, `end`, and `score` attributes.

    Returns
    -------
    list of InsertionSite
        A list of `InsertionSite` objects, including the original and newly created
        overlapping insertion sites with cumulative scores.
    """

    insertion_sites = sorted(insertion_sites, key=lambda x: x.start)
    result = []
    n = len(insertion_sites)

    for i in range(n):
        current_site = insertion_sites[i]
        current_start, current_end = current_site.start, current_site.end
        current_score = current_site.score

        # Always add the original interval first
        result.append(current_site)

        for j in range(i + 1, n):
            next_site = insertion_sites[j]
            next_start, next_end = next_site.start, next_site.end

            if next_start > current_end:
                break  # No more overlaps possible

            if current_start <= next_end and current_end >= next_start:
                overlap_start = max(current_start, next_start)
                overlap_end = min(current_end, next_end)
                overlap_score = current_score + next_site.score
                result.append(InsertionSite(overlap_start, overlap_end, overlap_score))

    return result


def filter_insertion_sites_by_max_score(
    insertion_sites: list[InsertionSite],
) -> list[InsertionSite]:
    """
    Filter a list of insertion sites to keep only one site with the maximum score
    if there are two or more sites with the exact same positions.

    Parameters
    ----------
    insertion_sites : list of InsertionSite
        A list of `InsertionSite` objects representing potential insertion sites.
        Each `InsertionSite` should have `start`, `end`, and `score` attributes.

    Returns
    -------
    list of InsertionSite
        A filtered list of `InsertionSite` objects where only one site with the maximum
        score is kept for each unique position.
    """

    # Dictionary to track the maximum score for each position
    position_to_max_site = {}

    for site in insertion_sites:
        key = (site.start, site.end)
        if key not in position_to_max_site:
            position_to_max_site[key] = site
        else:
            if site.score > position_to_max_site[key].score:
                position_to_max_site[key] = site

    # Extract the filtered list of insertion sites
    filtered_sites = list(position_to_max_site.values())

    return filtered_sites


def filter_insertion_sites_by_hairpins(
    insertion_sites: list[InsertionSite], hairpins: list[Hairpin]
) -> list[InsertionSite]:
    """
    Filter insertion sites based on coverage by hairpins.

    An insertion site will be dropped if it is fully covered by any hairpin.
    If an insertion site is partially covered by hairpins, it will be kept.

    Parameters
    ----------
    insertion_sites : list of InsertionSite
        List of `InsertionSite` objects representing potential insertion sites.
    hairpins : list of Hairpin
        List of `Hairpin` objects representing hairpin regions.

    Returns
    -------
    list of InsertionSite
        A filtered list of `InsertionSite` objects that are not fully covered
        by any `Hairpin` objects.
    """
    filtered_sites = []

    for site in insertion_sites:
        site_covered = False
        for hairpin in hairpins:
            if hairpin.start <= site.start and hairpin.end >= site.end:
                # Insertion site is fully covered by the hairpin
                site_covered = True
                break
        if not site_covered:
            filtered_sites.append(site)

    return filtered_sites
