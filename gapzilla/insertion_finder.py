from gapzilla.models import InsertionSite


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
