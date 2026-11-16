"""Arbitration between the dual-repeat and stem-loop results."""

import statistics

MODES = ("off", "unit", "covered", "intersect", "distance")


def dr_intervals(dr_arrays):
    units, spans = [], []
    for array in dr_arrays:
        block = []
        for unit in array.get("units", ()):
            start = unit["contig_index"]
            width = sum(len(unit.get(part, "")) for part in ("R1", "spacer_a", "R2", "spacer_b"))
            block.append((start, start + width))
        if block:
            units.extend(block)
            spans.append((min(s for s, _ in block), max(e for _, e in block)))
    return units, spans


def arbitrate_locus(
    sl_arrays, dr_arrays, window_offset, mode="off", gap=1200, stem=8, coverage=1.0
):
    if mode not in MODES:
        raise ValueError(f"Unknown arbitration mode: {mode}")
    if mode == "off" or not sl_arrays or not dr_arrays:
        return sl_arrays
    dr_units, dr_spans = dr_intervals(dr_arrays)
    if not dr_spans:
        return sl_arrays
    reach = 0 if mode == "intersect" else gap
    kept = []
    for array in sl_arrays:
        units = array.get("units", ())
        if not units:
            kept.append(array)
            continue
        if mode in ("unit", "covered"):
            mine = [(window_offset + u["start"] - 1, window_offset + u["end"]) for u in units]
            on = sum(1 for a, b in mine if any(min(b, d) - max(a, c) > 0 for c, d in dr_units))
            if (on > 0) if mode == "unit" else (on >= coverage * len(mine)):
                continue
            kept.append(array)
            continue
        if stem > 0 and statistics.median([u["stem_length"] for u in units]) >= stem:
            kept.append(array)
            continue
        start = window_offset + min(u["extended_start"] for u in units) - 1
        end = window_offset + max(u["extended_end"] for u in units)
        if any(
            (start < de and ds < end) if mode == "intersect" else max(ds - end, start - de) <= reach
            for ds, de in dr_spans
        ):
            continue
        kept.append(array)
    return kept
