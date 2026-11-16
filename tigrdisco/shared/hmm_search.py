"""Bundled Tas-core profile and interval helpers."""

from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "MITTas.hmm3"


def union_length(intervals):
    total = 0
    cs = ce = None
    for a, b in sorted(intervals):
        if cs is None:
            cs, ce = a, b
        elif a <= ce + 1:
            ce = max(ce, b)
        else:
            total += ce - cs + 1
            cs, ce = a, b
    return total + (0 if cs is None else ce - cs + 1)
