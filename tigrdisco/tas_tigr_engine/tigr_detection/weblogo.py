"""Information content of repeat columns."""

from __future__ import annotations

from collections import Counter
from math import log, log2

import regex as re

from .extractor import fuzzy_pattern, score_candidate

__all__ = ["column_ic", "array_bits", "candidate_frames", "information_stacks"]

_ACGT = frozenset("ACGT")


def column_ic(sequences: list[str]) -> float:
    sequences = [s for s in sequences if s]
    n = len(sequences)
    if n < 2:
        return 0.0
    width = min(len(s) for s in sequences)
    correction = 3.0 / (2.0 * log(2) * n)
    total = 0.0
    for i in range(width):
        entropy = 0.0
        for count in Counter(s[i] for s in sequences).values():
            p = count / n
            entropy -= p * log2(p)
        total += max(0.0, 2.0 - entropy - correction)
    return total


def array_bits(units) -> float:
    return column_ic([str(u["R1"]) for u in units]) + column_ic([str(u["R2"]) for u in units])


def candidate_frames(dna: str, config, limit: int) -> list[tuple[str, str]]:
    scores: dict[str, float] = {}
    for k in config.candidate_lengths:
        for seq, count in Counter(dna[i : i + k] for i in range(len(dna) - k + 1)).items():
            if set(seq) - _ACGT or count < config.min_candidate_count:
                continue
            scores[seq] = score_candidate(seq, count)

    pairs: list[tuple[float, str, str]] = []
    ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
    for r1, r1_score in ranked[: config.top_r1_candidates]:
        query = fuzzy_pattern(r1, config.discovery_mismatches, config.discovery_allow_indels)
        span = len(r1) + config.discovery_spacer_max
        r2s: list[str] = []
        for match in re.finditer(query, dna):
            lo = match.start() + len(r1) + config.discovery_spacer_min
            hi = match.start() + span
            if hi + len(r1) > len(dna):
                continue
            for offset in range(lo, hi + 1):
                r2 = dna[offset : offset + len(r1)]
                if not (set(r2) - _ACGT):
                    r2s.append(r2)
        for r2, count in Counter(r2s).most_common(3):
            if r2 != r1:
                pairs.append((r1_score + score_candidate(r2, count), r1, r2))

    pairs.sort(key=lambda p: -p[0])
    out: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for _, r1, r2 in pairs:
        if (r1, r2) not in seen:
            seen.add((r1, r2))
            out.append((r1, r2))
            if len(out) >= limit:
                break
    return out


def information_stacks(sequences: list[str]) -> list[list[tuple[str, float]]]:
    sample_size = len(sequences)
    correction = 3 / (2 * log(2) * sample_size)
    stacks: list[list[tuple[str, float]]] = []
    for column in zip(*sequences):
        counts = Counter(column)
        frequencies = {base: counts[base] / sample_size for base in "ACGT"}
        entropy = -sum(value * log2(value) for value in frequencies.values() if value)
        information = max(0.0, 2.0 - entropy - correction)
        stacks.append(
            sorted(
                ((base, information * frequencies[base]) for base in "ACGT" if frequencies[base]),
                key=lambda item: item[1],
            )
        )
    return stacks
