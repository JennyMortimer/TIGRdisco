"""Unit-by-unit reconstruction of a dual-repeat array."""

import regex as re
from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class TigrConfig:
    candidate_lengths: tuple = (10,)
    min_candidate_count: int = 2
    top_r1_candidates: int = 15
    discovery_mismatches: int = 2
    discovery_allow_indels: bool = False
    discovery_spacer_min: int = 7
    discovery_spacer_max: int = 13

    window_max_errors: int = 4
    window_allow_indels: bool = True

    unit_max_errors: int = 2
    unit_allow_indels: bool = True
    unit_spacer_min: int = 6
    unit_search_end_offset: int = 14
    terminal_spacer_b_length: int = 12

    reorganize_misaligned_unit: bool = True
    reorganize_unit_length: int = 36
    reorganize_r2_max_mismatches: int = 4

    prefix_max_errors: int = 2
    prefix_allow_indels: bool = False
    prefix_search_distance: int = 40
    prefix_nearest_distance: int = 20
    prefix_r1_min_length: int = 8
    prefix_r1_max_length: int = 13

    suffix_max_errors: int = 3
    suffix_allow_indels: bool = False
    suffix_search_distance: int = 40
    suffix_min_remaining_length: int = 10


def fuzzy_pattern(sequence, max_errors, allow_indels):
    error_type = "e" if allow_indels else "s"
    return f"({re.escape(sequence)}){{{error_type}<={max_errors}}}"


def count_substitution_mismatches(sequence, reference):
    return sum(a != b for a, b in zip(sequence, reference)) + abs(len(sequence) - len(reference))


def reorganize_unit_if_better(
    dna_seq,
    unit_start,
    next_r1_pos,
    actual_r1,
    core_r2,
    spacer_a,
    actual_r2,
    spacer_b,
    reference_spacer_a_length,
    reference_spacer_b_length,
    config,
):
    if not config.reorganize_misaligned_unit:
        return spacer_a, actual_r2, spacer_b
    if reference_spacer_a_length is None or reference_spacer_b_length is None:
        return spacer_a, actual_r2, spacer_b
    if next_r1_pos - unit_start != config.reorganize_unit_length:
        return spacer_a, actual_r2, spacer_b
    if (len(spacer_a), len(spacer_b)) == (
        reference_spacer_a_length,
        reference_spacer_b_length,
    ):
        return spacer_a, actual_r2, spacer_b

    r1_end = unit_start + len(actual_r1)
    new_r2_start = r1_end + reference_spacer_a_length
    new_r2_end = new_r2_start + len(core_r2)
    if new_r2_end > next_r1_pos:
        return spacer_a, actual_r2, spacer_b

    new_spacer_a = dna_seq[r1_end:new_r2_start]
    new_r2 = dna_seq[new_r2_start:new_r2_end]
    new_spacer_b = dna_seq[new_r2_end:next_r1_pos]
    if len(new_spacer_b) != reference_spacer_b_length:
        return spacer_a, actual_r2, spacer_b

    old_mismatches = count_substitution_mismatches(actual_r2, core_r2)
    new_mismatches = count_substitution_mismatches(new_r2, core_r2)
    if new_mismatches <= config.reorganize_r2_max_mismatches and new_mismatches <= old_mismatches:
        return new_spacer_a, new_r2, new_spacer_b

    return spacer_a, actual_r2, spacer_b


def find_best_match_in_window(
    dna, pattern, start_window, end_window, max_errors=4, allow_indels=True
):
    start_window = max(0, start_window)
    end_window = min(len(dna), end_window)
    search_area = dna[start_window:end_window]

    query = fuzzy_pattern(pattern, max_errors, allow_indels)
    match = re.search(query, search_area, flags=re.BESTMATCH)

    if match:
        return start_window + match.start(), sum(match.fuzzy_counts)
    return -1, 999


def score_candidate(seq, count):
    return count * len(seq)


def find_similar_r1_in_prefix(dna_seq, core_r1, first_r1_pos, config):
    prefix_r1_list = []
    current_front = first_r1_pos

    query = fuzzy_pattern(core_r1, config.prefix_max_errors, config.prefix_allow_indels)

    while current_front > 0:
        search_start = max(0, current_front - config.prefix_search_distance)
        search_end = current_front - config.prefix_nearest_distance

        if search_end - search_start < 5:
            break

        prefix_area = dna_seq[search_start:search_end]

        matches = list(re.finditer(query, prefix_area, flags=re.BESTMATCH))

        if not matches:
            break

        best_match = matches[-1]

        r1_pos_global = search_start + best_match.start()
        actual_r1 = best_match.group()

        if config.prefix_r1_min_length <= len(actual_r1) <= config.prefix_r1_max_length:
            prefix_r1_list.append((r1_pos_global, actual_r1))
            current_front = r1_pos_global
        else:
            current_front = search_start

    prefix_r1_list.sort(key=lambda x: x[0])
    return prefix_r1_list


def extract_tigr_array(dna_seq, r1, r2, config):
    units = []
    curr_pos, _ = find_best_match_in_window(
        dna_seq,
        r1,
        0,
        len(dna_seq),
        config.window_max_errors,
        config.window_allow_indels,
    )
    first_core_r1_pos = curr_pos
    last_core_r1_pos = curr_pos
    count = 0
    prefix_count = 0
    suffix_count = 0
    core_count = 0
    previous_spacer_b_length = None
    spacer_a_length_history = []
    spacer_b_length_history = []

    if first_core_r1_pos > 0:
        prefix_r1_list = find_similar_r1_in_prefix(dna_seq, r1, first_core_r1_pos, config)
        for prefix_r1_pos, actual_r1 in prefix_r1_list:
            len_r1 = len(actual_r1)
            r2_start = prefix_r1_pos + len_r1 + config.unit_spacer_min
            r2_end = prefix_r1_pos + len_r1 + config.unit_search_end_offset + len(r2)
            r2_idx, _ = find_best_match_in_window(
                dna_seq,
                r2,
                r2_start,
                r2_end,
                config.window_max_errors,
                config.window_allow_indels,
            )

            if r2_idx != -1:
                unit_r2_query = fuzzy_pattern(r2, config.unit_max_errors, config.unit_allow_indels)
                m2 = re.search(
                    unit_r2_query,
                    dna_seq[r2_idx : r2_idx + len(r2) + config.unit_max_errors],
                )
                if m2:
                    r2_idx += m2.start()
                    actual_r2 = m2.group()
                else:
                    actual_r2 = dna_seq[r2_idx : r2_idx + len(r2)]
                seq_sa = dna_seq[prefix_r1_pos + len_r1 : r2_idx]

                next_r1_idx, _ = find_best_match_in_window(
                    dna_seq,
                    r1,
                    r2_idx + len(actual_r2) + config.unit_spacer_min,
                    r2_idx + len(actual_r2) + config.unit_search_end_offset + len(r1),
                    config.window_max_errors,
                    config.window_allow_indels,
                )
                seq_sb = dna_seq[
                    r2_idx
                    + len(actual_r2) : (
                        next_r1_idx
                        if next_r1_idx != -1
                        else r2_idx + len(actual_r2) + config.terminal_spacer_b_length
                    )
                ]
                if next_r1_idx != -1:
                    previous_spacer_b_length = len(seq_sb)
                    spacer_a_length_history.append(len(seq_sa))
                    spacer_b_length_history.append(len(seq_sb))

                units.append(
                    {
                        "start": prefix_r1_pos,
                        "R1": actual_r1,
                        "spacer_a": seq_sa,
                        "R2": actual_r2,
                        "spacer_b": seq_sb,
                        "classification": "Prefix",
                    }
                )
                count += 1
                prefix_count += 1

    core_r1_pos_list = []
    while curr_pos != -1 and curr_pos < len(dna_seq):
        unit_r1_query = fuzzy_pattern(r1, config.unit_max_errors, config.unit_allow_indels)
        m1 = re.search(
            unit_r1_query,
            dna_seq[curr_pos : curr_pos + len(r1) + config.unit_max_errors],
        )
        if not m1:
            break
        unit_start = curr_pos + m1.start()
        actual_r1 = m1.group()
        len_r1 = len(actual_r1)
        core_r1_pos_list.append(unit_start)

        r2_start = unit_start + len_r1 + config.unit_spacer_min
        r2_end = unit_start + len_r1 + config.unit_search_end_offset + len(r2)
        r2_idx, _ = find_best_match_in_window(
            dna_seq,
            r2,
            r2_start,
            r2_end,
            config.window_max_errors,
            config.window_allow_indels,
        )

        if r2_idx == -1:
            break

        unit_r2_query = fuzzy_pattern(r2, config.unit_max_errors, config.unit_allow_indels)
        m2 = re.search(
            unit_r2_query,
            dna_seq[r2_idx : r2_idx + len(r2) + config.unit_max_errors],
        )
        if m2:
            r2_idx += m2.start()
            actual_r2 = m2.group()
        else:
            actual_r2 = dna_seq[r2_idx : r2_idx + len(r2)]
        seq_sa = dna_seq[unit_start + len_r1 : r2_idx]

        next_r1_idx, _ = find_best_match_in_window(
            dna_seq,
            r1,
            r2_idx + len(actual_r2) + config.unit_spacer_min,
            r2_idx + len(actual_r2) + config.unit_search_end_offset + len(r1),
            config.window_max_errors,
            config.window_allow_indels,
        )

        if next_r1_idx == -1:
            spacer_b_length = (
                previous_spacer_b_length
                if previous_spacer_b_length is not None
                else config.terminal_spacer_b_length
            )
            seq_sb = dna_seq[r2_idx + len(actual_r2) : r2_idx + len(actual_r2) + spacer_b_length]
            units.append(
                {
                    "start": unit_start,
                    "R1": actual_r1,
                    "spacer_a": seq_sa,
                    "R2": actual_r2,
                    "spacer_b": seq_sb,
                    "classification": "Core",
                }
            )
            count += 1
            core_count += 1
            last_core_r1_pos = unit_start + len(actual_r1)
            break
        else:
            seq_sb = dna_seq[r2_idx + len(actual_r2) : next_r1_idx]
            reference_spacer_a_length = (
                Counter(spacer_a_length_history).most_common(1)[0][0]
                if spacer_a_length_history
                else None
            )
            reference_spacer_b_length = (
                Counter(spacer_b_length_history).most_common(1)[0][0]
                if spacer_b_length_history
                else None
            )
            seq_sa, actual_r2, seq_sb = reorganize_unit_if_better(
                dna_seq,
                unit_start,
                next_r1_idx,
                actual_r1,
                r2,
                seq_sa,
                actual_r2,
                seq_sb,
                reference_spacer_a_length,
                reference_spacer_b_length,
                config,
            )
            previous_spacer_b_length = len(seq_sb)
            spacer_a_length_history.append(len(seq_sa))
            spacer_b_length_history.append(len(seq_sb))
            units.append(
                {
                    "start": unit_start,
                    "R1": actual_r1,
                    "spacer_a": seq_sa,
                    "R2": actual_r2,
                    "spacer_b": seq_sb,
                    "classification": "Core",
                }
            )
            curr_pos = next_r1_idx
            count += 1
            core_count += 1
            last_core_r1_pos = unit_start

    if core_r1_pos_list and last_core_r1_pos < len(dna_seq) - config.suffix_min_remaining_length:
        suffix_r1_list = find_similar_r1_in_suffix(dna_seq, r1, last_core_r1_pos, config)
        for suffix_r1_pos, actual_r1 in suffix_r1_list:
            len_r1 = len(actual_r1)
            r2_start = suffix_r1_pos + len_r1 + config.unit_spacer_min
            r2_end = suffix_r1_pos + len_r1 + config.unit_search_end_offset + len(r2)
            r2_idx, _ = find_best_match_in_window(
                dna_seq,
                r2,
                r2_start,
                r2_end,
                config.window_max_errors,
                config.window_allow_indels,
            )
            if r2_idx == -1:
                continue

            unit_r2_query = fuzzy_pattern(r2, config.unit_max_errors, config.unit_allow_indels)
            m2 = re.search(
                unit_r2_query,
                dna_seq[r2_idx : r2_idx + len(r2) + config.unit_max_errors],
            )
            if m2:
                r2_idx += m2.start()
                actual_r2 = m2.group()
            else:
                actual_r2 = dna_seq[r2_idx : r2_idx + len(r2)]
            seq_sa = dna_seq[suffix_r1_pos + len_r1 : r2_idx]

            next_r1_idx, _ = find_best_match_in_window(
                dna_seq,
                r1,
                r2_idx + len(actual_r2) + config.unit_spacer_min,
                r2_idx + len(actual_r2) + config.unit_search_end_offset + len(r1),
                config.window_max_errors,
                config.window_allow_indels,
            )
            if next_r1_idx == -1:
                spacer_b_length = (
                    previous_spacer_b_length
                    if previous_spacer_b_length is not None
                    else config.terminal_spacer_b_length
                )
                seq_sb = dna_seq[
                    r2_idx + len(actual_r2) : r2_idx + len(actual_r2) + spacer_b_length
                ]
            else:
                seq_sb = dna_seq[r2_idx + len(actual_r2) : next_r1_idx]
                reference_spacer_a_length = (
                    Counter(spacer_a_length_history).most_common(1)[0][0]
                    if spacer_a_length_history
                    else None
                )
                reference_spacer_b_length = (
                    Counter(spacer_b_length_history).most_common(1)[0][0]
                    if spacer_b_length_history
                    else None
                )
                seq_sa, actual_r2, seq_sb = reorganize_unit_if_better(
                    dna_seq,
                    suffix_r1_pos,
                    next_r1_idx,
                    actual_r1,
                    r2,
                    seq_sa,
                    actual_r2,
                    seq_sb,
                    reference_spacer_a_length,
                    reference_spacer_b_length,
                    config,
                )
                previous_spacer_b_length = len(seq_sb)
                spacer_a_length_history.append(len(seq_sa))
                spacer_b_length_history.append(len(seq_sb))

            seq_sb = seq_sb if seq_sb else "(Empty)"
            units.append(
                {
                    "start": suffix_r1_pos,
                    "R1": actual_r1,
                    "spacer_a": seq_sa,
                    "R2": actual_r2,
                    "spacer_b": seq_sb,
                    "classification": "Suffix",
                }
            )
            count += 1
            suffix_count += 1

    return {
        "total_count": count,
        "core_count": core_count,
        "prefix_count": prefix_count,
        "suffix_count": suffix_count,
        "units": tuple(sorted(units, key=lambda unit: unit["start"])),
    }


def find_similar_r1_in_suffix(dna_seq, core_r1, last_r1_pos, config):
    suffix_r1_list = []
    current_back = last_r1_pos

    max_errors = config.suffix_max_errors
    query = fuzzy_pattern(core_r1, max_errors, config.suffix_allow_indels)

    while current_back < len(dna_seq):
        search_start = current_back
        search_end = min(len(dna_seq), current_back + config.suffix_search_distance)

        if search_end - search_start < 5:
            break

        suffix_area = dna_seq[search_start:search_end]
        matches = list(re.finditer(query, suffix_area, flags=re.BESTMATCH))

        if not matches:
            break

        best_match = matches[0]

        r1_pos_global = search_start + best_match.start()
        actual_r1 = best_match.group()

        min_len = len(core_r1) - max_errors
        max_len = len(core_r1) + max_errors

        if min_len <= len(actual_r1) <= max_len:
            suffix_r1_list.append((r1_pos_global, actual_r1))
            current_back = r1_pos_global + len(actual_r1)
        else:
            current_back = r1_pos_global + 1

    suffix_r1_list.sort(key=lambda x: x[0])
    return suffix_r1_list
