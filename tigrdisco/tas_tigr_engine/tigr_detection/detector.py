"""Dual-repeat array detection."""

from dataclasses import dataclass
from statistics import median

from .weblogo import array_bits, candidate_frames
from .extractor import TigrConfig as LegacyConfig, extract_tigr_array

from ..models import FlankRegion, PredictedProtein, TigrArray, TigrUnit


@dataclass(frozen=True)
class TigrDetectionConfig:
    candidate_lengths: tuple[int, ...] = tuple(range(8, 14))
    minimum_units: int = 2
    minimum_core_units: int = 4
    minimum_core_spacing: float = 32.0
    maximum_core_spacing: float = 40.0
    core_spacing_tolerance: float = 3.0
    minimum_core_spacing_fraction: float = 0.75
    frame_candidates: int = 8
    minimum_bits: float = 18.0


def _legacy_config(config: TigrDetectionConfig) -> LegacyConfig:
    return LegacyConfig(
        candidate_lengths=config.candidate_lengths,
        discovery_mismatches=2,
        discovery_allow_indels=False,
        window_max_errors=4,
        window_allow_indels=True,
        unit_max_errors=3,
        unit_allow_indels=False,
        reorganize_misaligned_unit=True,
        reorganize_unit_length=36,
        reorganize_r2_max_mismatches=4,
        prefix_max_errors=2,
        prefix_allow_indels=False,
        prefix_search_distance=40,
        prefix_nearest_distance=20,
        suffix_max_errors=3,
        suffix_allow_indels=False,
        suffix_search_distance=40,
    )


def _contig_index(flank: FlankRegion, protein: PredictedProtein, unit: dict[str, object]) -> int:
    local_start = int(unit["start"])
    unit_length = len(_unit_sequence(unit))
    if protein.strand == "+":
        return flank.contig_start + local_start
    return flank.contig_end - (local_start + unit_length)


def _passes_periodic_unit_spacing(
    raw_units: list[dict[str, object]],
    config: TigrDetectionConfig,
) -> bool:
    positions = sorted(int(unit["start"]) for unit in raw_units)
    window_size = config.minimum_core_units
    if len(positions) < window_size:
        return False
    for window_start in range(len(positions) - window_size + 1):
        window = positions[window_start : window_start + window_size]
        gaps = [window[index + 1] - window[index] for index in range(len(window) - 1)]
        median_spacing = median(gaps)
        if not (config.minimum_core_spacing <= median_spacing <= config.maximum_core_spacing):
            continue
        within_tolerance = sum(
            abs(gap - median_spacing) <= config.core_spacing_tolerance for gap in gaps
        )
        if within_tolerance / len(gaps) >= config.minimum_core_spacing_fraction:
            return True
    return False


def _unit_sequence(unit: dict[str, object]) -> str:
    return "".join(
        "" if str(unit[key]) == "(Empty)" else str(unit[key])
        for key in ("R1", "spacer_a", "R2", "spacer_b")
    )


def _clean_units(units) -> list[dict[str, object]]:
    kept = [u for u in units if _unit_sequence(u) and not (set(_unit_sequence(u)) - set("ACGT"))]
    kept.sort(key=lambda unit: int(unit["start"]))
    return kept


def detect_tigr(
    flank: FlankRegion,
    protein: PredictedProtein,
    config: TigrDetectionConfig | None = None,
) -> TigrArray | None:
    config = config or TigrDetectionConfig()
    if flank.truncated or not flank.sequence:
        return None
    legacy = _legacy_config(config)
    best = None
    for R1, R2 in candidate_frames(flank.sequence, legacy, config.frame_candidates):
        result = extract_tigr_array(flank.sequence, R1, R2, legacy)
        if not result or result["total_count"] < config.minimum_units:
            continue
        units = _clean_units(result["units"])
        if len(units) < config.minimum_core_units:
            continue
        if not _passes_periodic_unit_spacing(units, config):
            continue
        bits = array_bits(units)
        if best is None or bits > best[0]:
            best = (bits, R1, R2, units)
    if best is None or best[0] < config.minimum_bits:
        return None
    bits, R1, R2, raw_units = best

    units = tuple(
        TigrUnit(
            contig_index=_contig_index(flank, protein, unit),
            R1="" if str(unit["R1"]) == "(Empty)" else str(unit["R1"]),
            spacer_a="" if str(unit["spacer_a"]) == "(Empty)" else str(unit["spacer_a"]),
            R2="" if str(unit["R2"]) == "(Empty)" else str(unit["R2"]),
            spacer_b="" if str(unit["spacer_b"]) == "(Empty)" else str(unit["spacer_b"]),
            classification=str(unit["classification"]),
        )
        for unit in raw_units
    )
    return TigrArray(
        side=flank.side,
        R1_consensus=R1,
        R2_consensus=R2,
        units=units,
        flank_contig_start=flank.contig_start,
        flank_contig_end=flank.contig_end,
        bits=round(bits, 3),
    )
