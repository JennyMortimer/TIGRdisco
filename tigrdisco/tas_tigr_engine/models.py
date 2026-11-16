"""Data models of the dual-repeat engine."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PredictedProtein:
    protein_id: str
    contig_id: str
    sequence: str
    start: int
    end: int
    strand: str
    partial_5prime: bool
    partial_3prime: bool

    def __post_init__(self) -> None:
        if self.start < 0 or self.end <= self.start:
            raise ValueError(f"Invalid interval for {self.protein_id}: {self.start}:{self.end}")
        if self.strand not in {"+", "-"}:
            raise ValueError(f"Invalid strand for {self.protein_id}: {self.strand}")
        if not self.sequence:
            raise ValueError(f"Empty protein sequence for {self.protein_id}")


@dataclass(frozen=True)
class HMMThresholds:
    use_ga_cutoff: bool = False
    max_full_evalue: float = 1e-15
    max_domain_evalue: float = 1e-8
    min_bit_score: float | None = 50.0
    min_hmm_coverage: float = 0.40
    min_protein_coverage: float = 0.45

    @classmethod
    def high_confidence_tas(cls) -> "HMMThresholds":
        return cls()

    @classmethod
    def broad_tas_screen(cls) -> "HMMThresholds":
        return cls(
            max_full_evalue=1e-5,
            max_domain_evalue=1.0,
            min_bit_score=5.0,
            min_hmm_coverage=0.20,
            min_protein_coverage=0.08,
        )


@dataclass(frozen=True)
class FlankRegion:
    side: str
    sequence: str
    oriented_start: int
    oriented_end: int
    contig_start: int
    contig_end: int
    truncated: bool


@dataclass(frozen=True)
class TigrUnit:
    contig_index: int
    R1: str
    spacer_a: str
    R2: str
    spacer_b: str
    classification: str = "Core"


@dataclass(frozen=True)
class TigrArray:
    side: str
    R1_consensus: str
    R2_consensus: str
    units: tuple[TigrUnit, ...]
    flank_contig_start: int
    flank_contig_end: int
    bits: float = 0.0
