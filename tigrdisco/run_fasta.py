"""FASTA mode: find Tas genes and run both engines."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from . import __version__
from .sl_tigr_engine import detect
from .tas_tigr_engine.tigr_detection.detector import TigrDetectionConfig
from .dr_tigr_engine import arbitration
from .run_database import work, file_hash

_COMPLEMENT = str.maketrans("ACGTNRYKMSWBDHVacgtnrykmswbdhv", "TGCANYRMKSWVHDBtgcanyrmkswvhdb")
TIERS = ("high", "broad")


def read_fasta(path: Path):
    name = None
    chunks: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, _check(name, "".join(chunks), path)
                token = line[1:].split()[0] if line[1:].strip() else ""
                if not token:
                    raise ValueError(f"{path}:{number}: empty FASTA identifier")
                if token in seen:
                    raise ValueError(f"{path}:{number}: duplicate FASTA identifier {token!r}")
                seen.add(token)
                name, chunks = token, []
            else:
                if name is None:
                    raise ValueError(f"{path}:{number}: sequence before the first header")
                chunks.append(line)
    if name is None:
        raise ValueError(f"No FASTA records in {path}")
    yield name, _check(name, "".join(chunks), path)


def _check(name: str, sequence: str, path: Path) -> str:
    sequence = sequence.upper().replace("U", "T")
    invalid = set(sequence) - set("ACGTNRYKMSWBDHV")
    if invalid:
        raise ValueError(
            f"{path}: record {name} contains invalid nucleotides: {''.join(sorted(invalid))}"
        )
    if not sequence:
        raise ValueError(f"{path}: record {name} is empty")
    return sequence


def call_orfs(records, mode: str, minimum_aa: int):
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=(mode == "meta"))
    if mode == "single":
        joined = "TTAATTAATTAA".join(sequence for _, sequence in records)
        finder.train(joined)
    out = []
    for contig, sequence in records:
        for gene in finder.find_genes(sequence):
            protein = gene.translate().rstrip("*")
            if len(protein) < minimum_aa or "*" in protein:
                continue
            out.append((contig, gene.begin - 1, gene.end, 1 if gene.strand > 0 else -1, protein))
    return out


def screen(orfs, hmm_path: Path, tier: str, cpus: int):
    import pyhmmer
    from .tas_tigr_engine.models import HMMThresholds

    limits = (
        HMMThresholds.high_confidence_tas() if tier == "high" else HMMThresholds.broad_tas_screen()
    )
    alphabet = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            pyhmmer.easel.TextSequence(name=f"orf{index}".encode(), sequence=protein).digitize(
                alphabet
            )
            for index, (_, _, _, _, protein) in enumerate(orfs)
        ],
    )
    with pyhmmer.plan7.HMMFile(hmm_path) as handle:
        profiles = list(handle)
    if len(profiles) != 1:
        raise ValueError("The Tas-core screen requires exactly one HMM profile")
    model = profiles[0]
    kept = {}
    for top in pyhmmer.hmmsearch([model], block, cpus=cpus, E=10.0, domE=10.0):
        for hit in top:
            name = hit.name if isinstance(hit.name, str) else hit.name.decode()
            index = int(name[3:])
            protein = orfs[index][4]
            domains = [d for d in hit.domains]
            if not domains:
                continue
            best = max(domains, key=lambda d: d.score)
            hmm_cov = (
                _coverage((d.alignment.hmm_from, d.alignment.hmm_to) for d in domains) / model.M
            )
            prot_cov = _coverage(
                (d.alignment.target_from, d.alignment.target_to) for d in domains
            ) / len(protein)
            if hit.evalue > limits.max_full_evalue or best.i_evalue > limits.max_domain_evalue:
                continue
            if limits.min_bit_score is not None and hit.score < limits.min_bit_score:
                continue
            if hmm_cov < limits.min_hmm_coverage or prot_cov < limits.min_protein_coverage:
                continue
            strict = HMMThresholds.high_confidence_tas()
            is_high = (
                hit.evalue <= strict.max_full_evalue
                and best.i_evalue <= strict.max_domain_evalue
                and hit.score >= strict.min_bit_score
                and hmm_cov >= strict.min_hmm_coverage
                and prot_cov >= strict.min_protein_coverage
            )
            kept[index] = {
                "tier": "high_confidence_tas" if is_high else "possible_or_remote_tas",
                "full_E": float(hit.evalue),
                "score": float(hit.score),
                "domain_E": float(best.i_evalue),
                "hmm_intervals": sorted(
                    {(d.alignment.target_from, d.alignment.target_to) for d in domains}
                ),
                "hmm_model": "MITTas HMM",
                "HMM_cov": round(hmm_cov, 6),
                "protein_cov": round(prot_cov, 6),
            }
    return kept


def _coverage(intervals) -> int:
    total = 0
    start = end = None
    for a, b in sorted(intervals):
        if start is None:
            start, end = a, b
        elif a <= end + 1:
            end = max(end, b)
        else:
            total += end - start + 1
            start, end = a, b
    return total + (0 if start is None else end - start + 1)


def tasks(records_by_id, anchors, sl_flank: int, dr_flank: int):
    for contig, start0, end0, strand, protein, meta in anchors:
        sequence = records_by_id[contig]
        a, b = start0, end0
        if strand < 0:
            sequence = sequence.translate(_COMPLEMENT)[::-1]
            a, b = len(sequence) - end0, len(sequence) - start0
        window_start = max(0, a - sl_flank)
        window_end = min(len(sequence), b + sl_flank)
        flanks = [
            ("upstream", max(0, a - dr_flank), a),
            ("downstream", b, min(len(sequence), b + dr_flank)),
        ]
        locus = f"{contig}__{start0 + 1}_{end0}_{'+' if strand > 0 else '-'}"
        yield (
            locus,
            sequence[window_start:window_end],
            window_start,
            a,
            b,
            protein,
            [(side, sequence[x:y], x, y) for side, x, y in flanks],
            len(sequence),
            strand < 0,
        ), meta


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fasta", type=Path, help="nucleotide FASTA: contigs, genomes or pre-cut windows")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--xlsx-style",
        choices=("template", "tables"),
        default="template",
        help="template: with domain diagrams, logos and Tas subtype; tables: compact, no images",
    )
    p.add_argument(
        "--xlsx",
        action="store_true",
        help="also write three Excel workbooks to OUT/xlsx",
    )
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--sl-flank", type=int, default=1000)
    p.add_argument("--dr-flank", type=int, default=600)
    p.add_argument(
        "--orf-mode",
        choices=("meta", "single"),
        default="meta",
        help="pyrodigal mode; meta suits mixed or short input, single trains on the input",
    )
    p.add_argument("--min-protein", type=int, default=60, help="shortest ORF kept, in amino acids")
    p.add_argument(
        "--tier",
        choices=TIERS,
        default="broad",
        help="anchor screen: broad keeps possible/remote Tas, high keeps only high-confidence",
    )
    p.add_argument("--limit", type=int, help="stop after this many input records")
    p.add_argument("--reread", choices=("on", "off"), default="on")
    p.add_argument("--dr-candidate-lengths", type=int, nargs="+", default=list(range(8, 14)))
    p.add_argument("--dr-minimum-units", type=int, default=2)
    p.add_argument("--dr-minimum-core-units", type=int, default=4)
    p.add_argument("--dr-minimum-core-spacing", type=float, default=32.0)
    p.add_argument("--dr-maximum-core-spacing", type=float, default=40.0)
    p.add_argument("--dr-core-spacing-tolerance", type=float, default=3.0)
    p.add_argument("--dr-minimum-core-spacing-fraction", type=float, default=0.75)
    p.add_argument("--dr-frame-candidates", type=int, default=8)
    p.add_argument("--dr-minimum-bits", type=float, default=18.0)
    p.add_argument(
        "--arbitrate",
        choices=arbitration.MODES,
        default="covered",
        help="how a stem-loop array that overlaps a dual-repeat array is removed (default: covered)",
    )
    p.add_argument(
        "--arbitrate-coverage",
        type=float,
        default=1.0,
        help="covered: fraction of units that must overlap a dual-repeat unit",
    )
    p.add_argument("--arbitrate-gap", type=int, default=1200, help="distance: largest gap in bp")
    p.add_argument(
        "--arbitrate-stem",
        type=int,
        default=8,
        help="intersect, distance: keep arrays whose median stem reaches this; 0 disables",
    )
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    dr_config = TigrDetectionConfig(
        **{
            k: (
                tuple(args.dr_candidate_lengths)
                if k == "candidate_lengths"
                else getattr(args, "dr_" + k)
            )
            for k in TigrDetectionConfig.__dataclass_fields__
        }
    )
    sl_config = detect.Config(reread=args.reread == "on")
    numeric = [v for v in asdict(dr_config).values() if isinstance(v, (int, float))]
    if (
        any(not math.isfinite(v) or v < 0 for v in numeric)
        or min(args.dr_candidate_lengths) < 1
        or args.dr_minimum_core_units < 2
        or args.dr_frame_candidates < 1
        or args.dr_minimum_units < 1
        or not 0 <= args.dr_minimum_core_spacing_fraction <= 1
        or args.dr_minimum_core_spacing > args.dr_maximum_core_spacing
    ):
        parser.error("Invalid DR parameter range")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.arbitrate_gap < 0 or args.arbitrate_stem < 0:
        parser.error("Arbitration gap and stem must be non-negative")
    if not 0 < args.arbitrate_coverage <= 1:
        parser.error("--arbitrate-coverage must be greater than 0 and at most 1")
    arbitrate = {
        "mode": args.arbitrate,
        "gap": args.arbitrate_gap,
        "stem": args.arbitrate_stem,
        "coverage": args.arbitrate_coverage,
    }
    if args.min_protein < 20:
        parser.error("--min-protein must be at least 20")
    if args.workers < 1 or min(args.sl_flank, args.dr_flank) < 1:
        raise SystemExit("Workers and flank lengths must be positive")
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit("Use a new/empty output directory")

    started = time.time()
    records = list(read_fasta(args.fasta))
    if args.limit is not None:
        records = records[: args.limit]
    by_id = dict(records)
    orfs = call_orfs(records, args.orf_mode, args.min_protein)
    hmm_path = Path(__file__).resolve().parent / "data" / "MITTas.hmm3"
    passed = screen(orfs, hmm_path, args.tier, max(1, args.workers))
    anchors = [
        (orfs[i][0], orfs[i][1], orfs[i][2], orfs[i][3], orfs[i][4], meta)
        for i, meta in sorted(passed.items())
    ]
    print(
        f"{len(records)} records, {len(orfs)} ORFs, {len(anchors)} Tas anchors "
        f"({round(time.time() - started, 1)} s)",
        flush=True,
    )

    args.out.mkdir(parents=True, exist_ok=True)
    pairs = list(tasks(by_id, anchors, args.sl_flank, args.dr_flank))
    metas = [m for _, m in pairs]
    total = sl_count = dr_count = both = 0
    with (
        ProcessPoolExecutor(max_workers=args.workers) as pool,
        (args.out / "loci.jsonl").open("w") as jsonl,
        (args.out / "loci.csv").open("w", newline="") as handle,
    ):
        writer = csv.writer(handle)
        writer.writerow(
            [
                "locus",
                "contig",
                "anchor_start1",
                "anchor_end1",
                "anchor_strand",
                "tier",
                "hmm_score",
                "SL",
                "DR",
                "architecture",
            ]
        )
        runner = partial(work, dr_config=dr_config, sl_config=sl_config, arbitrate=arbitrate)
        for result, meta, (task, _) in zip(
            pool.map(runner, [t for t, _ in pairs], chunksize=4), metas, pairs
        ):
            result["protein_sequence"] = task[5]
            result["anchor"] = meta
            result["contig"] = task[0].rsplit("__", 1)[0]
            jsonl.write(json.dumps(result) + "\n")
            sl = bool(result["SL"])
            dr = bool(result["DR"])
            total += 1
            sl_count += sl
            dr_count += dr
            both += sl and dr
            label = "SL+DR" if sl and dr else "SL" if sl else "DR" if dr else "neither"
            start1, end1, strand = task[0].rsplit("__", 1)[1].rsplit("_", 2)
            writer.writerow(
                [
                    result["locus"],
                    result["contig"],
                    start1,
                    end1,
                    strand,
                    meta["tier"],
                    f"{meta['score']:.1f}",
                    int(sl),
                    int(dr),
                    label,
                ]
            )
            if total % 500 == 0:
                print(
                    total,
                    "anchors",
                    sl_count,
                    "SL",
                    dr_count,
                    "DR",
                    round(time.time() - started, 1),
                    "seconds",
                    flush=True,
                )
    manifest = {
        "N": total,
        "SL": sl_count,
        "DR": dr_count,
        "both": both,
        "neither": total - sl_count - dr_count + both,
        "input_records": len(records),
        "orfs_called": len(orfs),
        "version": __version__,
        "seconds": time.time() - started,
        "input_sha256": file_hash(args.fasta),
        "sl_config": asdict(sl_config),
        "dr_config": asdict(dr_config),
        "sl_flank": args.sl_flank,
        "dr_flank": args.dr_flank,
        "arbitration": arbitrate,
        "anchor_mode": f"pyrodigal {args.orf_mode} ORF calling, Tas-core HMM screen at tier "
        f"'{args.tier}', minimum {args.min_protein} aa",
        "limit": args.limit,
        "dr_source": "tas_tigr_engine",
    }
    manifest["anchor_config"] = {
        "orf_mode": args.orf_mode,
        "min_protein": args.min_protein,
        "tier": args.tier,
    }
    manifest["xlsx_style"] = args.xlsx_style
    manifest["xlsx_requested"] = args.xlsx
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    if args.xlsx:
        from .export_xlsx import export_run

        export_run(args.out, args.out / "xlsx", style=args.xlsx_style, workers=args.workers)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
