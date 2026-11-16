"""GenBank mode: run both engines on records with a CENTRAL anchor."""

import argparse
import csv
import hashlib
import json
import re
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from functools import partial
from pathlib import Path

from Bio.Seq import Seq

from . import __version__
from .dr_tigr_engine import arbitration
from .sl_tigr_engine import detect
from .tas_tigr_engine.models import FlankRegion, PredictedProtein
from .tas_tigr_engine.tigr_detection.detector import TigrDetectionConfig, detect_tigr


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            digest.update(chunk)
    return digest.hexdigest()


def records(path):
    buf = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            buf.append(line)
            if line.strip() == "//":
                yield "".join(buf)
                buf = []
    if "".join(buf).strip():
        raise ValueError("Incomplete final GenBank record")


def prepare(path, sl_flank, dr_flank, limit):
    seen = set()
    for i, block in enumerate(records(path)):
        if limit is not None and i >= limit:
            break
        name = re.search(r"^LOCUS\s+(\S+)\s+(\d+)", block, re.M)
        central = re.findall(r"^     CENTRAL\s+(complement\()?<?(\d+)\.\.>?(\d+)\)?", block, re.M)
        if not name or len(central) != 1 or "ORIGIN" not in block:
            raise ValueError("Expected exactly one CENTRAL CDS and sequence per record")
        cid = name[1]
        seq = "".join(
            re.findall("[A-Za-z]+", block.split("ORIGIN", 1)[1].split("//", 1)[0])
        ).upper()
        if cid in seen or len(seq) != int(name[2]):
            raise ValueError(f"Duplicate ID or length mismatch: {cid}")
        seen.add(cid)
        strand, a, b = central[0]
        a = int(a) - 1
        b = int(b)
        if not 0 <= a < b <= len(seq):
            raise ValueError(f"Invalid CENTRAL: {cid}")
        if strand:
            seq = str(Seq(seq).reverse_complement())
            a, b = len(seq) - b, len(seq) - a
        ws = max(0, a - sl_flank)
        we = min(len(seq), b + sl_flank)
        nt = seq[a:b]
        aa = str(Seq(nt[: len(nt) // 3 * 3]).translate(table=11)).rstrip("*").replace("*", "X")
        flanks = [
            ("upstream", max(0, a - dr_flank), a),
            ("downstream", b, min(len(seq), b + dr_flank)),
        ]
        yield (
            cid,
            seq[ws:we],
            ws,
            a,
            b,
            aa,
            [(side, seq[x:y], x, y) for side, x, y in flanks],
            len(seq),
            bool(strand),
        )


def work(task, dr_config=None, sl_config=None, arbitrate=None):
    cid, window, offset, a, b, aa, flanks, L, rc = task
    sl = detect.process((f"{cid}|{offset}", window, sl_config or detect.Config()))
    protein = PredictedProtein(cid + "|CENTRAL", cid, aa, a, b, "+", False, False)
    dr = []
    for side, s, x, y in flanks:
        f = FlankRegion(side, s, x, y, x, y, len(s) < 100)
        arr = detect_tigr(f, protein, config=dr_config)
        if arr is not None:
            dr.append(asdict(arr))
    if arbitrate and arbitrate.get("mode", "off") != "off":
        sl = arbitration.arbitrate_locus(sl, dr, offset, **arbitrate)
    return {
        "locus": cid,
        "protein_sequence": aa,
        "record_length": L,
        "reverse_complemented": rc,
        "coordinate_system": "oriented record; DR zero-based half-open model, SL one-based window core",
        "SL_window_offset0": offset,
        "CENTRAL_start0": a,
        "CENTRAL_end0": b,
        "SL": sl,
        "DR": dr,
        "report_sl_stems": {
            f"{u['start']}:{u['end']}": {
                "left": window[u["extended_start"] - 1 : u["start"] - 1],
                "right": window[u["end"] : u["extended_end"]],
            }
            for ar in sl
            for u in ar["units"]
        },
    }


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genbank", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--sl-flank", type=int, default=1000)
    p.add_argument("--dr-flank", type=int, default=600)
    p.add_argument("--limit", type=int)
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


def main():
    p = build_parser()
    args = p.parse_args()
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
    import math

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
        p.error("Invalid DR parameter range")
    if args.limit is not None and args.limit < 1:
        p.error("--limit must be positive")
    if args.arbitrate_gap < 0 or args.arbitrate_stem < 0:
        p.error("Arbitration gap and stem must be non-negative")
    if not 0 < args.arbitrate_coverage <= 1:
        p.error("--arbitrate-coverage must be greater than 0 and at most 1")
    arbitrate = {
        "mode": args.arbitrate,
        "gap": args.arbitrate_gap,
        "stem": args.arbitrate_stem,
        "coverage": args.arbitrate_coverage,
    }
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit("Use a new/empty output directory")
    if args.workers < 1 or min(args.sl_flank, args.dr_flank) < 1:
        raise SystemExit("Workers and flank lengths must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    N = ns = nd = nb = 0
    with (
        ProcessPoolExecutor(max_workers=args.workers) as pool,
        (args.out / "loci.jsonl").open("w") as f,
        (args.out / "loci.csv").open("w", newline="") as c,
    ):
        w = csv.writer(c)
        w.writerow(["locus", "SL", "DR", "architecture"])
        for r in pool.map(
            partial(work, dr_config=dr_config, sl_config=sl_config, arbitrate=arbitrate),
            prepare(args.genbank, args.sl_flank, args.dr_flank, args.limit),
            chunksize=4,
        ):
            f.write(json.dumps(r) + "\n")
            s = bool(r["SL"])
            d = bool(r["DR"])
            N += 1
            ns += s
            nd += d
            nb += s and d
            label = "SL+DR" if s and d else "SL" if s else "DR" if d else "neither"
            w.writerow([r["locus"], int(s), int(d), label])
            if N % 500 == 0:
                print(
                    N, "loci", ns, "SL", nd, "DR", round(time.time() - t0, 1), "seconds", flush=True
                )
    m = {
        "N": N,
        "SL": ns,
        "DR": nd,
        "both": nb,
        "neither": N - ns - nd + nb,
        "version": __version__,
        "seconds": time.time() - t0,
        "input_sha256": file_hash(args.genbank),
        "sl_config": asdict(sl_config),
        "dr_config": asdict(dr_config),
        "sl_flank": args.sl_flank,
        "dr_flank": args.dr_flank,
        "arbitration": arbitrate,
        "anchor_mode": "CENTRAL annotations supplied, no Prodigal/HMM filtering performed",
        "limit": args.limit,
        "dr_source": "tas_tigr_engine",
    }
    m["xlsx_requested"] = args.xlsx
    m["xlsx_style"] = args.xlsx_style
    (args.out / "manifest.json").write_text(json.dumps(m, indent=2))
    print(json.dumps(m, indent=2))
    if args.xlsx:
        from .export_xlsx import export_run

        export_run(args.out, args.out / "xlsx", style=args.xlsx_style, workers=args.workers)


if __name__ == "__main__":
    main()
