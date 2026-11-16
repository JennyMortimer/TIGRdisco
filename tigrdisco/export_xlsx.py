"""Excel export of a finished run."""

import argparse
import hashlib
import json
from pathlib import Path

from . import __version__


def flatten(value, prefix=""):
    out = {}
    for key, item in value.items():
        key = prefix + key
        if isinstance(item, dict):
            out.update(flatten(item, key + "."))
        elif isinstance(item, list):
            out[key] = json.dumps(item, ensure_ascii=False)
        else:
            out[key] = item
    return out


def original_interval(record, start0, end0):
    if record["reverse_complemented"]:
        return record["record_length"] - end0 + 1, record["record_length"] - start0
    return start0 + 1, end0


def tables(records, kind):
    anchors, arrays, units = [], [], []
    for record in records:
        selected = (
            bool(record[kind]) if kind in ("SL", "DR") else not (record["SL"] or record["DR"])
        )
        if not selected:
            continue
        row = flatten({k: v for k, v in record.items() if k not in ("SL", "DR", "report_sl_stems")})
        row["SL"] = int(bool(record["SL"]))
        row["DR"] = int(bool(record["DR"]))
        row["architecture"] = (
            "SL+DR"
            if row["SL"] and row["DR"]
            else "SL" if row["SL"] else "DR" if row["DR"] else "neither"
        )
        row["anchor_input_start1"], row["anchor_input_end1"] = original_interval(
            record, record["CENTRAL_start0"], record["CENTRAL_end0"]
        )
        row["anchor_input_strand"] = "-" if record["reverse_complemented"] else "+"
        row["SL_array_count"] = len(record["SL"])
        row["DR_array_count"] = len(record["DR"])
        anchors.append(row)
        for ai, array in enumerate(record.get(kind, []), 1):
            identity = {"locus": record["locus"], "engine": kind, "array_index": ai}
            arrays.append(
                {
                    **identity,
                    **flatten({k: v for k, v in array.items() if k != "units"}),
                    "exported_unit_count": len(array["units"]),
                }
            )
            for ui, unit in enumerate(array["units"], 1):
                detail = {**identity, "unit_index": ui, **flatten(unit)}
                if kind == "DR":
                    sequence = "".join(unit[k] for k in ("R1", "spacer_a", "R2", "spacer_b"))
                    start = unit["contig_index"]
                    end = start + len(sequence)
                else:
                    if array.get("strand", "+") != "+":
                        raise ValueError(
                            "Unsupported SL window strand; refusing ambiguous coordinate conversion"
                        )
                    sequence = "".join(
                        unit[k]
                        for k in (
                            "box_C",
                            "spacer_A",
                            "box_D",
                            "connector",
                            "box_C_prime",
                            "spacer_B",
                            "box_D_prime",
                        )
                    )
                    start = record["SL_window_offset0"] + unit["start"] - 1
                    end = record["SL_window_offset0"] + unit["end"]
                detail["unit_sequence_oriented"] = sequence
                detail["input_start1"], detail["input_end1"] = original_interval(record, start, end)
                detail["input_strand"] = row["anchor_input_strand"]
                units.append(detail)
    return anchors, arrays, units


def write_table(book, name, rows, fallback):
    columns = list(dict.fromkeys(key for row in rows for key in row)) or fallback
    for part, offset in enumerate(range(0, max(1, len(rows)), 1048575), 1):
        chunk = rows[offset : offset + 1048575]
        sheet = book.add_worksheet(name if part == 1 else f"{name}_{part}")
        header = book.add_format({"bold": True, "bg_color": "#214F70", "font_color": "white"})
        sheet.write_row(0, 0, columns, header)
        sheet.freeze_panes(1, 0)
        sheet.set_column(0, len(columns) - 1, 20)
        sheet.set_column(0, 0, 55)
        for ri, row in enumerate(chunk, 1):
            for ci, key in enumerate(columns):
                value = row.get(key)
                if isinstance(value, str) and len(value) > 32767:
                    raise ValueError(f"Excel cell limit exceeded: {name}/{key}; use original JSONL")
                sheet.write(ri, ci, value)
        sheet.autofilter(0, 0, len(chunk), len(columns) - 1)


def export_run(run, out, style="template", workers=2):
    import xlsxwriter

    run, out = Path(run), Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new/empty XLSX output directory")
    with (run / "loci.jsonl").open(encoding="utf-8") as handle:
        records = [json.loads(line) for line in handle if line.strip()]
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    counts = {
        "N": len(records),
        "SL": sum(bool(r["SL"]) for r in records),
        "DR": sum(bool(r["DR"]) for r in records),
        "both": sum(bool(r["SL"]) and bool(r["DR"]) for r in records),
        "neither": sum(not (r["SL"] or r["DR"]) for r in records),
    }
    for key, value in counts.items():
        if manifest.get(key) != value:
            raise ValueError(f"Manifest count mismatch: {key}")
    if style not in ("template", "tables") or workers < 1:
        raise ValueError("Invalid XLSX style or workers")
    labels, recovery = {}, {}
    if style == "template":
        from .template_report import prepare, write_main

        labels = prepare(records, workers)
    out.mkdir(parents=True, exist_ok=True)
    notes = {
        "source_run": str(run.resolve()),
        "exporter_version": __version__,
        "classification": "Post-arbitration results from this run; both-positive anchors appear in both positive workbooks.",
        "negative_meaning": "Neither architecture detected under this run configuration; not evidence of biological absence.",
        "coordinates": "input_start1/end1 and anchor_input_start1/end1: original input record, 1-based inclusive. Native fields retain JSONL conventions.",
        "sequences": "Unit and protein sequences follow the Tas anchor orientation, including reverse-strand anchors.",
        "protein_sequence": "Present for new FASTA/GenBank runs. Empty or absent in JSONL from a run that did not store it; never inferred by exporter.",
        "annotations": "Detection scores preserved. Template mode adds HMM domain diagrams and NTD profile annotations; tables mode does not. TasA-like means neither NTD profile passed, not demonstrated domain absence.",
    }
    notes["xlsx_style"] = style
    notes["logo_method"] = (
        "Local DNA information-content renderer, 2-H-3/(2 ln(2) n); DR dominant length Core group >=4, SL template has no logo columns or images; exclusions listed for DR arrays. Not external WebLogo software."
    )
    notes["recovery"] = json.dumps(recovery)
    images = {}
    import tempfile

    for kind, filename in [
        ("SL", "sl_tas_with_tigr.xlsx"),
        ("DR", "dr_tas_with_tigr.xlsx"),
        ("neither", "tas_without_tigr.xlsx"),
    ]:
        anchors, arrays, units = tables(records, kind)
        with (
            tempfile.TemporaryDirectory(prefix="tigrdisco_report_") as tmp,
            xlsxwriter.Workbook(
                str(out / filename),
                {"constant_memory": True, "strings_to_formulas": False, "strings_to_urls": False},
            ) as book,
        ):
            if style == "template":
                images[kind] = write_main(book, records, kind, labels, Path(tmp))
            write_table(book, "Anchors", anchors, ["locus", "SL", "DR", "architecture"])
            if kind != "neither":
                write_table(book, "Arrays", arrays, ["locus", "engine", "array_index"])
                write_table(book, "Units", units, ["locus", "engine", "array_index", "unit_index"])
            write_table(
                book,
                "Run_info",
                [
                    {"field": k, "value": v}
                    for k, v in {
                        **notes,
                        **{"manifest." + k: v for k, v in flatten(manifest).items()},
                    }.items()
                ],
                ["field", "value"],
            )
    audit = {
        "exporter_version": __version__,
        "counts": counts,
        "xlsx_style": style,
        "image_counts": images,
        "recovery": recovery,
        "source_sha256": {
            name: hashlib.sha256((run / name).read_bytes()).hexdigest()
            for name in ("loci.jsonl", "manifest.json")
        },
    }
    (out / "export_manifest.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        type=Path,
        required=True,
        help="result directory of tigrdisco-fasta or tigrdisco-db",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="new or empty output directory",
    )
    parser.add_argument("--style", choices=("template", "tables"), default="template")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    print(
        json.dumps(
            export_run(args.run, args.out, style=args.style, workers=args.workers),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
