"""Template-style Excel reports with domain diagrams and sequence logos."""

from collections import defaultdict
from functools import lru_cache

from PIL import Image, ImageDraw
from . import report
from .shared import hmm_diagram, hmm_search, ntd

DR_KEYS = ("R1", "spacer_a", "R2", "spacer_b")
SL_KEYS = ("box_C", "spacer_A", "box_D", "connector", "box_C_prime", "spacer_B", "box_D_prime")


def add_missing_hmm(records, workers):
    import pyhmmer

    needed = [r for r in records if "hmm_intervals" not in r.get("anchor", {})]
    if not needed:
        return
    alphabet = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            pyhmmer.easel.TextSequence(
                name=f"p{i}".encode(), sequence=r["protein_sequence"]
            ).digitize(alphabet)
            for i, r in enumerate(needed)
        ],
    )
    with pyhmmer.plan7.HMMFile(hmm_search.DATA) as f:
        models = list(f)
    for top in pyhmmer.hmmsearch(models, block, cpus=workers, E=10.0, domE=10.0):
        for hit in top:
            name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            r = needed[int(name[1:])]
            domains = list(hit.domains)
            if not domains:
                continue
            intervals = sorted({(d.alignment.target_from, d.alignment.target_to) for d in domains})
            meta = {
                "score": float(hit.score),
                "full_E": float(hit.evalue),
                "domain_E": float(min(d.i_evalue for d in domains)),
                "protein_cov": hmm_search.union_length(intervals) / len(r["protein_sequence"]),
                "HMM_cov": hmm_search.union_length(
                    (d.alignment.hmm_from, d.alignment.hmm_to)
                    for d in sorted(domains, key=lambda d: d.alignment.hmm_from)
                )
                / models[0].M,
                "hmm_intervals": intervals,
                "hmm_model": "MITTas HMM",
                "tier": "supplied_anchor_report_only",
            }
            r.setdefault("anchor", {}).update(
                {k: v for k, v in meta.items() if k not in r.get("anchor", {})}
            )
    for r in needed:
        r.setdefault("anchor", {}).setdefault("tier", "supplied_anchor_no_reportable_hmm")
        r["anchor"].setdefault("hmm_intervals", [])


def array_text(record, kind):
    blocks = []
    for ai, array in enumerate(record[kind], 1):
        units = array["units"]
        if kind == "DR":
            keys = DR_KEYS
            heading = f"Array={ai}; Side={array['side']}; IC={array['bits']:.1f} bits"
            labels = ["Idx0", "R1", "Spacer A", "R2", "Spacer B", "Class"]
            data = [
                [str(u["contig_index"]), *[u[k] for k in keys], u["classification"]] for u in units
            ]
        else:
            a = min(u["start"] for u in units) - 1 + record["SL_window_offset0"]
            b = max(u["end"] for u in units) + record["SL_window_offset0"]
            side = (
                "upstream"
                if b <= record["CENTRAL_start0"]
                else "downstream" if a >= record["CENTRAL_end0"] else "overlapping_CDS"
            )
            heading = f"Array={ai}; Side={side}; C={array['box_C']}; D={array['box_D']}; class={array.get('array_class','')}; support={array.get('support_fraction','')}"
            labels = [
                "Idx1",
                "L-stem",
                "box C",
                "Spacer A",
                "box D",
                "Connector",
                "box C'",
                "Spacer B",
                "box D'",
                "R-stem",
            ]
            data = []
            for u in units:
                stem = record.get("report_sl_stems", {}).get(f"{u['start']}:{u['end']}", {})
                data.append(
                    [
                        str(u["start"] + record["SL_window_offset0"]),
                        stem.get("left", "NA"),
                        *[u[k] for k in SL_KEYS],
                        stem.get("right", "NA"),
                    ]
                )
        width = [max(len(row[i]) for row in [labels] + data) for i in range(len(labels))]
        lines = [
            " | ".join(value.ljust(width[i]) for i, value in enumerate(row))
            for row in [labels] + data
        ]
        blocks.append(
            "\n".join(
                [
                    heading,
                    "Coordinates: Tas-oriented record (see Units for original input coordinates)",
                    *lines,
                    f"{len(units)} units in total",
                ]
            )
        )
    return "\n\n".join(blocks)


def logo_group(array, kind):
    groups = defaultdict(list)
    keys = DR_KEYS if kind == "DR" else SL_KEYS
    eligible_units = [u for u in array["units"] if kind != "DR" or u["classification"] == "Core"]
    for u in eligible_units:
        seq = "".join(u[k] for k in keys)
        if seq and not (set(seq) - set("ACGT")):
            signature = (len(seq),) if kind == "DR" else tuple(len(u[k]) for k in keys)
            groups[signature].append(seq)
    minimum = 4 if kind == "DR" else 3
    eligible = [k for k, v in groups.items() if len(v) >= minimum]
    if not eligible and kind == "SL":
        box_keys = ("box_C", "box_D", "box_C_prime", "box_D_prime")
        box_groups = defaultdict(list)
        for u in eligible_units:
            seq = "".join(u[k] for k in box_keys)
            if seq and not (set(seq) - set("ACGT")):
                box_groups[tuple(len(u[k]) for k in box_keys)].append(seq)
        choices = [k for k, v in box_groups.items() if len(v) >= minimum]
        if choices:
            key = min(choices, key=lambda k: (-len(box_groups[k]), k))
            return (
                box_groups[key],
                f'BOX-ONLY C/D/C-prime/D-prime; used={len(box_groups[key])}/{len(array["units"])}; widths={key}; variable spacers excluded',
            )
    if not eligible:
        return (
            [],
            f'No logo: fewer than {minimum} unambiguous units in a compatible group; total={len(array["units"])}',
        )
    key = min(eligible, key=lambda k: (-len(groups[k]), k))
    selected = groups[key]
    return (
        selected,
        f'used={len(selected)}/{len(array["units"])}; length={len(selected[0])}; group={key}; '
        + ("Core only" if kind == "DR" else "matched segment lengths; mismatches retained"),
    )


@lru_cache(maxsize=2048)
def logo_glyph(base, width, height):
    colors = {"A": "#FF9900", "C": "#0000FF", "G": "#0000FF", "T": "#FF6600"}
    tile = Image.new("RGBA", (160, 160))
    ImageDraw.Draw(tile).text((5, 0), base, font=hmm_diagram._font("bold", 110), fill=colors[base])
    return tile.crop(tile.getbbox()).resize((width, height), Image.Resampling.LANCZOS)


def render_logo(sequences, output):
    from .tas_tigr_engine.tigr_detection.weblogo import information_stacks

    width, height = max(700, 32 * len(sequences[0]) + 85), 180
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    left, right, top, bottom = 52, width - 15, 14, 135
    step = (right - left) / len(sequences[0])
    for i, stack in enumerate(information_stacks(sequences)):
        y = bottom
        for base, bits in stack:
            if bits <= 0:
                continue
            h = max(1, round(bits / 2 * (bottom - top)))
            tile = logo_glyph(base, max(1, int(step * 0.94)), h)
            y -= h
            image.paste(tile, (round(left + i * step), y), tile)
    draw.line((left, top, left, bottom, right, bottom), fill="black", width=2)
    for value in (0, 1, 2):
        draw.text(
            (15, bottom - value * (bottom - top) / 2 - 7),
            str(value),
            font=hmm_diagram._font("regular", 14),
            fill="black",
        )
    draw.text(
        (left, 163),
        "Information (bits)",
        font=hmm_diagram._font("regular", 12),
        fill="black",
        anchor="la",
    )
    for pos in range(5, len(sequences[0]) + 1, 5):
        draw.text(
            (left + (pos - 0.5) * step, bottom + 8),
            str(pos),
            anchor="mt",
            font=hmm_diagram._font("regular", 13),
            fill="black",
        )
    draw.text(
        (right, 163),
        f"Sequence logo | n={len(sequences)}",
        anchor="ra",
        font=hmm_diagram._font("regular", 12),
        fill="black",
    )
    image.save(output)


def prepare(records, workers=2):
    if any(not r.get("protein_sequence") for r in records):
        raise ValueError("Template XLSX needs protein sequences; rerun with --xlsx.")
    add_missing_hmm(records, workers)
    print(f"Annotating NTD profiles for {len(records)} anchors (report only)...", flush=True)
    labels = ntd.annotate([(r["locus"], r["protein_sequence"]) for r in records])
    return labels


def write_main(book, records, kind, labels, tmp):
    selected = [
        r for r in records if (bool(r[kind]) if kind in ("SL", "DR") else not (r["SL"] or r["DR"]))
    ]
    fields = (
        report.WITHOUT_FIELDS
        if kind == "neither"
        else (
            report.WITH_FIELDS
            if kind == "SL"
            else report.WITH_FIELDS + ("weblogo_array_usage", "weblogo")
        )
    )
    sheet = book.add_worksheet(
        {"SL": "sl-tigr_engine", "DR": "Tas with TIGR", "neither": "Tas without TIGR"}[kind]
    )
    lengths = (
        sorted(max(map(len, array_text(r, kind).splitlines()), default=0) for r in selected)
        if kind != "neither"
        else []
    )
    array_width = (
        min(255, lengths[min(len(lengths) - 1, int(0.99 * len(lengths)))] + 3) if lengths else 80
    )
    fmt = report._formats(book)
    report._configure(sheet, fields, fmt, array_width)
    if kind == "DR":
        sheet.set_column(13, 13, 36, fmt["body"])
        sheet.set_column(14, 14, 80, fmt["body"])
    image_counts = {"hmm": 0, "logo": 0, "missing_hmm": 0}
    for index, r in enumerate(selected, 1):
        anchor = r.get("anchor", {})
        row = {
            "protein_name": r["locus"],
            "sequence_name": r.get("contig", r["locus"]),
            "strand": "-" if r["reverse_complemented"] else "+",
            "protein_sequence": r["protein_sequence"],
            **anchor,
            "ntd_annotation": labels[r["locus"]],
        }
        text = array_text(r, kind) if kind != "neither" else ""
        row["tigr_array"] = text
        logos = []
        usage = []
        if kind == "DR":
            for ai, array in enumerate(r[kind], 1):
                sequences, explanation = logo_group(array, kind)
                usage.append(f"Array {ai}: {explanation}")
                if sequences:
                    logo = tmp / f"{kind}_{index}_{ai}.png"
                    render_logo(sequences, logo)
                    logos.append((ai, logo))
                    image_counts["logo"] += 1
            row["weblogo_array_usage"] = "\n".join(usage)
        for col, field in enumerate(fields):
            value = row.get(field, "")
            if isinstance(value, str) and len(value) > 32767:
                if field != "tigr_array":
                    raise ValueError(f'Excel cell too long: {r["locus"]} {field}')
                value = (
                    "Array text exceeds Excel cell capacity. See complete Arrays and Units sheets."
                )
            if field in report.NUMERIC:
                report._write_numeric(sheet, index, col, value)
            else:
                sheet.write(index, col, value)
        if kind != "neither":
            intervals = anchor.get("hmm_intervals", [])
            if intervals:
                png = tmp / f"{kind}_{index}_hmm.png"
                hmm_diagram.render(
                    len(r["protein_sequence"]),
                    intervals,
                    hmm_search.union_length(intervals) / len(r["protein_sequence"]),
                    anchor.get("hmm_model", "MITTas HMM"),
                    png,
                )
                sheet.insert_image(
                    index,
                    12,
                    str(png),
                    {
                        "x_scale": 0.36,
                        "y_scale": 0.36,
                        "x_offset": 2,
                        "y_offset": 2,
                        "object_position": 1,
                    },
                )
                image_counts["hmm"] += 1
            else:
                sheet.write(index, 12, "No reportable Tas-core HMM domain")
                image_counts["missing_hmm"] += 1
            height = min(500 / max(1, len(logos)), 112)
            for li, (ai, png) in enumerate(logos):
                with Image.open(png) as im:
                    scale = min(550 / im.width, height / im.height)
                sheet.insert_image(
                    index,
                    14,
                    str(png),
                    {
                        "x_scale": scale,
                        "y_scale": scale,
                        "x_offset": 2,
                        "y_offset": round(li * height) + 2,
                        "object_position": 1,
                    },
                )
            sheet.set_row(
                index,
                min(409, max(94, (text.count("\n") + 1) * 12.8, len(logos) * height * 0.75 + 5)),
            )
        if index % 500 == 0:
            print(f"Excel {kind}: {index}/{len(selected)} anchors", flush=True)
    sheet.autofilter(0, 0, len(selected), len(fields) - 1)
    if kind != "DR":
        sheet.freeze_panes(1, 0)
    return image_counts
