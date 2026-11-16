"""Column layout and cell formats of the Excel reports."""

WITH_FIELDS = (
    "protein_name",
    "sequence_name",
    "strand",
    "protein_sequence",
    "tier",
    "score",
    "full_E",
    "domain_E",
    "HMM_cov",
    "protein_cov",
    "ntd_annotation",
    "tigr_array",
    "hmm_alignment",
)
WITHOUT_FIELDS = WITH_FIELDS[:11]
NUMERIC = {"score", "full_E", "domain_E", "HMM_cov", "protein_cov"}
WIDTHS = {
    "protein_name": 22,
    "sequence_name": 22,
    "strand": 8,
    "protein_sequence": 42,
    "tier": 24,
    "score": 12,
    "full_E": 14,
    "domain_E": 14,
    "HMM_cov": 12,
    "protein_cov": 12,
    "ntd_annotation": 22,
    "hmm_alignment": 76,
}


def _formats(book):
    return {
        "header": book.add_format(
            {"bold": True, "align": "center", "valign": "vcenter", "text_wrap": True}
        ),
        "body": book.add_format({"valign": "top", "text_wrap": True}),
        "mono": book.add_format(
            {"valign": "top", "text_wrap": True, "font_name": "Consolas", "font_size": 10}
        ),
        "score": book.add_format({"valign": "top", "num_format": "0.000"}),
        "scientific": book.add_format({"valign": "top", "num_format": "0.0E+00"}),
        "coverage": book.add_format({"valign": "top", "num_format": "0.000"}),
    }


def _configure(sheet, fields, fmt, array_width):
    sheet.set_row(0, 26)
    for index, field in enumerate(fields):
        style = fmt["body"]
        if field == "score":
            style = fmt["score"]
        elif field in ("full_E", "domain_E"):
            style = fmt["scientific"]
        elif field in ("HMM_cov", "protein_cov"):
            style = fmt["coverage"]
        elif field == "tigr_array":
            style = fmt["mono"]
        sheet.set_column(
            index, index, array_width if field == "tigr_array" else WIDTHS.get(field, 18), style
        )
        sheet.write(0, index, field, fmt["header"])


def _write_numeric(sheet, row, col, value):
    import math

    if value is None or str(value).strip() == "":
        sheet.write_blank(row, col, None)
        return
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite numerical annotation at row {row + 1}")
    sheet.write_number(row, col, number)
