"""Excel report tests."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile
from openpyxl import load_workbook
from tigrdisco.export_xlsx import export_run
from tigrdisco.template_report import logo_group


class TemplateReportTests(unittest.TestCase):
    def test_all_three_reference_schemas_and_panes(self):
        from xml.etree import ElementTree as ET

        spec = json.loads(
            (Path(__file__).resolve().parent / "data/excel_template_spec.json").read_text()
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r = self.record()
            r["SL"] = [
                {
                    "box_C": "CCA",
                    "box_D": "TG",
                    "support_fraction": 1.0,
                    "array_class": "core",
                    "strand": "+",
                    "units": [
                        {
                            "start": 10 + i * 40,
                            "end": 41 + i * 40,
                            "extended_start": 5 + i * 40,
                            "extended_end": 46 + i * 40,
                            "box_C": "CCA",
                            "spacer_A": "A" * 8,
                            "box_D": "TG",
                            "connector": "AA",
                            "box_C_prime": "CCA",
                            "spacer_B": "T" * 8,
                            "box_D_prime": "TG",
                            "motif_errors": 0,
                        }
                        for i in range(3)
                    ],
                }
            ]
            negative = {**r, "locus": "negative", "SL": [], "DR": []}
            (root / "loci.jsonl").write_text("\n".join(map(json.dumps, [r, negative])))
            (root / "manifest.json").write_text(
                json.dumps({"N": 2, "SL": 1, "DR": 1, "both": 1, "neither": 1})
            )
            with patch(
                "tigrdisco.template_report.ntd.annotate",
                return_value={r["locus"]: "test", "negative": "test"},
            ):
                export_run(root, root / "out")
            ns = {
                "s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
                "x": "http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing",
            }
            for name, expected in spec.items():
                path = root / "out" / name
                book = load_workbook(path, read_only=True)
                self.assertEqual(book.sheetnames[0], expected["sheet"])
                self.assertEqual(
                    list(next(book.worksheets[0].iter_rows(values_only=True))), expected["columns"]
                )
                self.assertEqual(book.worksheets[0].max_column, len(expected["columns"]))
                book.close()
                with zipfile.ZipFile(path) as z:
                    for sheet_path in [
                        n
                        for n in z.namelist()
                        if n.startswith("xl/worksheets/sheet") and n.endswith(".xml")
                    ]:
                        xml = ET.fromstring(z.read(sheet_path))
                        pane = xml.find("s:sheetViews/s:sheetView/s:pane", ns)
                        if pane is not None:
                            self.assertEqual(float(pane.get("xSplit", "0")), 0)
                        if sheet_path == "xl/worksheets/sheet1.xml":
                            self.assertEqual(
                                None if pane is None else pane.attrib, expected["original_pane"]
                            )
                    if name.startswith("sl_"):
                        drawing = ET.fromstring(z.read("xl/drawings/drawing1.xml"))
                        cols = [int(a.find("x:from/x:col", ns).text) for a in drawing]
                        self.assertEqual(cols, [12])

    def record(self):
        return {
            "locus": "example__101_400_-",
            "contig": "example",
            "record_length": 1000,
            "reverse_complemented": True,
            "CENTRAL_start0": 600,
            "CENTRAL_end0": 900,
            "SL_window_offset0": 0,
            "protein_sequence": "M" + "A" * 99,
            "anchor": {
                "hmm_intervals": [[10, 80]],
                "score": 42,
                "full_E": 1e-8,
                "domain_E": 1e-9,
                "tier": "test",
            },
            "SL": [],
            "DR": [
                {
                    "side": "upstream",
                    "bits": 20,
                    "units": [
                        {
                            "contig_index": 100 + 36 * i,
                            "R1": "ACGTACGTA",
                            "spacer_a": "GGGGGGGGG",
                            "R2": "TGCATGCAT",
                            "spacer_b": "AAAAAAAAA",
                            "classification": "Core",
                        }
                        for i in range(4)
                    ],
                }
            ],
        }

    def test_template_columns_images_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r = self.record()
            (root / "loci.jsonl").write_text(json.dumps(r) + "\n")
            (root / "manifest.json").write_text(
                json.dumps({"N": 1, "SL": 0, "DR": 1, "both": 0, "neither": 0})
            )
            with patch(
                "tigrdisco.template_report.ntd.annotate",
                return_value={r["locus"]: "not_evaluated_test"},
            ):
                export_run(root, root / "out")
            book = load_workbook(root / "out/dr_tas_with_tigr.xlsx", read_only=True)
            self.assertEqual(book.sheetnames[0], "Tas with TIGR")
            row = list(book.worksheets[0].iter_rows(values_only=True))
            self.assertEqual(row[0][12:], ("hmm_alignment", "weblogo_array_usage", "weblogo"))
            self.assertEqual(row[1][5], 42)
            self.assertIn("used=4/4", row[1][13])
            book.close()
            with zipfile.ZipFile(root / "out/dr_tas_with_tigr.xlsx") as z:
                self.assertEqual(len([n for n in z.namelist() if n.startswith("xl/media/")]), 2)
                drawing = z.read("xl/drawings/drawing1.xml").decode()
                self.assertIn("<xdr:col>12</xdr:col>", drawing)
                self.assertIn("<xdr:col>14</xdr:col>", drawing)

    def test_missing_protein_fails_instead_of_fake_hmm(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r = self.record()
            r.pop("protein_sequence")
            (root / "loci.jsonl").write_text(json.dumps(r) + "\n")
            (root / "manifest.json").write_text(
                json.dumps({"N": 1, "SL": 0, "DR": 1, "both": 0, "neither": 0})
            )
            with self.assertRaisesRegex(ValueError, "protein sequences"):
                export_run(root, root / "out")
            self.assertFalse((root / "out").exists())

    def test_sl_variable_spacers_do_not_misalign_logo(self):
        units = [
            {
                "box_C": "CCA",
                "box_D": "TG",
                "box_C_prime": "CCA",
                "box_D_prime": "TG",
                "spacer_A": "A" * n,
                "spacer_B": "T" * n,
                "connector": "AC",
            }
            for n in (8, 9, 10)
        ]
        sequences, usage = logo_group({"units": units}, "SL")
        self.assertEqual(sequences, ["CCATGCCATG"] * 3)
        self.assertIn("BOX-ONLY", usage)

    def test_no_eligible_dr_logo_is_explicit(self):
        array = self.record()["DR"][0]
        array["units"] = array["units"][:2]
        sequences, usage = logo_group(array, "DR")
        self.assertFalse(sequences)
        self.assertIn("No logo", usage)

    def test_zero_anchor_template(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "loci.jsonl").write_text("")
            (root / "manifest.json").write_text(
                json.dumps(dict.fromkeys(("N", "SL", "DR", "both", "neither"), 0))
            )
            export_run(root, root / "out")
            for name in ("sl_tas_with_tigr", "dr_tas_with_tigr", "tas_without_tigr"):
                book = load_workbook(root / "out" / f"{name}.xlsx", read_only=True)
                self.assertEqual(book.worksheets[0].max_row, 1)
                book.close()


if __name__ == "__main__":
    unittest.main()
