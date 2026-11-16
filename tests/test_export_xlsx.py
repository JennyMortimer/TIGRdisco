"""Excel export tests."""

import json
import tempfile
import unittest
from pathlib import Path
from openpyxl import load_workbook
from tigrdisco.export_xlsx import export_run, tables


class ExcelExportTest(unittest.TestCase):
    def record(self):
        return {
            "locus": "=unsafe",
            "record_length": 1000,
            "reverse_complemented": True,
            "CENTRAL_start0": 200,
            "CENTRAL_end0": 500,
            "SL_window_offset0": 50,
            "protein_sequence": "MAAA",
            "SL": [],
            "DR": [
                {
                    "units": [
                        {
                            "contig_index": 100,
                            "R1": "AC",
                            "spacer_a": "GG",
                            "R2": "TT",
                            "spacer_b": "CA",
                        }
                    ]
                }
            ],
        }

    def test_reverse_coordinates(self):
        a, _, u = tables([self.record()], "DR")
        self.assertEqual((a[0]["anchor_input_start1"], a[0]["anchor_input_end1"]), (501, 800))
        self.assertEqual((u[0]["input_start1"], u[0]["input_end1"]), (893, 900))
        self.assertEqual(u[0]["unit_sequence_oriented"], "ACGGTTCA")

    def test_empty_and_both_and_literal_cells(self):
        for both in (False, True):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                r = self.record()
                records = []
                if both:
                    r["SL"] = [
                        {
                            "strand": "+",
                            "units": [
                                {
                                    "start": 1,
                                    "end": 7,
                                    **{
                                        k: "A"
                                        for k in (
                                            "box_C",
                                            "spacer_A",
                                            "box_D",
                                            "connector",
                                            "box_C_prime",
                                            "spacer_B",
                                            "box_D_prime",
                                        )
                                    },
                                }
                            ],
                        }
                    ]
                    records = [r]
                (root / "loci.jsonl").write_text("".join(json.dumps(x) + "\n" for x in records))
                (root / "manifest.json").write_text(
                    json.dumps(
                        {
                            "N": int(both),
                            "SL": int(both),
                            "DR": int(both),
                            "both": int(both),
                            "neither": 0,
                        }
                    )
                )
                export_run(root, root / "xlsx", style="tables")
                for prefix in ("sl_", "dr_"):
                    book = load_workbook(
                        root / "xlsx" / f"{prefix}tas_with_tigr.xlsx", read_only=True
                    )
                    self.assertEqual(book["Anchors"].max_row, 2 if both else 1)
                    if both:
                        self.assertEqual(book["Anchors"]["A2"].value, "=unsafe")
                        self.assertEqual(book["Anchors"]["A2"].data_type, "s")
                    book.close()
                with self.assertRaises(ValueError):
                    export_run(root, root / "xlsx", style="tables")


if __name__ == "__main__":
    unittest.main()
