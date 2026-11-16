"""FASTA mode and arbitration tests."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

from tigrdisco.dr_tigr_engine import arbitration


class FastaModeTests(unittest.TestCase):
    def _run(self, *extra):
        from tigrdisco import run_fasta

        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            run_fasta.main(
                [
                    str(root / "examples" / "example_windows.fa"),
                    "--out",
                    str(out),
                    "--workers",
                    "2",
                    "--limit",
                    "8",
                    *extra,
                ]
            )
            rows = list(csv.DictReader((out / "loci.csv").open()))
            manifest = json.loads((out / "manifest.json").read_text())
            jsonl = [json.loads(line) for line in (out / "loci.jsonl").open()]
        return rows, manifest, jsonl

    def test_fasta_mode_finds_anchors_and_runs_both_engines(self):
        rows, manifest, jsonl = self._run()
        self.assertEqual(manifest["input_records"], 8)
        self.assertGreater(manifest["orfs_called"], 0)
        self.assertEqual(len(rows), manifest["N"])
        self.assertEqual(len(jsonl), manifest["N"])
        self.assertGreater(manifest["N"], 0, "no Tas anchor found in the bundled example")
        self.assertIn("pyrodigal", manifest["anchor_mode"])
        for row in rows:
            self.assertIn(row["architecture"], {"SL", "DR", "SL+DR", "neither"})
            self.assertIn(row["tier"], {"high_confidence_tas", "possible_or_remote_tas"})
            self.assertIn(row["anchor_strand"], {"+", "-"})
        for record in jsonl:
            self.assertIn("anchor", record)
            self.assertIn("contig", record)
            self.assertLessEqual(record["CENTRAL_start0"], record["CENTRAL_end0"])

    def test_high_tier_is_a_subset_of_broad(self):
        broad, _, _ = self._run()
        high, _, _ = self._run("--tier", "high")
        self.assertLessEqual(len(high), len(broad))
        self.assertTrue({r["locus"] for r in high} <= {r["locus"] for r in broad})

    def test_output_directory_must_be_empty(self):
        from tigrdisco import run_fasta

        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            out.mkdir()
            (out / "existing.txt").write_text("x")
            with self.assertRaises(SystemExit):
                run_fasta.main(
                    [
                        str(root / "examples" / "example_windows.fa"),
                        "--out",
                        str(out),
                        "--limit",
                        "2",
                    ]
                )

    def test_fasta_validation_rejects_duplicates_and_bad_bases(self):
        from tigrdisco import run_fasta

        with tempfile.TemporaryDirectory() as tmp:
            dup = Path(tmp) / "dup.fa"
            dup.write_text(">a\nACGT\n>a\nACGT\n")
            with self.assertRaises(ValueError):
                list(run_fasta.read_fasta(dup))
            bad = Path(tmp) / "bad.fa"
            bad.write_text(">a\nACGTZZ\n")
            with self.assertRaises(ValueError):
                list(run_fasta.read_fasta(bad))


class InRunArbitrationTests(unittest.TestCase):
    DR = [
        {
            "units": [
                {
                    "contig_index": 100,
                    "R1": "AAA",
                    "spacer_a": "CCCCCCCCC",
                    "R2": "GGG",
                    "spacer_b": "TTTTTTTTT",
                },
                {
                    "contig_index": 136,
                    "R1": "AAA",
                    "spacer_a": "CCCCCCCCC",
                    "R2": "GGG",
                    "spacer_b": "TTTTTTTTT",
                },
            ]
        }
    ]

    @staticmethod
    def _sl(start, end, stem):
        return [
            {
                "units": [
                    {
                        "start": start,
                        "end": end,
                        "extended_start": start - 6,
                        "extended_end": end + 10,
                        "stem_length": stem,
                    }
                ]
            }
        ]

    def test_dr_intervals_span_the_four_segments(self):
        units, spans = arbitration.dr_intervals(self.DR)
        self.assertEqual(units, [(100, 124), (136, 160)])
        self.assertEqual(spans, [(100, 160)])

    def test_off_returns_the_stem_loop_result_untouched(self):
        sl = self._sl(101, 130, 5)
        self.assertIs(arbitration.arbitrate_locus(sl, self.DR, 0, "off"), sl)
        self.assertIs(arbitration.arbitrate_locus(sl, [], 0, "distance"), sl)

    def test_unit_mode_drops_an_overlapping_array_whatever_its_stem(self):
        self.assertEqual(
            arbitration.arbitrate_locus(self._sl(101, 130, 12), self.DR, 0, "unit"), []
        )
        self.assertEqual(
            len(arbitration.arbitrate_locus(self._sl(300, 330, 3), self.DR, 0, "unit")), 1
        )

    def test_intersect_spares_a_stable_hairpin(self):
        self.assertEqual(
            arbitration.arbitrate_locus(self._sl(101, 130, 5), self.DR, 0, "intersect"), []
        )
        self.assertEqual(
            len(arbitration.arbitrate_locus(self._sl(101, 130, 8), self.DR, 0, "intersect")), 1
        )
        self.assertEqual(
            arbitration.arbitrate_locus(self._sl(101, 130, 8), self.DR, 0, "intersect", stem=0), []
        )

    def test_distance_reaches_further_than_intersect(self):
        far = self._sl(900, 930, 5)
        self.assertEqual(len(arbitration.arbitrate_locus(far, self.DR, 0, "intersect")), 1)
        self.assertEqual(arbitration.arbitrate_locus(far, self.DR, 0, "distance", gap=1200), [])
        self.assertEqual(len(arbitration.arbitrate_locus(far, self.DR, 0, "distance", gap=100)), 1)

    def test_window_offset_is_applied_to_the_stem_loop_side(self):
        sl = self._sl(1, 30, 5)
        self.assertEqual(len(arbitration.arbitrate_locus(sl, self.DR, 0, "intersect")), 1)
        self.assertEqual(arbitration.arbitrate_locus(sl, self.DR, 100, "intersect"), [])

    def test_covered_only_settles_an_array_that_is_entirely_on_the_dual_repeat(self):
        def arrays(*iv):
            return [
                {
                    "units": [
                        {
                            "start": a,
                            "end": b,
                            "extended_start": a - 5,
                            "extended_end": b + 5,
                            "stem_length": 5,
                        }
                        for a, b in iv
                    ]
                }
            ]

        partly = arrays((101, 120), (500, 520))
        wholly = arrays((101, 120), (137, 156))
        self.assertEqual(arbitration.arbitrate_locus(partly, self.DR, 0, "unit"), [])
        self.assertEqual(arbitration.arbitrate_locus(wholly, self.DR, 0, "unit"), [])
        self.assertEqual(len(arbitration.arbitrate_locus(partly, self.DR, 0, "covered")), 1)
        self.assertEqual(arbitration.arbitrate_locus(wholly, self.DR, 0, "covered"), [])
        self.assertEqual(
            arbitration.arbitrate_locus(partly, self.DR, 0, "covered", coverage=0.5), []
        )

    def test_coverage_must_be_a_fraction(self):
        from tigrdisco import run_fasta

        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "x.fa"
            fa.write_text(">a\nACGT\n")
            for bad in ("0", "1.5", "-0.2"):
                with self.assertRaises(SystemExit):
                    run_fasta.main(
                        [str(fa), "--out", str(Path(tmp) / "o"), "--arbitrate-coverage", bad]
                    )

    def test_unknown_mode_is_rejected(self):
        with self.assertRaises(ValueError):
            arbitration.arbitrate_locus(self._sl(1, 30, 5), self.DR, 0, "sometimes")

    def test_default_is_covered_in_both_pipeline_entry_points(self):
        from tigrdisco import run_fasta, run_database

        for parser in (run_fasta.build_parser(), run_database.build_parser()):
            args = parser.parse_args(["x", "--out", "y"])
            self.assertEqual(args.arbitrate, "covered")
            self.assertEqual(args.arbitrate_coverage, 1.0)
            self.assertEqual(args.arbitrate_gap, 1200)
            self.assertEqual(args.arbitrate_stem, 8)

    def test_negative_arbitration_parameters_are_rejected(self):
        from tigrdisco import run_fasta

        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "x.fa"
            fa.write_text(">a\nACGT\n")
            with self.assertRaises(SystemExit):
                run_fasta.main([str(fa), "--out", str(Path(tmp) / "o"), "--arbitrate-gap", "-1"])

    def test_default_run_records_and_applies_covered(self):
        from tigrdisco import run_fasta

        root = Path(__file__).resolve().parents[1]
        res = {}
        for tag, extra in (("default", []), ("off", ["--arbitrate", "off"])):
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "out"
                run_fasta.main(
                    [
                        str(root / "examples" / "example_windows.fa"),
                        "--out",
                        str(out),
                        "--workers",
                        "2",
                        "--limit",
                        "8",
                        *extra,
                    ]
                )
                res[tag] = (
                    list(csv.DictReader((out / "loci.csv").open())),
                    json.loads((out / "manifest.json").read_text()),
                )
        self.assertEqual(
            res["default"][1]["arbitration"],
            {"mode": "covered", "gap": 1200, "stem": 8, "coverage": 1.0},
        )
        self.assertEqual(res["off"][1]["arbitration"]["mode"], "off")
        d, o = res["default"][0], res["off"][0]
        self.assertEqual(len(d), len(o))
        for a, b in zip(d, o):
            self.assertEqual(a["locus"], b["locus"])
            self.assertEqual(a["DR"], b["DR"])
            self.assertLessEqual(int(a["SL"]), int(b["SL"]))


if __name__ == "__main__":
    unittest.main()
