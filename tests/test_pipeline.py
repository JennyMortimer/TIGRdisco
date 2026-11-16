"""Stem-loop engine tests."""

import inspect
import unittest
from pathlib import Path

from tigrdisco import run_database, run_fasta
from tigrdisco.shared import ntd
from tigrdisco.sl_tigr_engine import detect

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "examples" / "example_windows.fa"
SEGMENTS = ("box_C", "spacer_A", "box_D", "connector", "box_C_prime", "spacer_B", "box_D_prime")


def scan(**overrides):
    config = detect.Config(**overrides)
    arrays = []
    for name, sequence in run_fasta.read_fasta(EXAMPLE):
        arrays.extend(detect.process((name, sequence, config)))
    return arrays


class FixedStepsTest(unittest.TestCase):
    def test_defaults_are_all_on(self):
        cfg = detect.Config()
        self.assertTrue(cfg.system)
        self.assertTrue(cfg.annotate_system)
        self.assertTrue(cfg.annotate_locus)
        self.assertTrue(cfg.qc3)

    def test_no_switch_exposes_them(self):
        for module in (run_fasta, run_database):
            text = inspect.getsource(module)
            for flag in (
                "--no-system",
                "--no-annotate-system",
                "--no-annotate-locus",
                "--qc3",
                "--box",
            ):
                self.assertNotIn(flag, text, f"{flag} must not be a command-line option")


class DetectionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.arrays = scan()
        cls.units = [unit for array in cls.arrays for unit in array["units"]]

    def test_example_produces_arrays(self):
        self.assertGreater(len(self.arrays), 0)
        self.assertGreater(len(self.units), 0)

    def test_field_names(self):
        for name in ("box_C", "box_D", "box_C_length", "box_D_length"):
            self.assertIn(name, self.arrays[0])
        for name in SEGMENTS:
            self.assertIn(name, self.units[0])

    def test_seven_fields_reconstruct_the_unit(self):
        for unit in self.units:
            width = sum(len(unit[k]) for k in SEGMENTS)
            self.assertEqual(width, unit["end"] - unit["start"] + 1)

    def test_frame_length_formula(self):
        for unit in self.units:
            c, d = len(unit["box_C"]), len(unit["box_D"])
            self.assertEqual(len(unit["box_C_prime"]), c)
            self.assertEqual(len(unit["box_D_prime"]), d)
            width = unit["end"] - unit["start"] + 1
            self.assertEqual(
                width,
                2 * (c + d)
                + len(unit["spacer_A"])
                + len(unit["spacer_B"])
                + len(unit["connector"]),
            )

    def test_stems_are_symmetric(self):
        for unit in self.units:
            left = unit["start"] - unit["extended_start"]
            right = unit["extended_end"] - unit["end"]
            self.assertEqual(left, right)
            self.assertEqual(left, unit["stem_length"])

    def test_no_duplicate_arrays(self):
        seen = set()
        for array in self.arrays:
            key = (
                array["contig_id"],
                array["strand"],
                array["units"][0]["start"],
                array["units"][-1]["end"],
                array["box_C"],
                array["box_D"],
            )
            self.assertNotIn(key, seen, f"duplicate array: {key}")
            seen.add(key)

    def test_one_motif_per_locus(self):
        by_locus = {}
        for array in self.arrays:
            by_locus.setdefault(array["contig_id"], set()).add((array["box_C"], array["box_D"]))
        for locus, motifs in by_locus.items():
            self.assertEqual(len(motifs), 1, f"{locus} reports {motifs}")


class SharedTest(unittest.TestCase):
    def test_bundled_profiles_exist(self):
        from tigrdisco.shared import hmm_search

        for path in (hmm_search.DATA, ntd.DEFAULT_TASH_HMM, ntd.DEFAULT_TASR_HMM):
            self.assertTrue(Path(path).is_file(), f"missing bundled profile: {path}")

    def test_short_proteins_are_not_subtyped(self):
        labels = ntd.annotate([("short", "M" * 50)])
        self.assertEqual(labels["short"], "unable_to_classify_<200aa")


class RereadTest(unittest.TestCase):
    def test_reread_is_on_by_default(self):
        self.assertTrue(detect.Config().reread)
        self.assertEqual(detect.REREAD_FRAME, (3, 2))
        self.assertEqual(detect.Config().frame, ())

    def test_frame_pins_both_widths(self):
        self.assertEqual(list(detect._WX(detect.Config())), [2, 3, 4])
        self.assertEqual(list(detect._WY(detect.Config())), [2, 3, 4])
        cfg = detect.Config(frame=(3, 2))
        self.assertEqual(list(detect._WX(cfg)), [3])
        self.assertEqual(list(detect._WY(cfg)), [2])

    def test_reread_keeps_the_locus_set_and_prefers_the_3_2_frame(self):
        def summary(arrays):
            loci = {a["contig_id"] for a in arrays}
            frames = [(a["box_C_length"], a["box_D_length"]) for a in arrays]
            return loci, frames

        on_loci, on_frames = summary(scan())
        off_loci, off_frames = summary(scan(reread=False))
        self.assertEqual(on_loci, off_loci)
        self.assertGreater(on_frames.count((3, 2)), off_frames.count((3, 2)))


if __name__ == "__main__":
    unittest.main()
