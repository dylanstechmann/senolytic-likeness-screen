from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

try:
    import rdkit  # noqa: F401

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

import screen  # noqa: E402


def make_row(smiles: str):
    descriptors = screen.compute_descriptors(smiles)
    if descriptors is None:
        raise AssertionError(f"failed to parse {smiles}")
    return {"name": smiles, "descriptors": descriptors}


class DescriptorTests(unittest.TestCase):
    @unittest.skipUnless(RDKIT_AVAILABLE, "RDKit not installed")
    def test_ethanol_and_benzene_descriptors(self) -> None:
        ethanol = screen.compute_descriptors("CCO")
        benzene = screen.compute_descriptors("c1ccccc1")
        self.assertAlmostEqual(ethanol["mw"], 46.07, delta=0.05)
        self.assertEqual(ethanol["hbd"], 1)
        self.assertEqual(benzene["rings"], 1)
        self.assertAlmostEqual(benzene["arom_prop"], 1.0, delta=0.01)
        self.assertEqual(benzene["fsp3"], 0.0)

    @unittest.skipUnless(RDKIT_AVAILABLE, "RDKit not installed")
    def test_salt_strips_to_largest_fragment(self) -> None:
        plain = screen.compute_descriptors("CCO")
        salted = screen.compute_descriptors("CCO.Cl")
        self.assertEqual(plain["mw"], salted["mw"])

    def test_invalid_smiles_returns_none(self) -> None:
        if not RDKIT_AVAILABLE:
            self.skipTest("RDKit not installed")
        self.assertIsNone(screen.compute_descriptors("not-a-smiles"))


class ScoringTests(unittest.TestCase):
    def test_score_is_monotone_and_bounded(self) -> None:
        self.assertEqual(screen.proximity_score(0.0), 1.0)
        self.assertLess(screen.proximity_score(3.0), screen.proximity_score(1.0))
        self.assertGreater(screen.proximity_score(10.0), 0.0)

    def test_identical_vectors_have_zero_distance(self) -> None:
        a = {key: 1.5 for key in screen.DESCRIPTOR_KEYS}
        self.assertEqual(screen.z_distance(a, a), 0.0)

    def test_nearest_reference_is_selected(self) -> None:
        near = {key: 0.0 for key in screen.DESCRIPTOR_KEYS}
        far = {key: 9.0 for key in screen.DESCRIPTOR_KEYS}
        candidate = {key: 0.5 for key in screen.DESCRIPTOR_KEYS}
        candidate_row = {"name": "cand", "z": candidate}
        references = [
            {"name": "far", "z": far},
            {"name": "near", "z": near},
        ]
        screen.score_candidates([candidate_row], references)
        self.assertEqual(candidate_row["nearest_reference"], "near")
        self.assertEqual(candidate_row["score"], screen.proximity_score(candidate_row["distance"]))

    def test_standardize_adds_z_columns(self) -> None:
        rows = []
        for value in (1.0, 2.0, 3.0):
            row = {"descriptors": {key: value for key in screen.DESCRIPTOR_KEYS}}
            rows.append(row)
        screen.standardize(rows)
        self.assertIn("z", rows[0])
        # population std over (1, 2, 3) is sqrt(2/3); z-scores are +/- sqrt(3/2)
        self.assertAlmostEqual(rows[0]["z"]["mw"], -(3 / 2) ** 0.5, delta=0.001)
        self.assertAlmostEqual(rows[2]["z"]["mw"], (3 / 2) ** 0.5, delta=0.001)


class FileTests(unittest.TestCase):
    def test_read_tsv_skips_comments_and_header(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ref.tsv"
            path.write_text("name\tquery\n# comment\nNAVITOCLAX\tnavitoclax\n")
            rows = screen.read_tsv(path)
        self.assertEqual(rows, [("NAVITOCLAX", "navitoclax")])

    def test_finish_creates_provenance_for_new_output_directory(self) -> None:
        reference = {
            "chembl_id": "CHEMBL1", "name": "reference", "smiles": "CCO",
            "descriptors": {key: 1.0 for key in screen.DESCRIPTOR_KEYS},
        }
        candidate = {
            "chembl_id": "CHEMBL2", "name": "candidate", "smiles": "CCC",
            "descriptors": {key: 2.0 for key in screen.DESCRIPTOR_KEYS},
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "fresh-output"
            result = screen.finish(
                root, {"started": "2026-09-26"}, [reference], [candidate], 1
            )
            self.assertEqual(result, 0)
            self.assertTrue((root / "data" / "results.tsv").is_file())
            self.assertEqual(len(list((root / "data" / "provenance").glob("*_run.json"))), 1)
            self.assertTrue((root / "REPORT.md").is_file())


if __name__ == "__main__":
    unittest.main()
