"""Scientific source/ingestion regression tests; never run native FEFF.

Run with ``python -m unittest -v test_round1`` from this source directory.
Mutation tests use disposable copies. Database inheritance is checked read-only
when round1-candidates.db exists and FEFFDB_ROUND1_BASE_DB names its pinned base.
Missing calculation/database fixtures are reported as skips, not as validation.
"""
from __future__ import annotations

import copy
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

import pilot


ABSORBERS = {"Na", "Mg", "Al", "Si", "P", "S", "Cl", "K", "Sc", "Br"}
EXPECTED_CN = {element: (4 if element in {"Si", "P", "S"} else 6)
               for element in ABSORBERS}


def reference(element="Mg"):
    return copy.deepcopy(next(r for r in pilot.records() if r["absorber"] == element))


def card(text, name):
    values = [line.split("*", 1)[0].split()[1:] for line in text.splitlines()
              if re.match(rf"^{re.escape(name)}\s", line)]
    if len(values) != 1:
        raise AssertionError(f"Expected one active {name} card; found {values}")
    return values[0]


class NoNativeFEFF(unittest.TestCase):
    def setUp(self):
        guard = patch.object(pilot, "run_reference", side_effect=AssertionError(
            "These tests must not start native FEFF"))
        guard.start()
        self.addCleanup(guard.stop)


class PreparedSourceTests(NoNativeFEFF):
    def test_complete_requested_absorber_inventory(self):
        records = pilot.records()
        elements = [r["absorber"] for r in records]
        self.assertEqual(set(elements), ABSORBERS)
        self.assertEqual(len(records), 10)
        self.assertEqual(len({r["id"] for r in records}), 10)
        self.assertEqual(len(elements), len(set(elements)))
        self.assertEqual(len({r["cif_sha256"] for r in records}), 9)
        self.assertEqual(sum(len(r["absorber_site_indices"]) for r in records), 11)

    def test_every_real_cif_site_has_consistent_input_and_population(self):
        for record in pilot.records():
            weights, multiplicities = [], []
            for index in record["absorber_site_indices"]:
                with self.subTest(reference=record["id"], site=index):
                    text, metadata = pilot.prepare(record, index)
                    repeat, repeated_metadata = pilot.prepare(record, index)
                    self.assertEqual(text, repeat)
                    self.assertEqual(metadata, repeated_metadata)
                    recipe = metadata["recipe"]
                    self.assertEqual(card(text, "EDGE"), [record["edge"]])
                    self.assertEqual(int(card(text, "NLEG")[0]), recipe["max_legs"])
                    self.assertEqual(float(card(text, "EXAFS")[0]), recipe["kmax_invA"])
                    self.assertEqual(float(card(text, "S02")[0]), recipe["s02"])
                    self.assertEqual(int(card(text, "EXCHANGE")[0]), recipe["exchange"])
                    self.assertEqual(float(card(text, "RPATH")[0]), recipe["path_radius_A"])
                    scf = card(text, "SCF")
                    self.assertEqual(float(scf[0]), recipe["scf_radius_A"])
                    self.assertEqual(int(scf[2]), recipe["scf_iterations"])
                    self.assertEqual(float(scf[3]), recipe["scf_mixing"])
                    self.assertEqual(card(text, "CRITERIA"), ["0", "0"])
                    atoms = [row for _, row in pilot.card_rows(text, "ATOMS")]
                    self.assertEqual(len(atoms), metadata["atom_count"])
                    extent = max(math.dist([0, 0, 0], list(map(float, a[:3]))) for a in atoms)
                    self.assertEqual(extent, metadata["actual_cluster_radius_A"])
                    self.assertGreaterEqual(extent, max(recipe["path_radius_A"], recipe["scf_radius_A"]))
                    self.assertEqual(len(metadata["first_shell_distances_A"]), EXPECTED_CN[record["absorber"]])
                    self.assertTrue(all(0 < r < record["shell_cutoff_A"] for r in metadata["first_shell_distances_A"]))
                    weights.append(metadata["site_weight"])
                    multiplicities.append(metadata["site_multiplicity"])
            self.assertAlmostEqual(sum(weights), 1.0)
            self.assertEqual(weights, [n / sum(multiplicities) for n in multiplicities])

    def test_supported_radius_and_scf_changes_are_written_to_input(self):
        recipe = dict(pilot.RECIPE, cluster_radius_A=8.5, scf_radius_A=4.5,
                      scf_iterations=40, scf_mixing=0.15, path_radius_A=3.5)
        text, metadata = pilot.prepare(reference("Na"), 1, recipe=recipe)
        self.assertEqual(metadata["recipe"], recipe)
        self.assertEqual(card(text, "SCF"), ["4.5", "0", "40", "0.15", "1"])
        self.assertEqual(card(text, "RPATH"), ["3.5"])

    def test_unsupported_recipe_values_cannot_be_reported_as_executed(self):
        unsupported = {"max_legs": 4, "kmax_invA": 12, "s02": 0.75,
                       "exchange": 1, "polarization": [1, 0, 0],
                       "with_h": True, "rng_seed": 987}
        for name, value in unsupported.items():
            with self.subTest(card=name), self.assertRaisesRegex(ValueError, "Unsupported recipe"):
                pilot.prepare(reference(), 1, recipe=dict(pilot.RECIPE, **{name: value}))
        for recipe in ({k: v for k, v in pilot.RECIPE.items() if k != "s02"},
                       dict(pilot.RECIPE, undocumented_option=True)):
            with self.subTest(recipe=recipe), self.assertRaisesRegex(ValueError, "Unsupported recipe"):
                pilot.prepare(reference(), 1, recipe=recipe)


class SourceIntegrityTests(NoNativeFEFF):
    def setUp(self):
        super().setUp()
        temporary = tempfile.TemporaryDirectory(prefix="round1-cif-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.record = reference()
        self.path = self.root / "references" / self.record["id"] / "source.cif"
        self.path.parent.mkdir(parents=True)
        shutil.copyfile(pilot.ROOT / "references" / self.record["id"] / "source.cif", self.path)

    def test_changed_cif_bytes_and_wrong_pin_are_refused(self):
        original = self.path.read_bytes()
        for changed in (b"not a CIF\n", original + b"\n# changed after review\n"):
            self.path.write_bytes(changed)
            with self.subTest(changed=changed[:20]), self.assertRaisesRegex(ValueError, "CIF hash mismatch"):
                pilot.load_structure(self.record, self.root)
        self.path.write_bytes(original)
        self.record["cif_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "CIF hash mismatch"):
            pilot.load_structure(self.record, self.root)

    def test_partial_occupancy_is_rejected_even_with_an_updated_hash(self):
        original = self.path.read_text()
        changed = original.replace("O1 O2- 4 b 0.5 0.5 0.5 1. 0 d", "O1 O2- 4 b 0.5 0.5 0.5 0.5 0 d")
        self.assertNotEqual(original, changed)
        self.path.write_text(changed)
        self.record["cif_sha256"] = pilot.sha(self.path)
        with self.assertRaisesRegex(ValueError, "Disordered/partial occupancy"):
            pilot.load_structure(self.record, self.root)

    def test_changed_inequivalent_absorber_inventory_is_rejected(self):
        self.record["absorber_site_indices"] = [2]
        with self.assertRaisesRegex(ValueError, "Absorber site inventory changed"):
            pilot.load_structure(self.record, self.root)


class CompleteShellTests(NoNativeFEFF):
    def setUp(self):
        super().setUp()
        self.distances = [2.00004] * 2 + [2.00496] * 4
        self.paths = [dict(reff=2.0000, degen=2), dict(reff=2.0050, degen=4)]

    def test_close_split_shells_keep_multiplicity_and_tolerate_native_rounding(self):
        self.assertAlmostEqual(pilot.match_coordination(self.distances, self.paths[::-1]), 0.00004)

    def test_missing_duplicate_and_averaged_paths_cannot_replace_the_full_shell(self):
        cases = [(self.paths[:1], "Incomplete coordination"),
                 (self.paths + self.paths[:1], "Incomplete coordination"),
                 ([dict(reff=sum(self.distances) / 6, degen=6)], "distance mismatch")]
        for paths, message in cases:
            with self.subTest(paths=paths), self.assertRaisesRegex(ValueError, message):
                pilot.match_coordination(self.distances, paths)

    def test_nonphysical_degeneracy_or_distance_is_rejected(self):
        for value in (0, -1, 1.5, float("nan"), float("inf")):
            with self.subTest(degen=value), self.assertRaisesRegex(ValueError, "degeneracy"):
                pilot.match_coordination([2.0], [dict(reff=2.0, degen=value)])
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(reff=value), self.assertRaises(ValueError):
                pilot.match_coordination([2.0], [dict(reff=value, degen=1)])
            with self.subTest(cif_distance=value), self.assertRaises(ValueError):
                pilot.match_coordination([value], [dict(reff=2.0, degen=1)])


class ArchivedRunTests(NoNativeFEFF):
    def setUp(self):
        super().setUp()
        self.record = reference()
        index = self.record["absorber_site_indices"][0]
        source = pilot.ROOT / "calculations" / self.record["id"] / f"site-{index}"
        if not (source / "validation.json").is_file():
            self.skipTest("Completed MgO archive not available")
        temporary = tempfile.TemporaryDirectory(prefix="round1-run-test-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name) / source.name
        shutil.copytree(source, self.directory)

    def test_every_archived_site_has_full_shell_and_usable_larch_chi(self):
        pending = []
        selected_paths = 0
        for record in pilot.records():
            for index in record["absorber_site_indices"]:
                directory = pilot.ROOT / "calculations" / record["id"] / f"site-{index}"
                if not (directory / "validation.json").is_file():
                    pending.append(f"{record['id']}/site-{index}")
                    continue
                with self.subTest(reference=record["id"], site=index):
                    report = pilot.validate_site(record, directory)
                    self.assertEqual(report["coordination_sum"], EXPECTED_CN[record["absorber"]])
                    self.assertEqual(report["cif_coordination"], EXPECTED_CN[record["absorber"]])
                    self.assertEqual(report["larch_chi_evaluation"], "passed")
                    self.assertEqual(report["review_status"], "candidate")
                    self.assertFalse(report["experimental_validation"])
                    self.assertFalse(report["scf_hit_iteration_limit"])
                    selected_paths += len(report["selected_paths"])
        if pending:
            self.skipTest(f"Validated completed sites; pending: {', '.join(pending)}")
        self.assertEqual(selected_paths, 16)

    def test_corrupted_structural_metadata_is_rejected_before_rehashing(self):
        path = self.directory / "run.json"
        original = json.loads(path.read_text())
        corruptions = {
            "site_weight": original["site_weight"] / 2,
            "site_multiplicity": original["site_multiplicity"] + 1,
            "first_shell_distances_A": [r + 0.001 for r in original["first_shell_distances_A"]],
            "actual_cluster_radius_A": original["actual_cluster_radius_A"] + 0.5,
        }
        for field, value in corruptions.items():
            path.write_text(json.dumps(dict(original, **{field: value})))
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "structural metadata"):
                pilot.validate_site(self.record, self.directory)

    def test_failed_or_missing_native_module_cannot_be_a_candidate(self):
        path = self.directory / "run.json"
        original = json.loads(path.read_text())
        for failure in ("failed", "missing", "nonzero"):
            changed = copy.deepcopy(original)
            if failure == "failed":
                changed["execution_status"] = "failed"
            elif failure == "missing":
                changed["module_results"].pop()
            else:
                changed["module_results"][-1]["returncode"] = 1
            path.write_text(json.dumps(changed))
            with self.subTest(failure=failure), self.assertRaisesRegex(ValueError, "six modules"):
                pilot.validate_site(self.record, self.directory)

    def test_changed_input_cannot_reuse_archived_validation(self):
        path = self.directory / "feff.inp"
        text = path.read_text()
        changed = text.replace("S02 1.0", "S02 0.5")
        self.assertNotEqual(text, changed)
        path.write_text(changed)
        with self.assertRaisesRegex(ValueError, "Input has changed"):
            pilot.validate_site(self.record, self.directory)

    def test_removed_first_shell_path_cannot_be_silently_pruned(self):
        report = json.loads((self.directory / "validation.json").read_text())
        (self.directory / report["selected_paths"][0]["filename"]).unlink()
        with self.assertRaisesRegex(ValueError, "Incomplete coordination|No FEFF path files"):
            pilot.validate_site(self.record, self.directory)

    def test_native_fatal_warning_and_iteration_limit_are_not_validated(self):
        path = self.directory / "feff.log"
        original = path.read_text()
        for suffix, message in [("Results might be meaningless", "Fatal FEFF diagnostic"),
                                ("SCF ITERATION NUMBER 30 OUT OF 30", "SCF reached the iteration limit")]:
            path.write_text(original + "\n" + suffix + "\n")
            with self.subTest(diagnostic=suffix), self.assertRaisesRegex(ValueError, message):
                pilot.validate_site(self.record, self.directory)

    def test_modified_raw_path_bytes_cannot_silently_get_a_new_checksum(self):
        report = json.loads((self.directory / "validation.json").read_text())
        path = self.directory / report["selected_paths"][0]["filename"]
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "Archived calculation bytes changed"):
            pilot.validate_site(self.record, self.directory)


class DatabaseInheritanceTests(NoNativeFEFF):
    def test_candidate_preserves_pinned_base_rows_and_blobs_read_only(self):
        candidate = pilot.ROOT / "round1-candidates.db"
        supplied = os.environ.get("FEFFDB_ROUND1_BASE_DB")
        if not candidate.is_file() or not supplied:
            self.skipTest("Need candidate database and FEFFDB_ROUND1_BASE_DB for inheritance verification")
        base = Path(supplied).resolve()
        self.assertTrue(base.is_file())
        config = json.loads((pilot.ROOT / "references.json").read_text())
        manifest = json.loads(candidate.with_suffix(".manifest.json").read_text())
        before = {p: pilot.sha(p) for p in (base, candidate)}
        self.assertEqual(before[base], config["upstream_base"]["sha256"])
        self.assertEqual(before[base], manifest["base_sha256"])
        self.assertEqual(before[candidate], manifest["candidate_sha256"])
        for database in (base, candidate):
            self.assertFalse(Path(str(database) + "-wal").exists())
            self.assertFalse(Path(str(database) + "-shm").exists())
        with sqlite3.connect(base.as_uri() + "?mode=ro&immutable=1", uri=True) as old, \
                sqlite3.connect(candidate.resolve().as_uri() + "?mode=ro&immutable=1", uri=True) as new:
            self.assertEqual(new.execute("PRAGMA integrity_check").fetchall(), [("ok",)])
            self.assertEqual(new.execute("PRAGMA foreign_key_check").fetchall(), [])
            for (table,) in old.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"):
                self.assertRegex(table, r"^[a-z_]+$")
                old_rows = old.execute(f'SELECT * FROM "{table}"').fetchall()
                new_rows = set(new.execute(f'SELECT * FROM "{table}"').fetchall())
                self.assertTrue(all(row in new_rows for row in old_rows), table)
            self.assertEqual(old.execute("SELECT count(*) FROM feffdat").fetchone()[0], manifest["upstream_paths"])
            self.assertEqual(new.execute("SELECT count(*) FROM feffdat").fetchone()[0], manifest["total_paths"])
            self.assertEqual(manifest["upstream_paths"] + manifest["added_paths"], manifest["total_paths"])
            new_ids = {row["id"] for row in manifest["paths"]}
            self.assertEqual(len(new_ids), manifest["added_paths"])
            for ident in new_ids:
                self.assertIsNone(old.execute("SELECT id FROM feffdat WHERE id=?", (ident,)).fetchone())
                self.assertEqual(new.execute("SELECT status FROM drxas_path_review WHERE feffdat_id=?", (ident,)).fetchone(), ("candidate",))
        self.assertEqual({p: pilot.sha(p) for p in before}, before)


if __name__ == "__main__":
    unittest.main()
