"""Scientific source/ingestion regressions; never execute native FEFF.

Run ``python -m unittest -v test_round2`` in this directory. Mutations occur only
in disposable copies. Pending sources, native archives and database inheritance
are explicit skips. FEFFDB_ROUND2_BASE_DB supplies the pinned 374-path base for
the read-only database test; calculations remain candidates, not sample fits.
"""
from __future__ import annotations

import copy
import hashlib
import json
import lzma
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


RARE_EARTHS = {"La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"}
CHEMISTRY = {"Ga", "In", "Hf", "S", "Al", "Si"}
ABSORBERS = RARE_EARTHS | CHEMISTRY
SITE_CN = {"Ga": [4, 6], "In": [6, 6], "Hf": [7], "S": [4],
           "La": [11], "Pr": [9], "Nd": [9], "Sm": [6, 7, 7],
           "Eu": [6, 6], "Gd": [9], "Tb": [7, 7, 6], "Dy": [8],
           "Ho": [6, 6], "Er": [6, 6], "Tm": [8], "Yb": [6, 6], "Lu": [8], "Al": [4], "Si": [4]}


def reference(element="Eu"):
    return copy.deepcopy(next(r for r in pilot.records() if r["absorber"] == element))


def card(text, name):
    values = [line.split("*", 1)[0].split()[1:] for line in text.splitlines()
              if re.match(rf"^{re.escape(name)}\s", line)]
    if len(values) != 1:
        raise AssertionError(f"Expected one active {name} card; found {values}")
    return values[0]


def expected_cn(record, index):
    return record["expected_cn_by_site"][str(index)]


def archive(record, index):
    return pilot.ROOT / "calculations" / record["id"] / f"site-{index}"


class NoNativeFEFF(unittest.TestCase):
    def setUp(self):
        # Avoid patching subprocess globally: import-time platform detection uses it.
        guard = patch.object(pilot, "run_reference", side_effect=AssertionError(
            "These tests must not start native FEFF"))
        guard.start()
        self.addCleanup(guard.stop)


class PreparedSourceTests(NoNativeFEFF):
    def test_complete_requested_absorber_inventory(self):
        records = pilot.records()
        elements = [r["absorber"] for r in records]
        self.assertEqual(len(elements), len(set(elements)))
        self.assertEqual(len({r["id"] for r in records}), len(records))
        self.assertFalse(set(elements) - ABSORBERS)
        missing = ABSORBERS - set(elements)
        if missing and missing <= {"Al", "Si"}:
            self.skipTest(f"Final experimental source selection pending: {sorted(missing)}")
        self.assertEqual(set(elements), ABSORBERS)
        self.assertEqual(len(records), 19)
        self.assertEqual(len(set(elements) & RARE_EARTHS), 13)
        self.assertEqual(len(set(elements) & CHEMISTRY), 6)
        for record in records:
            self.assertEqual(record["edge"], "L3" if record["absorber"] in RARE_EARTHS | {"Hf"} else "K")
            if record["absorber"] in {"Al", "Si"}:
                self.assertEqual(record["scatterer"], {"Al": "N", "Si": "N"}[record["absorber"]])
                self.assertTrue(all(n == 4 for n in record["expected_cn_by_site"].values()))

    def test_every_real_cif_site_has_consistent_input_and_population(self):
        from pymatgen.core import Composition
        from xraydb import atomic_number
        for record in pilot.records():
            weights, multiplicities = [], []
            _, cluster, sites, _ = pilot.load_structure(record)
            self.assertTrue(cluster.struct.is_ordered)
            self.assertEqual([s[0] for s in sites], record["absorber_site_indices"])
            self.assertEqual(set(record["expected_cn_by_site"]), {str(s[0]) for s in sites})
            counts = Composition(cluster.struct.composition).element_composition.reduced_composition.get_el_amt_dict()
            stoichiometry = {atomic_number(e): n for e, n in counts.items()}
            observed_cn = []
            for index in record["absorber_site_indices"]:
                with self.subTest(reference=record["id"], site=index):
                    text, metadata = pilot.prepare(record, index)
                    repeat, repeated_metadata = pilot.prepare(record, index)
                    self.assertEqual((text, metadata), (repeat, repeated_metadata))
                    self.assertEqual(metadata["cif_sha256"], pilot.sha(pilot.ROOT / "references" / record["id"] / "source.cif"))
                    recipe = metadata["recipe"]
                    self.assertEqual(card(text, "EDGE"), [record["edge"]])
                    for name, key, convert in [("NLEG", "max_legs", int), ("EXAFS", "kmax_invA", float),
                                               ("S02", "s02", float), ("EXCHANGE", "exchange", int),
                                               ("RPATH", "path_radius_A", float)]:
                        self.assertEqual(convert(card(text, name)[0]), recipe[key])
                    scf = card(text, "SCF")
                    self.assertEqual([float(scf[0]), int(scf[1]), int(scf[2]), float(scf[3]), int(scf[4])],
                                     [recipe["scf_radius_A"], 0, recipe["scf_iterations"], recipe["scf_mixing"], recipe.get("scf_mixing_iterations", 1)])
                    self.assertEqual(card(text, "CRITERIA"), ["0", "0"])
                    self.assertFalse(re.search(r"^\s*(UNFREEZEF|NOHOLE)\b", text, flags=re.M))
                    atoms = [row for _, row in pilot.card_rows(text, "ATOMS")]
                    self.assertEqual(len(atoms), metadata["atom_count"])
                    extent = max(math.sqrt(sum(float(v)**2 for v in row[:3])) for row in atoms)
                    self.assertEqual(extent, metadata["actual_cluster_radius_A"])
                    self.assertGreaterEqual(extent, max(recipe["path_radius_A"], recipe["scf_radius_A"]))
                    self.assertLessEqual(len(atoms), 400)
                    potentials = [row for _, row in pilot.card_rows(text, "POTENTIALS")]
                    self.assertEqual({int(row[1]) for row in potentials}, set(stoichiometry))
                    for row in potentials[1:]:
                        self.assertAlmostEqual(float(row[5]), stoichiometry[int(row[1])])
                    distances = metadata["first_shell_distances_A"]
                    self.assertEqual(len(distances), expected_cn(record, index))
                    self.assertTrue(all(0 < r < record["shell_cutoff_A"] for r in distances))
                    self.assertEqual(distances, [r for r in metadata["scatterer_distances_to_rpath_A"] if r < record["shell_cutoff_A"]])
                    observed_cn.append(len(distances))
                    weights.append(metadata["site_weight"])
                    multiplicities.append(metadata["site_multiplicity"])
            self.assertAlmostEqual(sum(weights), 1.0)
            self.assertEqual(weights, [n / sum(multiplicities) for n in multiplicities])
            if record["absorber"] in SITE_CN:
                self.assertEqual(observed_cn, SITE_CN[record["absorber"]])

    def test_bixbyite_and_mixed_coordination_weights_follow_actual_sites(self):
        for element, expected in {"Eu": [.25, .75], "Ho": [.25, .75], "Er": [.25, .75],
                                  "Yb": [.75, .25], "Sm": [1/3]*3, "Tb": [1/3]*3}.items():
            r = reference(element)
            with self.subTest(element=element):
                self.assertEqual([pilot.prepare(r, i)[1]["site_weight"] for i in r["absorber_site_indices"]], expected)
        # Equal-weight averaging would incorrectly reduce the distinct Yb populations.
        self.assertNotEqual(pilot.prepare(reference("Yb"), 1)[1]["site_weight"], .5)

    def test_laf3_keeps_two_long_neighbors_and_hfo2_keeps_close_splits(self):
        la = pilot.prepare(reference("La"), 1)[1]["first_shell_distances_A"]
        self.assertEqual(len(la), 11)
        self.assertEqual(sum(r > 2.9 for r in la), 2)
        self.assertTrue(all(abs(r - 2.99925) < .00002 for r in la[-2:]))
        hf = pilot.prepare(reference("Hf"), 1)[1]["first_shell_distances_A"]
        self.assertEqual(len(hf), 7)
        close = [r for r in hf if 2.17 < r < 2.18]
        self.assertEqual(len(close), 2)
        self.assertGreater(abs(close[1] - close[0]), 2 * pilot.DISTANCE_TOLERANCE_A)
        separate = [dict(reff=round(r, 4), degen=1) for r in hf]
        self.assertLessEqual(pilot.match_coordination(hf, separate), pilot.DISTANCE_TOLERANCE_A)
        collapsed = [dict(reff=round(r, 4), degen=1) for r in hf if r not in close]
        collapsed.append(dict(reff=sum(close)/2, degen=2))
        with self.assertRaisesRegex(ValueError, "distance mismatch"):
            pilot.match_coordination(hf, collapsed)

    def test_l3_next_edge_markers_are_physical_boundaries_not_fit_limits(self):
        from xraydb import xray_edge
        from larch.xafs import etok
        boundaries = json.loads((pilot.ROOT / "next-edge-boundaries.json").read_text())
        indexed = {item["absorber"]: item for item in boundaries}
        l3 = [r for r in pilot.records() if r["edge"] == "L3"]
        self.assertEqual(set(indexed), {r["absorber"] for r in l3})
        self.assertEqual(len(indexed), len(boundaries))
        for r in l3:
            with self.subTest(element=r["absorber"]):
                marker = r["next_edge_interference"]
                self.assertEqual(marker, indexed[r["absorber"]])
                lower, upper = [xray_edge(r["absorber"], e).energy for e in ("L3", "L2")]
                self.assertEqual((marker["lower_edge"], marker["next_edge"]), ("L3", "L2"))
                self.assertEqual((marker["lower_energy_eV"], marker["next_energy_eV"]), (lower, upper))
                self.assertEqual(marker["delta_energy_eV"], upper - lower)
                self.assertAlmostEqual(marker["k_at_next_edge_invA"], float(etok(upper - lower)))
                self.assertLess(marker["k_at_next_edge_invA"], pilot.RECIPE["kmax_invA"])
                warnings = " ".join(r["calculation_cautions"])
                self.assertIn("L2", warnings)
                self.assertIn("not a usable-data limit", warnings)
                self.assertEqual(pilot.prepare(r, r["absorber_site_indices"][0])[1]["recipe"]["kmax_invA"], 20)

    def test_supported_scf_changes_and_optional_nmix_are_written_to_input(self):
        recipe = dict(pilot.RECIPE, scf_radius_A=4.5, scf_iterations=40,
                      scf_mixing=0.15, scf_mixing_iterations=3, path_radius_A=3.5)
        text, metadata = pilot.prepare(reference(), 1, recipe=recipe)
        self.assertEqual(metadata["recipe"], recipe)
        self.assertEqual(card(text, "SCF"), ["4.5", "0", "40", "0.15", "3"])
        self.assertEqual(card(text, "RPATH"), ["3.5"])
        self.assertEqual(card(pilot.prepare(reference(), 1, recipe=pilot.RECIPE)[0], "SCF")[-1], "1")

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

    def test_publication_aln_transcription_preserves_measured_geometry_and_identity(self):
        from pymatgen.io.cif import CifParser
        r = reference("Al")
        directory = pilot.ROOT / "references" / r["id"]
        source = json.loads((directory / "provenance.json").read_text())
        transcription = json.loads((directory / "transcription.json").read_text())
        self.assertEqual(pilot.source_identity(r), ("publication", "10.1016/0038-1098(77)90959-0"))
        self.assertNotIn("cod_id", r)
        self.assertFalse(source["original_file_unmodified"])
        self.assertFalse(transcription["original_file_unmodified"])
        self.assertIn("transcribed", source["derivation"])
        self.assertEqual(transcription["source_record_id"], r["doi"])
        self.assertIsNone(source["measurement_temperature_K"])
        self.assertIsNone(source["measurement_pressure_GPa"])
        raw = next(iter(CifParser.from_str((directory / "source.cif").read_text()).as_dict().values()))
        self.assertEqual(raw["_cell_length_a"], "3.110(1)")
        self.assertEqual(raw["_cell_length_c"], "4.980(1)")
        self.assertEqual(raw["_atom_site_fract_z"], ["0", "0.3821(3)"])
        self.assertEqual(raw["_atom_site_U_iso_or_equiv"], ["0.0038(2)", "0.0044(2)"])
        self.assertTrue(all(abs(float(x) + 1/3) < 1.e-14 for x in raw["_atom_site_fract_x"]))
        self.assertTrue(all(abs(float(y) + 2/3) < 1.e-14 for y in raw["_atom_site_fract_y"]))
        text, metadata = pilot.prepare(r, 1)
        self.assertIn("publication", text)
        self.assertIn(r["doi"], text)
        self.assertNotIn("COD", text)
        # Printed Table 3 provides an independent check of the Table 1/2 transcription.
        expected = [1.8891]*3 + [1.9029]
        actual = metadata["first_shell_distances_A"]
        self.assertEqual(len(actual), 4)
        self.assertTrue(all(abs(a-b) < .00005 for a, b in zip(sorted(actual), expected)))
        self.assertGreater(actual[-1] - actual[0], .01)
        self.assertEqual(metadata["cell_composition"], {"Al": 2.0, "N": 2.0})
        for cod_record in [record for record in pilot.records() if record.get("cod_id")]:
            self.assertEqual(pilot.source_identity(cod_record), ("COD", cod_record["cod_id"]))

    def test_invalid_scf_and_nmix_controls_are_refused(self):
        for value in (0, 31, 1.5, True, "3", None):
            with self.subTest(nmix=value), self.assertRaisesRegex(ValueError, "SCF mixing iterations"):
                pilot.prepare(reference(), 1, recipe=dict(pilot.RECIPE, scf_mixing_iterations=value))
        for key, values in {"scf_iterations": (1, 101, 30.0, True), "scf_mixing": (0, 1, -0.1, float("nan"))}.items():
            for value in values:
                with self.subTest(key=key, value=value), self.assertRaisesRegex(ValueError, "Invalid SCF"):
                    pilot.prepare(reference(), 1, recipe=dict(pilot.RECIPE, **{key: value}))


class SourceIntegrityTests(NoNativeFEFF):
    def setUp(self):
        super().setUp()
        temporary = tempfile.TemporaryDirectory(prefix="round2-cif-test-")
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
        changed = original.replace("Eu1 Eu+3 0.25 0.25 0.25 1 0.0", "Eu1 Eu+3 0.25 0.25 0.25 0.5 0.0")
        self.assertNotEqual(original, changed)
        self.path.write_text(changed)
        self.record["cif_sha256"] = pilot.sha(self.path)
        with self.assertRaisesRegex(ValueError, "Disordered/partial occupancy"):
            pilot.load_structure(self.record, self.root)

    def test_changed_inequivalent_absorber_inventory_is_rejected(self):
        self.record["absorber_site_indices"] = [1]
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
        source = archive(self.record, self.record["absorber_site_indices"][0])
        if not (source / "validation.json").is_file():
            self.skipTest("Completed Eu2O3 archive not available")
        temporary = tempfile.TemporaryDirectory(prefix="round2-run-test-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name) / source.name
        shutil.copytree(source, self.directory)

    def test_every_archived_site_has_full_shell_and_native_controls(self):
        pending = []
        completed = 0
        for record in pilot.records():
            for index in record["absorber_site_indices"]:
                directory = archive(record, index)
                if not (directory / "validation.json").is_file():
                    pending.append(f"{record['id']}/site-{index}")
                    continue
                with self.subTest(reference=record["id"], site=index):
                    report = pilot.validate_site(record, directory)
                    self.assertEqual(report["coordination_sum"], expected_cn(record, index))
                    self.assertEqual(report["cif_coordination"], expected_cn(record, index))
                    self.assertEqual(report["larch_chi_evaluation"], "passed")
                    self.assertEqual(report["review_status"], "candidate")
                    self.assertFalse(report["experimental_validation"])
                    self.assertFalse(report["scf_hit_iteration_limit"])
                    controls = json.loads((directory / "pot.json").read_text())
                    recipe = json.loads((directory / "run.json").read_text())["recipe"]
                    proof = report["native_potential_controls"]
                    self.assertEqual(proof, {key: controls[key] for key in proof})
                    self.assertEqual(set(proof), {"iunf", "ixc", "nscmt", "lfms1", "nmix", "ca1", "rfms1", "ihole"})
                    self.assertEqual((proof["iunf"], proof["ixc"], proof["lfms1"]), (0, 0, 0))
                    self.assertEqual(proof["ihole"], {"K": 1, "L3": 4}[record["edge"]])
                    self.assertEqual(proof["nmix"], recipe.get("scf_mixing_iterations", 1))
                    self.assertEqual(proof["nscmt"], recipe["scf_iterations"])
                    self.assertAlmostEqual(proof["ca1"], recipe["scf_mixing"], places=7)
                    self.assertAlmostEqual(proof["rfms1"], recipe["scf_radius_A"], places=6)
                    self.assertEqual(sum(p["degen"] for p in report["selected_paths"]), expected_cn(record, index))
                    completed += 1
        if pending:
            self.skipTest(f"Validated {completed} completed sites; pending: {', '.join(pending)}")
        self.assertEqual(completed, sum(len(r["absorber_site_indices"]) for r in pilot.records()))

    def test_native_control_tampering_is_rejected_before_artifact_rehash(self):
        path = self.directory / "pot.json"
        original = json.loads(path.read_text())
        mutations = {"iunf": 1, "ihole": 1, "ixc": 1, "lfms1": 1,
                     "nscmt": original["nscmt"] + 1, "rfms1": original["rfms1"] + .1,
                     "nmix": original["nmix"] + 1, "ca1": original["ca1"] + .01}
        for field, value in mutations.items():
            path.write_text(json.dumps(dict(original, **{field: value})))
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "Native potential controls"):
                pilot.validate_site(self.record, self.directory)

    def test_corrupted_structural_metadata_is_rejected_before_rehashing(self):
        path = self.directory / "run.json"
        original = json.loads(path.read_text())
        corruptions = {"site_weight": original["site_weight"] / 2,
                       "site_multiplicity": original["site_multiplicity"] + 1,
                       "first_shell_distances_A": [r + .001 for r in original["first_shell_distances_A"]],
                       "actual_cluster_radius_A": original["actual_cluster_radius_A"] + .5}
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
    def test_candidate_preserves_all_374_base_rows_and_blobs_read_only(self):
        candidate = pilot.ROOT / "round2-candidates.db"
        supplied = os.environ.get("FEFFDB_ROUND2_BASE_DB")
        if not candidate.is_file() or not supplied:
            self.skipTest("Need candidate database and FEFFDB_ROUND2_BASE_DB for inheritance verification")
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
            self.assertEqual(old.execute("SELECT count(*) FROM feffdat").fetchone()[0], 374)
            self.assertEqual(manifest["upstream_paths"], 374)
            self.assertEqual(new.execute("SELECT count(*) FROM feffdat").fetchone()[0], manifest["total_paths"])
            self.assertEqual(manifest["upstream_paths"] + manifest["added_paths"], manifest["total_paths"])
            new_ids = {row["id"] for row in manifest["paths"]}
            self.assertEqual(len(new_ids), manifest["added_paths"])
            expected_paths = {(r["id"], i, p["filename"]): p["sha256"]
                              for r in pilot.records() for i in r["absorber_site_indices"]
                              for p in json.loads((archive(r, i) / "validation.json").read_text())["selected_paths"]}
            self.assertEqual(len(new_ids), len(expected_paths))
            self.assertEqual({row["reference"] for row in manifest["paths"]}, {r["id"] for r in pilot.records()})
            manifest_by_id = {row["id"]: row for row in manifest["paths"]}
            record_by_id = {r["id"]: r for r in pilot.records()}
            for ident in new_ids:
                self.assertIsNone(old.execute("SELECT id FROM feffdat WHERE id=?", (ident,)).fetchone())
                self.assertEqual(new.execute("SELECT status FROM drxas_path_review WHERE feffdat_id=?", (ident,)).fetchone(), ("candidate",))
                row = new.execute("SELECT d.feffdat,i.inptext,c.ciftext,r.reference_id,r.source_json,r.validation_json,c.source_db,c.source_id "
                                  "FROM feffdat d JOIN feffinp i ON d.feffinp_id=i.id "
                                  "JOIN cif c ON i.cif_id=c.id JOIN drxas_path_review r ON r.feffdat_id=d.id "
                                  "WHERE d.id=?", (ident,)).fetchone()
                path_blob, input_blob, cif_blob, reference_id, source_json, validation_json, source_db, source_id = row
                evidence = json.loads(validation_json)
                record = record_by_id[reference_id]
                self.assertEqual((source_db, source_id), pilot.source_identity(record))
                run, selected = evidence["run"], evidence["path"]
                directory = archive(record, run["site_index"])
                raw_path = lzma.decompress(path_blob)
                self.assertEqual(raw_path, (directory / selected["filename"]).read_bytes())
                self.assertEqual(hashlib.sha256(raw_path).hexdigest(), manifest_by_id[ident]["sha256"])
                self.assertEqual(manifest_by_id[ident]["reference"], reference_id)
                self.assertEqual(manifest_by_id[ident]["filename"], selected["filename"])
                self.assertEqual(lzma.decompress(input_blob), (directory / "feff.inp").read_bytes())
                self.assertEqual(lzma.decompress(cif_blob), (pilot.ROOT / "references" / reference_id / "source.cif").read_bytes())
                self.assertEqual(json.loads(source_json), json.loads((pilot.ROOT / "references" / reference_id / "provenance.json").read_text()))
                self.assertEqual(evidence["validation"]["review_status"], "candidate")
                self.assertFalse(evidence["validation"]["experimental_validation"])
                self.assertIn("native_potential_controls", evidence["validation"])
        self.assertEqual({p: pilot.sha(p) for p in before}, before)


if __name__ == "__main__":
    unittest.main()
