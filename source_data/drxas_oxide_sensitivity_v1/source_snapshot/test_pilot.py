"""Scientific ingestion guardrails; no native FEFF process or database writes."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import pilot


class CoordinationTests(unittest.TestCase):
    def setUp(self):
        # A 0.005-A splitting is physically distinct from FEFF output rounding.
        self.expected = [2.005] * 4 + [2.000] * 2
        self.paths = [dict(filename="feff0001.dat", reff=2.000, degen=2),
                      dict(filename="feff0002.dat", reff=2.005, degen=4)]

    def test_preserves_close_split_shells_and_their_multiplicities(self):
        self.assertEqual(pilot.match_coordination(self.expected, self.paths[::-1]), 0.0)

    def test_accepts_only_native_output_rounding(self):
        expected = [2.00004] * 2 + [2.00496] * 4
        self.assertAlmostEqual(pilot.match_coordination(expected, self.paths), 0.00004)

    def test_missing_or_duplicate_path_cannot_change_coordination(self):
        for paths in (self.paths[:1], self.paths + [self.paths[0]]):
            with self.subTest(paths=paths), self.assertRaisesRegex(ValueError, "Incomplete coordination"):
                pilot.match_coordination(self.expected, paths)

    def test_averaging_split_shells_is_rejected_even_with_correct_total_cn(self):
        merged = [dict(reff=(2 * 2.000 + 4 * 2.005) / 6, degen=6)]
        with self.assertRaisesRegex(ValueError, "distance mismatch"):
            pilot.match_coordination(self.expected, merged)

    def test_wrong_path_distance_is_rejected_with_correct_total_cn(self):
        paths = copy.deepcopy(self.paths)
        paths[1]["reff"] += 0.001
        with self.assertRaisesRegex(ValueError, "distance mismatch"):
            pilot.match_coordination(self.expected, paths)

    def test_nonintegral_nonpositive_and_nonfinite_degeneracy_are_rejected(self):
        for degen in (1.5, 0, -1, float("nan"), float("inf")):
            with self.subTest(degen=degen), self.assertRaisesRegex(ValueError, "degeneracy"):
                pilot.match_coordination([2.0], [dict(reff=2.0, degen=degen)])

    def test_empty_coordination_is_not_a_success(self):
        with self.assertRaisesRegex(ValueError, "Invalid CIF neighbor distances"):
            pilot.match_coordination([], [])

    def test_nonfinite_or_nonpositive_distances_cannot_pass(self):
        for value in (float('nan'), float('inf'), 0, -1):
            with self.subTest(value=value), self.assertRaises(ValueError):
                pilot.match_coordination([2.0], [dict(reff=value, degen=1)])
            with self.subTest(expected=value), self.assertRaises(ValueError):
                pilot.match_coordination([value], [dict(reff=2.0, degen=1)])


POTENTIAL_INPUT = """TITLE deliberately nonstoichiometric finite sphere
POTENTIALS
* ipot Z tag lmax1 lmax2 xnatph
0 73 Ta -1 -1 99
1 73 Ta -1 -1 99
2 8 O -1 -1 99
ATOMS
* x y z ipot tag
0 0 0 0 Ta
1 0 0 1 Ta
0 1 0 1 Ta
0 0 1 1 Ta
2 0 0 2 O
END
"""


class PotentialTests(unittest.TestCase):
    def test_stoichiometry_comes_from_cell_not_finite_sphere_counts(self):
        result = pilot.complete_potentials(POTENTIAL_INPUT, {"Ta": 4, "O": 10})
        rows = {int(row[0]): row for _, row in pilot.card_rows(result, "POTENTIALS")}
        self.assertEqual({i: float(row[5]) for i, row in rows.items()}, {0: 0.01, 1: 4, 2: 10})
        self.assertEqual(list(pilot.card_rows(result, "ATOMS")),
                         list(pilot.card_rows(POTENTIAL_INPUT, "ATOMS")))

    def test_absorber_without_same_element_background_keeps_cell_weight(self):
        text = "POTENTIALS\n0 77 Ir -1 -1\n1 8 O -1 -1\nATOMS\n0 0 0 0 Ir\n1 0 0 1 O\nEND\n"
        result = pilot.complete_potentials(text, {"Ir": 2, "O": 4})
        self.assertEqual([float(row[5]) for _, row in pilot.card_rows(result, "POTENTIALS")], [2, 4])

    def test_atoms_must_use_known_and_present_potentials(self):
        for text in (POTENTIAL_INPUT.replace("2 0 0 2 O", "2 0 0 3 O"),
                     POTENTIAL_INPUT.replace("2 0 0 2 O\n", "")):
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, "ATOMS/POTENTIALS mismatch"):
                pilot.complete_potentials(text, {"Ta": 4, "O": 10})

    def test_missing_element_potential_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "CIF composition"):
            pilot.complete_potentials(POTENTIAL_INPUT, {"Ta": 4, "O": 10, "H": 1})

    def test_nonconsecutive_potentials_are_rejected(self):
        text = POTENTIAL_INPUT.replace("2 8 O", "3 8 O")
        with self.assertRaisesRegex(ValueError, "nonconsecutive potentials"):
            pilot.complete_potentials(text, {"Ta": 4, "O": 10})

    def test_two_background_potentials_for_same_element_need_explicit_weights(self):
        text = POTENTIAL_INPUT.replace("ATOMS\n", "3 8 Osecond -1 -1 99\nATOMS\n")
        with self.assertRaisesRegex(ValueError, "Multiple background potentials"):
            pilot.complete_potentials(text, {"Ta": 4, "O": 10})

    def test_invalid_cell_amounts_are_rejected(self):
        for amount in (0, -1, float("nan"), float("inf")):
            with self.subTest(amount=amount), self.assertRaisesRegex(ValueError, "stoichiometry"):
                pilot.complete_potentials(POTENTIAL_INPUT, {"Ta": 4, "O": amount})


class SourceIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="pilot-unit-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.record = copy.deepcopy(next(r for r in pilot.records()
                                         if r["cod_id"] == "2101854"))
        self.path = self.root / "references" / self.record["id"] / "source.cif"
        self.path.parent.mkdir(parents=True)

    def copy_original(self):
        source = pilot.ROOT / "references" / self.record["id"] / "source.cif"
        shutil.copyfile(source, self.path)

    def test_fake_cif_with_wrong_hash_is_refused_before_parsing(self):
        self.path.write_text("This is intentionally not a CIF.\n")
        with self.assertRaisesRegex(ValueError, "CIF hash mismatch"):
            pilot.load_structure(self.record, self.root)

    def test_real_cif_with_changed_bytes_is_refused(self):
        self.copy_original()
        self.path.write_bytes(self.path.read_bytes() + b"\n# changed after source review\n")
        with self.assertRaisesRegex(ValueError, "CIF hash mismatch"):
            pilot.load_structure(self.record, self.root)

    def test_real_cif_with_tampered_recorded_hash_is_refused(self):
        self.copy_original()
        self.record["cif_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "CIF hash mismatch"):
            pilot.load_structure(self.record, self.root)

    def test_verified_real_reference_has_one_ir_site_in_ordered_six_atom_cell(self):
        self.copy_original()
        text, cluster, sites, _ = pilot.load_structure(self.record, self.root)
        self.assertEqual(text.encode(), self.path.read_bytes())
        self.assertTrue(cluster.struct.is_ordered)
        self.assertEqual(len(cluster.struct), 6)
        self.assertEqual(cluster.struct.composition.get_el_amt_dict(), {"Ir": 2, "O": 4})
        self.assertEqual([index for index, *_ in sites], [1])
        self.assertEqual([int(multiplicity) for _, _, multiplicity, _ in sites], [2])

    def test_changed_absorber_inventory_requires_review(self):
        self.copy_original()
        self.record["absorber_site_indices"] = [2]
        with self.assertRaisesRegex(ValueError, "Absorber site inventory changed"):
            pilot.load_structure(self.record, self.root)

    def test_partial_occupancy_cannot_be_normalized_into_an_ordered_reference(self):
        self.record = copy.deepcopy(next(r for r in pilot.records() if r['cod_id'] == '1530633'))
        path = self.root/'references'/self.record['id']/'source.cif'
        path.parent.mkdir(parents=True)
        original = (pilot.ROOT/'references'/self.record['id']/'source.cif').read_text()
        changed = original.replace('O1 O-2 0.267 0.35 0 1 0.0', 'O1 O-2 0.267 0.35 0 0.5 0.0')
        self.assertNotEqual(original, changed)
        path.write_text(changed)
        self.record['cif_sha256'] = pilot.sha(path)
        with self.assertRaisesRegex(ValueError, 'Disordered/partial occupancy'):
            pilot.load_structure(self.record, self.root)


class ArchivedRunTests(unittest.TestCase):
    def setUp(self):
        self.record = pilot.records()[0]
        source = pilot.ROOT/'calculations'/self.record['id']/'site-1'
        if not source.exists():
            self.skipTest('Run the pilot to create calculation fixtures')
        self.temp = tempfile.TemporaryDirectory(prefix='pilot-run-test-')
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)/'site-1'
        shutil.copytree(source, self.directory)

    def test_archived_calculation_is_usable_by_larch(self):
        result = pilot.validate_site(self.record, self.directory)
        self.assertEqual(result['coordination_sum'], 6)
        self.assertEqual(result['larch_chi_evaluation'], 'passed')
        self.assertEqual(result['review_status'], 'candidate')
        self.assertFalse(result['experimental_validation'])

    def test_failed_missing_or_unsuccessful_module_is_not_a_valid_archive(self):
        path = self.directory/'run.json'
        original = json.loads(path.read_text())
        for kind in ('failed', 'missing', 'nonzero'):
            record = copy.deepcopy(original)
            if kind == 'failed':
                record['execution_status'] = 'failed'
            elif kind == 'missing':
                record['module_results'].pop()
            else:
                record['module_results'][-1]['returncode'] = 1
            pilot.write_json(path, record)
            with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, 'six modules'):
                pilot.validate_site(self.record, self.directory)

    def test_changed_input_cannot_reuse_old_validation(self):
        path = self.directory/'feff.inp'
        path.write_text(path.read_text().replace('S02 1.0', 'S02 0.5'))
        with self.assertRaisesRegex(ValueError, 'Input has changed'):
            pilot.validate_site(self.record, self.directory)

    def test_missing_first_shell_output_is_not_silently_pruned(self):
        (self.directory/'feff0002.dat').unlink()
        with self.assertRaisesRegex(ValueError, 'Incomplete coordination'):
            pilot.validate_site(self.record, self.directory)

    def test_modified_raw_file_cannot_silently_get_a_new_checksum(self):
        path = self.directory/'feff0001.dat'
        path.write_bytes(path.read_bytes()+b'\n')
        with self.assertRaisesRegex(ValueError, 'Archived calculation bytes changed'):
            pilot.validate_site(self.record, self.directory)

    def test_native_potential_error_blocks_a_candidate(self):
        path = self.directory/'feff.log'
        path.write_text(path.read_text()+'\nResults might be meaningless\n')
        with self.assertRaisesRegex(ValueError, 'Fatal FEFF diagnostic'):
            pilot.validate_site(self.record, self.directory)

    def test_scf_iteration_limit_is_not_a_converged_candidate(self):
        path = self.directory/'feff.log'
        path.write_text(path.read_text()+'\nSCF ITERATION NUMBER 30 OUT OF 30\n')
        with self.assertRaisesRegex(ValueError, 'SCF reached the iteration limit'):
            pilot.validate_site(self.record, self.directory)


if __name__ == "__main__":
    unittest.main()
