"""Analytic checks for comparison metrics; never executes native FEFF."""
import unittest
import numpy as np
from sensitivity import metrics


class MetricTests(unittest.TestCase):
    def setUp(self):
        self.k = np.linspace(0, 20, 401)
        self.signal = np.exp(1j * self.k)

    def test_identical_complex_signal_is_zero(self):
        result = metrics(self.k, self.signal, self.signal, (3, 12))
        for key in ('complex_relative_rms', 'amplitude_relative_rms', 'phase_weighted_rms_rad'):
            self.assertLess(result[key], 1.e-14)
        self.assertEqual(result['review_flags'], [])

    def test_ten_percent_amplitude_has_known_relative_difference(self):
        result = metrics(self.k, self.signal, 1.1*self.signal, (3, 12))
        self.assertAlmostEqual(result['complex_relative_rms'], .1)
        self.assertAlmostEqual(result['amplitude_relative_rms'], .1)
        self.assertLess(result['phase_weighted_rms_rad'], 1.e-14)
        self.assertIn('amplitude_relative_rms', result['review_flags'])

    def test_constant_phase_has_known_complex_difference(self):
        phase = .1
        result = metrics(self.k, self.signal, self.signal*np.exp(1j*phase), (3, 12))
        self.assertAlmostEqual(result['complex_relative_rms'], 2*np.sin(phase/2))
        self.assertAlmostEqual(result['phase_weighted_rms_rad'], phase)
        self.assertLess(result['amplitude_relative_rms'], 1.e-14)

    def test_two_pi_phase_is_not_a_physical_change(self):
        result = metrics(self.k, self.signal, self.signal*np.exp(2j*np.pi), (3, 12))
        self.assertLess(result['phase_weighted_rms_rad'], 1.e-14)
        self.assertEqual(result['review_flags'], [])

    def test_window_excludes_changes_outside_requested_range(self):
        changed = self.signal.copy()
        changed[self.k > 12] *= 10
        result = metrics(self.k, self.signal, changed, (3, 12))
        self.assertEqual(result['complex_relative_rms'], 0)

    def test_phase_at_vanishing_baseline_does_not_dominate(self):
        low = self.k < 6
        baseline = self.signal.copy()
        baseline[low] *= 1.e-6
        changed = baseline.copy()
        changed[low] *= -1
        result = metrics(self.k, baseline, changed, (3, 12))
        self.assertLess(result['phase_weighted_rms_rad'], 1.e-14)
        self.assertLess(result['phase_points'], result['points'])

    def test_zero_reference_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'Zero reference'):
            metrics(self.k, self.signal*0, self.signal, (3, 12))


if __name__ == '__main__':
    unittest.main()
