# Oxide pilot sensitivity findings — 2026-11-16

All 15 native calculations completed all six modules successfully. Every run
retained the six-oxygen coordination multiset, with the original close split
paths. SCF stopped after 6–9 iterations, below the cap of 30. Three fresh
baselines reproduced the published inputs and all 36 raw FEFF path files byte
for byte; package versions and all six executable fingerprints also matched.

The larger-radius comparisons are locally stable under the predeclared numerical
review triggers. The original SCF=5 Å / cluster=8 Å recipe can remain the
candidate baseline. **All references remain candidate. This is neither a proof
of global convergence nor experimental validation.**

| Reference | SCF 5→6 Å: maximum path complex NRMS | Maximum amplitude NRMS | Maximum phase RMS | Cluster 8→9 Å: maximum path complex NRMS |
|---|---:|---:|---:|---:|
| Rutile IrO₂ | 2.3125% | 0.1281% | 0.023105 rad | 0.0528% |
| β-PtO₂ | 2.0460% | 0.1213% | 0.020412 rad | 0.0213% |
| B-Ta₂O₅ | 0.1575% | 0.1357% | 0.001102 rad | 0.0824% |

Each maximum is taken independently over the individual first-shell paths and
the reported k windows; the maxima need not occur on the same path/window.
For the complete shell at k=3–12 Å⁻¹, SCF 5→6 Å complex NRMS is 2.1428%, 1.8895%,
and 0.1181% for Ir, Pt, and Ta respectively. Ir/Pt differences mainly involve
phase, while their amplitude changes are much smaller. No fitted E0 or distance
adjustment was used to absorb these changes.

The only review trigger occurred in the **smaller SCF=4 Å β-PtO₂ control**:
phase RMS reached 0.035743 rad (3–8 Å⁻¹), versus the 0.03 rad trigger; complex
NRMS reached 3.6298%. Both split paths trigger phase review in the 3–8 and 3–12
windows. This argues against reducing the Pt SCF radius to 4 Å. The smaller
cluster=7 Å controls stayed below all triggers. Lower-radius controls are
sensitivity probes, not evidence that smaller settings should replace the baseline.

## Method and limits

- Five settings per reference: (SCF, cluster)=(5,8), (4,8), (6,8), (5,7), (5,9) Å.
  Only one radius varies at a time. RPATH=4 Å, NLEG=2, explicit unit-cell potential
  stoichiometry, SCF mixing=0.2, and all other pilot settings remain fixed.
- Native-degeneracy complex Larch paths are reconstructed as
  `Z = -chi_imag + i*chi`, with S02=1, sigma2=0, E0=0 and deltaR=0.
  Sigma2=0 deliberately avoids damping away high-k changes. Each full shell is
  the complex sum of its paths. This is numerical path comparison, not a sample fit.
- `||k²(Z_variant-Z_baseline)|| / ||k²Z_baseline||` is the complex relative RMS.
  Amplitude relative RMS uses the magnitudes with the same denominator. Circular
  phase RMS uses squared baseline k² amplitude as weight, excluding baseline
  amplitudes below 5% of their window maximum; this controls phase instability
  near nodes. Real chi(k) relative RMS is also provided in the machine-readable results.
- Windows: 3–8, 8–12, 12–16, 16–20, and combined 3–12 Å⁻¹ on a 0.05 Å⁻¹ grid.
  High-k windows test numerical behavior; they do not imply that an experimental
  sample has usable data to k=20.
- Predeclared review triggers: 5% complex RMS, 3% amplitude RMS, 0.03 rad phase RMS.
  These are explicit engineering review criteria, **not published universal
  convergence tolerances or experimental uncertainty**. A 0.03 rad shift would
  correspond roughly to 0.0025 Å at k=6 if modeled entirely as 2k deltaR; no such
  fitted distance or bias is claimed here.
- No simultaneous SCF+cluster increase, SCF>6 Å, cluster>9 Å, iteration-tolerance
  sweep, alternative exchange/core-hole model, thermal disorder, multiple
  scattering, or comparison with measured standards was evaluated.
- Native underflow/denormal and radius-counter (`ratx=0`) warnings remain archived.
  They were not suppressed or reclassified. The early SCF stop has no explicit
  native convergence flag, so iteration count alone is not a convergence claim.
- Phase/source caveats are unchanged, especially the high-pressure B-Ta₂O₅
  average structure with uncertain oxygen positions, defects/intergrowth and a
  short Ta–O bond. Numerical stability cannot establish the sample phase.

## Reproducibility and evidence

`source_snapshot/` preserves the exact 95-file pilot from Dr-XAS/feffdb commit
`547463406c812c2ea70333cb3861d80be8ceeb6c`. The source inventory is checked before
execution and analysis. Each new run under `calculations/` stores the input,
all raw FEFF paths, native logs, versions, executable SHA-256 fingerprints,
module exit codes and validation. The parent source and XAScraft source were not edited.

- `sensitivity.py`: standalone runner and analysis.
- `protocol.json`: settings and review triggers recorded before native execution.
- `results.json`: full path/shell metrics and native provenance summary.
- `metrics.csv`, `curves/*.csv`: numerical results and real/imaginary/amplitude/phase curves.
- `sensitivity-summary.png` and `.svg`: maximum-per-path window comparisons.
- `test_sensitivity.py`, `metric-tests.log`: seven independent analytic metric checks, passed.
- `README.md`: complete tabulated results; `artifact-sha256.json`: delivery inventory.

To reanalyse without native calculations: `python sensitivity.py analyse`.
To repeat native calculations, copy `sensitivity.py`, `test_sensitivity.py`,
`plot_summary.py`, the license and unchanged `source_snapshot/` into a fresh
directory, then run `python sensitivity.py run --workers 3`. Existing completed
runs are revalidated and reused, never overwritten. Run
`python plot_summary.py` after analysis to regenerate the figure.

Primary documentation distinguishes the SCF sphere from the input cluster and
supports inspecting the complex path components:
[FEFF SCF card](https://feff.phys.washington.edu/FEFF/Docs/feff8/feff84/feff84/SCF.html),
[FEFF potential/cluster cautions](https://feff.phys.washington.edu/feff/Docs/feff8/feff8web/node49.html),
[Larch FEFF path documentation](https://xraypy.github.io/xraylarch/xafs_feffpaths.html).
