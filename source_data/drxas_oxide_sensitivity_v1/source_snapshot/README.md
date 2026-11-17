# Dr-XAS oxide reference pilot v1

Three experimental crystal structures supply ten new **candidate** first-shell
oxygen single-scattering paths at the Ir, Pt and Ta L3 edges. This is a local
curation pilot, not a validated reference release or an XAScraft deployment.

The original `feffdb/feffdat.db` remains unchanged. `pilot-candidates.db` is a
separate, explicitly named evaluation database: 348 preserved upstream paths plus
10 candidates = 358 paths. Its extra `drxas_path_review` table carries status,
source citation, native diagnostics and validation evidence. Upstream FeffDB
clients do not enforce this additional status table: do not substitute this file
for the production database before scientific review.

| Reference | COD | Edge | Selected FEFF R_eff (angstrom), degeneracy | Total O CN |
|---|---|---|---|---|
| Rutile IrO2 | 2101854 | Ir L3 | 1.9604 (2), 1.9988 (4) | 6 |
| beta-PtO2, CaCl2 type | 1530633 | Pt L3 | 1.9891 (2), 2.0034 (4) | 6 |
| B-Ta2O5, high-pressure polymorph | 2102123 | Ta L3 | 1.6699, 1.9030, 1.9085, 2.1049, 2.1688, 2.2405 (1 each) | 6 |

Each structure has one inequivalent metal site. Site multiplicity is provenance;
it is **not** an additional multiplier of FEFF degeneracy. The oxygen selection
cutoff is explicitly reviewed as 2.5 angstrom for these three six-coordinate
polyhedra; it is not a general first-shell rule. All 36 generated SS paths out to
4 angstrom are retained in `calculations`; only the ten selected O paths enter the
candidate database. Nearby split shells are not deduplicated by distance.

## Sources and phase limits

- IrO2: Bolzan, Fong, Kennedy and Howard (1997), neutron powder refinement,
  [DOI 10.1107/S0108768197001468](https://doi.org/10.1107/S0108768197001468),
  [COD 2101854](https://www.crystallography.net/cod/2101854.html).
  Occupancy is implicit in this CIF, and symmetry is supplied by the declared
  space group. Measurement temperature and pressure are unknown here.
- beta-PtO2: Range, Rau, Klement and Heyns (1987), single-crystal refinement,
  [DOI 10.1016/0025-5408(87)90220-0](https://doi.org/10.1016/0025-5408(87)90220-0),
  [COD 1530633](https://www.crystallography.net/cod/1530633.html).
  The reported 40 kbar / 1500 C conditions describe synthesis, not diffraction.
  Measurement temperature and pressure are unknown here.
- B-Ta2O5: Zibrov, Filonenko, Sundberg and Werner (2000), powder refinement,
  [DOI 10.1107/S0108768100005462](https://doi.org/10.1107/S0108768100005462),
  [COD 2102123](https://www.crystallography.net/cod/2102123.html).
  Diffraction temperature is 293 K; synthesis was 8 GPa / 1470 K; measurement
  pressure is unspecified. The source sample had B/Z phase intergrowth and
  defects. This ordered average structure is not generic ambient L-Ta2O5.
  The short 1.67(6)-angstrom Ta-O distance and uncertain oxygen coordinates are
  retained, not repaired. COD documents its correction of malformed occupancy
  fields: these are unchanged COD-distributed bytes, not claimed identical to the
  original IUCr supplement.

The source CIFs and their original notices are preserved byte-for-byte under
`references/*/source.cif`. Each adjacent `provenance.json` records download URL,
retrieval date, SHA-256, article citation, conditions and caveats. COD's public
data policy and the original attribution notices are retained. Upstream code and
database retain their existing license; this pilot does not redistribute FEFF
executables or claim ownership of the source structures.

Deferred references are documented in `source-selection.json`: historical PtO
has a reported reproducibility concern; alpha-PtO2 has coarse coordinates;
Z-Ta2O5 needs a separate seven-coordinate-shell interpretation. No idealized
replacement, occupancy normalization or structure substitution was made.

## Calculation and validation

Recipe: FEFF8L, L3, isotropic/unpolarized, S02=1, Hedin-Lundqvist EXCHANGE=0,
EXAFS kmax=20 inverse angstrom, cluster radius=8 angstrom, RPATH=4 angstrom,
NLEG=2, CRITERIA=0/0, and `SCF 5.0 0 30 0.2 1`. Native core-hole and AFOLP defaults
are retained and echoed in the raw path headers. This recipe has not been
established as a converged or universally transferable reference calculation.

`pilot.py` uses larixite's CIF conversion with a fixed RNG seed and removes the
wall-clock input header. It explicitly fills every potential's stoichiometry
from the symmetry-expanded cell, checks the actual sphere, and rejects partial
or mixed occupancy. All six native modules must finish successfully. Inputs,
path files, stdout/stderr, package versions, executable hashes and module exit
codes are archived separately for every absorber site.

The archived runs stopped SCF iterations at 6 (Ir), 9 (Pt) and 8 (Ta), below the
30-iteration limit. The log has no explicit SCF-convergence flag; early native
termination is recorded without claiming radius or experimental convergence.

Technical checks passed:

- Full CIF oxygen-distance multiset matches FEFF degeneracies and distances;
  maximum difference is below 0.00005 angstrom (tolerance 0.0002 angstrom).
- All 36 raw paths have valid geometry, finite seven-column arrays and increasing
  k from 0 through 20 inverse angstrom. All ten candidates evaluate to finite,
  nonzero chi(k) with Larch; the smoke-test sigma2 is not an experimental value.
- A fresh IrO2 repeat produced byte-identical input and all six raw FEFF files.
  `reproducibility.json` records hashes. This tests repeatability, not changing
  the SCF/cluster radii or comparison with measured reference spectra.
- 28 unit/archive tests cover close split shells, missing/duplicate paths,
  wrong distances/occupancies/source hashes, stoichiometry, input mutation,
  incomplete module execution, altered archived bytes and SCF iteration-limit refusal.
- Candidate SQLite integrity and foreign keys pass; every pre-existing row and
  BLOB payload is preserved. Decompressed new files match the raw source hashes.

The original FeffDB Python package could not be imported in the available
environment: installed `pyshortcuts 1.9.9` lacks `hash_password` and
`check_password_rules`, which upstream FeffDB imports. This is recorded as a
consumer-environment blocker, not reported as a passing API test. No upstream
package code or shared dependencies were changed. SQLite and Larch consumption
of the ten new database BLOBs were verified independently; see `verification.json`.
Also, upstream's `feffpath` helper defaults to degeneracy 1: when that helper is
used, pass the database's native degeneracy explicitly. This source-level finding
does not replace the blocked FeffDB API execution check.

Native warnings are deliberately retained: IEEE underflow/denormal flags and a
pathfinder radius warning reporting `ratx=0`. The input clusters actually extend
beyond 7.93 angstrom and cover RPATH=4. This disagreement is consistent with the
known FEFF8L radius-counter diagnostic issue; it is not proof of convergence.
No invalid/divide-by-zero/overflow diagnostic was accepted. Each candidate's
database metadata contains the warnings and CN evidence, not merely "passed".

Before promotion, review the phase/temperature/pressure match, examine native
diagnostics, test SCF/cluster-radius sensitivity, and compare with suitable
measured standards. Keep all three references marked `candidate` until then.

## Reproduce locally

Use Python with `larixite`, `pymatgen`, `xraylarch` (with FEFF8L executables),
`numpy` and `xraydb`. Exact versions and binary hashes are in `run.json`.
No XAScraft app import, network request or default user database is required.
The currently installed scientific environment was reused without installing or
changing its packages. Input comparison is deliberately strict: changed converter
versions require review rather than silently updating the reference.

From the repository root:

```sh
python source_data/drxas_oxide_pilot_v1/pilot.py validate
python -m unittest discover -s source_data/drxas_oxide_pilot_v1 -p 'test_*.py'

# New output directory; existing runs are never overwritten.
python source_data/drxas_oxide_pilot_v1/pilot.py run --output /tmp/oxide-rebuild
python source_data/drxas_oxide_pilot_v1/pilot.py build-db \
  --base feffdb/feffdat.db --calculations /tmp/oxide-rebuild \
  --output /tmp/oxide-candidates.db
```

The builder accepts only the pinned upstream SQLite SHA in `references.json`,
opens it read-only, rejects active WAL/SHM files, writes a separate candidate,
and verifies preservation before publishing that local output. The upstream
first-shell helper and its distance-based duplicate suppression are not used.
