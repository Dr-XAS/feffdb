# Dr-XAS element and coordination expansion, round 2

Local candidate source package, 2026-11-16 (America/Chicago). This draft has no publication commit. It extends the content-pinned 374-row candidate database from round 1 in a separate SQLite file; previous source packages and the original upstream database remain unchanged.

Nineteen experimental reference models cover thirteen new absorbing elements at L3 (La, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu) and six missing environments for existing elements (Ga-O, In-O, Hf-O, Al-N, Si-N and S-O). All 29 inequivalent absorber-site calculations completed, yielding 106 selected first-shell paths and 281 archived native paths. The separate candidate SQLite contains 480 paths (374 earlier + 106 new); all 24 source tests pass with zero skips. Each site retains its full selected single-scattering first-shell distance multiset and native FEFF degeneracies. All generated paths out to RPATH are archived, including paths not selected for the database.

## Sources and interpretation

Eighteen CIFs retain exact COD-distributed bytes and original attribution notices. AlN is an explicitly marked transcription of Schulz and Thiemann's 1977 experimental single-crystal refinement, DOI 10.1016/0038-1098(77)90959-0. Its delivered transcription.json records published values, uncertainties, page/table locations and independent checking. transcribe.py recreates the CIF. The original paper does not print a deposited CIF: the standard wurtzite P63mc encoding and unit occupancies are documented modeling encodings, not independently refined additional measurements. The resulting 3x1.889119 +1x1.902858 A shell agrees with the independently printed bond table. No ideal u=.375 substitution or lattice rescaling was made.

The initial Si-C target is deferred because available historical cells could not be verified reliably and newer candidates had doping or composition concerns. The replacement is measured beta-Si3N4 at 295 K (COD2102550; DOI10.1107/S010876810401393X), with Si-N scattering. Its source paper also includes ab-initio calculations; this CIF is the experimental refinement. Source Z=1.5 conflicts with expanded Si6N8 (Z=2); the original bytes are retained and FEFF stoichiometry uses the expanded composition. The 1900 K synthesis temperature is distinct from the 295 K diffraction temperature.

Er2O3 COD1534469 was measured at 514 K and Yb2O3 COD1534460 at 598 K, verified against the original 1985 paper tables. Their labels and metadata retain those temperatures. Unknown temperature/pressure stays unknown. B-type Sm2O3 and Tb2O3 are separate monoclinic phases; neither is labeled a generic ambient cubic oxide. Source exclusions and their reasons are in the source-selection JSON files. Original journal PDFs/screenshots are not redistributed; citations, URLs, source hashes and factual transcription evidence are retained.

## Sites, coordination and fit boundaries

- Ga2O3 has tetrahedral CN4 and octahedral CN6 Ga sites, with populations1/2 each.
- In/Eu/Ho/Er bixbyite sites have populations1/4 and3/4. In the selected Yb source the site ordering is reversed:3/4 then1/4. Sm/Tb have three equal-population sites with CN6/7/7 and7/7/6, respectively.
- Native degeneracy is per absorbing atom. An ensemble model must apply normalized population weights explicitly; adding inequivalent sites at full amplitude double-counts populations. Multiplicity is not an extra amplitude multiplier.
- LaF3 retains all11 selected F neighbors, including two longer contacts near2.99925 A. HfO2 keeps all7 O distances, including the close2.171033/2.172095 A pair. No distance-based deduplication removes split paths.
- Every L3 source records the next L2 edge from XrayDB and its equivalent k using Larch etok. For example, La encounters L2 near k10.35 inverse angstrom. These are interference markers, not recommendations for an exact experimental kmax. Chemical shifts, other edges, noise and the analysis window still require sample-specific review.

## Recipe and numerical evidence

The default recipe is unchanged from round1: FEFF8L, isotropic EXAFS, cluster8 A, SCF5 A, RPATH4 A, NLEG2, CRITERIA0/0, kmax20 inverse angstrom, S02=1, EXCHANGE0, SCF30 iterations/mixing0.2/nmix1. A validated optional nmix setting is supported for an explicit future retry; actual recipe values and native controls are recorded for every delivered calculation.

Native default frozen atomic f-electron density is retained (iunf=0; no UNFREEZEF card). This is not an independently validated rare-earth valence, multiplet, magnetic or XANES model. See the FEFF8 manual (https://feff.phys.washington.edu/feff/Docs/feff8/feff82.pdf), SCF guidance (https://feff.phys.washington.edu/feff/Docs/feff8/feff85/feff85/SCF.html) and FEFF8L scope (https://github.com/xraypy/feff85exafs).

Validation regenerates every input and structural metadata record from its pinned CIF; compares the full selected-shell distance/degeneracy multiset within0.0002 A; verifies finite seven-column path arrays, monotonic0-20 k, positive mean-free paths and finite nonzero Larch chi(k); and requires all six native modules to succeed. pot.json is retained and its actual iunf/ixc/ihole/rfms1/nscmt/lfms1/ca1/nmix values are checked against the recorded recipe and edge. Native warnings and executable/package fingerprints remain in the archive.

All new records stay candidate/review. Early SCF termination is recorded without claiming a native convergence flag. No radius-sensitivity study or measured-standard fit has yet been performed for this batch. The earlier three-oxide radius study does not validate these new references. Theoretical k20 is not an experimental usable range. The original source Uiso values, where transcribed, are not substituted for fitted EXAFS sigma2.

## Reproduction and ingestion

```sh
python source_data/drxas_element_round2_v1/pilot.py validate
python -m unittest discover -s source_data/drxas_element_round2_v1 -p 'test_*.py'
python source_data/drxas_element_round2_v1/pilot.py run --output /tmp/round2-rerun
python source_data/drxas_element_round2_v1/pilot.py build-db \
  --base source_data/drxas_element_round1_v1/round1-candidates.db \
  --calculations /tmp/round2-rerun --output /tmp/round2-candidates.db
```

Exact scientific package versions and native executable hashes are in run.json. The importer uses read-only SQLite and never imports the FeffDB Python package or runs FEFF. New manifest IDs alone are imported, existing rows/BLOBs and original app path bytes are verified unchanged, and every path remains review. XAScraft pins the complete deterministic evidence ZIP by SHA256; a future publication commit must be assigned explicitly rather than borrowing an older commit. General FeffDB consumers do not enforce drxas_path_review, so the candidate SQLite is not a validated production replacement.

The source license is copied unchanged from the existing fork. Raw scientific data retain their original notices; the AlN CIF is a documented derivation from published factual parameters. FEFF executables are not redistributed.
