# Dr-XAS absorber expansion, round 3

Local candidate source draft, 2026-11-17. This round has no publication commit and has not been published or deployed. Every new reference remains **candidate/review**, with **no experimental XAS validation**. The commit in `references.json/upstream_base` identifies the existing round-2 database used as input; it is not a citation for these new calculations.

The source manifest contains 14 absorber/edge records from 11 distinct structural models, requiring all 19 inequivalent absorber-site calculations. All 19 site calculations completed and passed native/technical validation, producing **90 selected paths and 179 archived native paths**. The separate candidate SQLite contains **570 paths (480 inherited + 90 new)**. All 26 source tests passed with no skips. Application integration is recorded separately after freezing this source package. `coverage-audit.json` records the observed validation artifacts and keeps source completion separate from app publication or deployment.

## What coverage means

Coverage counts an element as an **absorber** only when a retained path has that element in `metadata.absorber` and an explicit absorption edge. Its occurrence as a scatterer or in a source CIF does not establish absorber coverage. One element/edge reference does not cover all phases, oxidation states, coordinations, scattering pairs or experimental conditions.

The audited local pre-round-3 catalogue is version `5-da3f0b7fa940caf8`: 603 entries, 62 absorbing elements including review (464 available entries and 139 review entries). Only 39 absorbing elements have default-visible available entries. The catalogue's 43 scatterer elements are a separate metric. If all 14 new absorbers below pass final technical validation and are ingested as review, total candidate-inclusive absorber coverage becomes **76**, while default-visible absorber coverage stays **39**. Registration or successful calculation alone does not establish that app ingestion or remote deployment occurred.

| New absorber | Edge | Reference model | Inequivalent sites | Selected neighbors per site |
|---|---|---|---:|---|
| Li | K | LiF, NBS measured 26°C | 1 | F: 6 |
| Be | K | BeO, crystal no. 1 in air | 1 | O: 4, split 3+1 |
| B | K | hBN, Pease AA′ model at 35°C | 1 | N: 3 |
| C | K | diamond, 25°C | 1 | C: 4 |
| N | K | same hBN model | 1 | B: 3 |
| O | K | same BeO model | 1 | Be: 4, split 3+1 |
| F | K | same LiF model | 1 | Li: 6 |
| Rb | K | RbNO3 phase IV, 296 K before heating | 3 | O: 12/12/12 envelopes |
| I | K | molecular iodine, conflicting low-temperature metadata | 1 | I: 1 intramolecular neighbor |
| Cs | K | CsNO3 phase II, 296 K | 3 | O: selected nearest-12 envelopes |
| Os | L3 | rutile OsO2 | 1 | O: 6, split 2+4 |
| Tl | L3 | bixbyite Tl2O3 | 2 | O: 6/6 |
| Bi | L3 | elemental rhombohedral Bi, 298 K | 1 | Bi: 3+3 local envelope |
| Th | L3 | pure fluorite ThO2 endmember | 1 | O: 8 |

## Sources and unresolved evidence

Ten distinct models retain exact COD-distributed CIF bytes and their original notices. LiF is a locally authored, explicitly marked transcription of [Swanson and Tatge, NBS Circular 539, volume 1 (1953), pp. 61–62, Table 35](https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular539v1.pdf). The measured column is 26°C and a=4.0270 Å; the separate 25°C comparison gives 4.0269 Å and is not substituted. `transcription.json` records the original table values, units, page locations, absence of a printed lattice uncertainty and source-PDF hash. Rocksalt Fm-3m coordinates and full unit occupancies encode the explicitly reported stoichiometric model; they are not an additional experimental coordinate or occupancy refinement. No COD ID or DOI is invented for this derived CIF.

The hBN CIF directly cites Wyckoff's 1963 compilation. Its exact cell and AA′ B-over-N stacking were independently checked against [Pease's 1952 powder-diffraction study](https://doi.org/10.1107/S0365110X52001064); those two source identities remain distinct. The experimental specimen had stacking faults. The ordered periodic model does not simulate those faults, impurities, thermal layer puckering or arbitrary BN polytypes. Diamond uses current corrected COD9011997, with a=b=c=3.566986 Å and an explicit 25°C measurement; no older mistyped cell axis is retained.

RbNO3 uses the original 296 K phase-IV structure, not a heated branch from the same paper. CsNO3 uses ordered phase II at 296 K. Unknown numeric diffraction temperature or pressure stays null and is not inferred from synthesis, the phase name or a nominal ambient value.

For I2 COD4511285, the unchanged COD and original publisher supporting CIF state **100(2) K**, while the [primary article](https://doi.org/10.1021/cg5005159) states **110(1) K** for the matching cell. This remains unresolved: the normalized measurement-temperature field is null, both claims remain explicit, and no temperature correction is applied. FEFF uses the exact deposited nuclear coordinates, not optimized structures from the same paper. The source formula/Z bookkeeping and differing printed bond values are retained and explained in provenance. The selected CN1 intramolecular bond does not represent the longer crystal-packing contacts.

OsO2 follows the documented Boman-1970 → Baur/Khan-1971 table → AMCSD/COD chain, not substitution into another rutile oxide. Tl2O3 and Bi retain their own experimental refinements. ThO2 is a pure fluorite endmember with the source's limited precision and unknown numeric measurement conditions; mixed La occupancy is not imported. Exact sources, citations, hashes, parser observations, exclusions and measurement-versus-synthesis distinctions belong in each reference's `provenance.json` and the source-selection records.

## Coordination and spectral boundaries

Native degeneracy is per absorbing atom. Rb and Cs each require all three inequivalent sites with normalized populations 1/3. Tl sites have populations 1/4 and 3/4. An ensemble model must apply these weights explicitly; multiplying each site by its cell multiplicity or summing every site at full amplitude double-counts populations. Split distances remain separate and the full distance/degeneracy multiset is checked, including BeO's 3+1 shell.

Cs uses an explicitly chosen **nearest-12 oxygen envelope**, not an asserted unique chemical coordination shell. Its common 3.95 Å cutoff preserves each actual distance. The per-site twelfth-to-thirteenth O gaps are 0.137535, 0.462309 and 0.264545 Å. For site 3, O neighbor 13 at 3.972503 Å is inside native RPATH=4 Å but outside the selected envelope. Archiving all generated paths and selecting a database shell are different operations. The longer Rb-O contacts and Bi's shorter-plus-longer 3+3 envelope likewise must not be called equal bonds.

Li through F are low-energy K-edge candidate scattering references. Core excitons, charge transfer, nonspherical potentials, band structure and anisotropy can matter near these edges; this workflow does not validate XANES or excitonic peaks. [Primary BSE work on diamond, hBN and lithium halides](https://doi.org/10.1103/PhysRevB.106.155133) illustrates why a valid structural model is insufficient to establish near-edge accuracy. Native SCF completion or a finite chi(k) is not experimental validation.

The manifest records tabulated higher absorption edges of all constituent elements within the calculated k range, using XrayDB energies and Larch `etok`. Examples are B K → N K in BN at k≈7.632 Å^-1, Be K → O K at 10.643 Å^-1, Li K → F K at 12.981 Å^-1, and Os L3 → L2 at 19.934 Å^-1. These are interference landmarks, not automatic clipping limits or recommended experimental kmax. Chemical shifts, noise, sample composition and the actual analysis window still matter. I and Cs use their K edges here; their L3 alternatives have close L2 boundaries and would be separate references. The native grid to k=20 Å^-1 is never a claim of experimental usable data to that limit.

## Recipe and validation levels

The default recipe is FEFF8L, isotropic EXAFS, cluster radius 8 Å, SCF radius 5 Å, RPATH 4 Å, NLEG 2, CRITERIA 0/0, kmax 20 Å^-1, S02=1 and EXCHANGE=0. SCF uses 30 iterations, mixing 0.2 and nmix=1 unless a separately recorded supported retry explicitly changes nmix. Actual native controls, recipe, package versions and executable hashes are retained for every run. Default frozen atomic f density is retained (`iunf=0`, no UNFREEZEF). For Th this does not validate chemical valence, actinide correlations, magnetic behavior, multiplets or XANES; consult the [FEFF8 manual](https://feff.phys.washington.edu/feff/Docs/feff8/feff82.pdf).

The evidence has three distinct levels:

1. **Source and geometry:** exact raw/derived CIF hashes, experimental source identity, ordered composition, all absorber sites, multiplicities and complete selected-neighbor distances; unknown or conflicting conditions remain explicit.
2. **Native and numerical:** regenerated input and structural metadata must match; all six native modules must exit successfully; logs must contain SCF evidence without reaching the iteration limit or fatal diagnostics. `pot.json` controls (`iunf`, `ixc`, `ihole`, `rfms1`, `nscmt`, `lfms1`, `ca1`, `nmix`) must match the recipe. Path arrays must be finite, have the expected seven columns, monotonic 0–20 k and positive mean-free paths. Selected distance/degeneracy multisets must match the source within 0.0002 Å and produce finite nonzero Larch chi(k).
3. **Experimental applicability:** not established. No measured-standard fit, independent radius-sensitivity study, charge-state test or sample-specific XANES validation is supplied for this round. Earlier oxide tests do not validate these new references. SCF early stopping is recorded without inventing a convergence flag.

Passing levels 1 and 2 does not promote a record beyond candidate/review. Source displacement parameters are not a validated EXAFS sigma2. All native paths out to RPATH, including unselected paths, remain archived with the input, module outputs, logs, fingerprints and technical-validation report.

## Remaining absorber gaps

Even after all 14 proposed absorbers are successfully ingested, **42 elements remain without absorber coverage**: H, He, Ne, Ar, Kr, Tc, Xe, Pm, Po, At, Rn, Fr, Ra, Ac, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts and Og.

H/He have no agreed source and core-edge validation workflow in this batch; they are not treated as ordinary hard-X-ray EXAFS gaps. Ne/Ar/Kr/Xe require measured solid structures with explicit cryogenic or pressure conditions. UO2 COD2310001 has an explicit 993 K cell but unresolved primary oxide-measurement context in a beta-uranium study, so it is deferred rather than labeled ambient UO2. The other listed radioactive or synthetic elements have no individually verified source/edge/native-validation plan in this round; Rn/Og also need noble-gas phase assessment. These are documented research gaps, not claims that all such calculations are impossible. Th's inclusion does not establish coverage or accuracy for the remaining actinides. No statement that the entire periodic table is filled is warranted.

## Reproduction and packaging

From the source-package directory, with the required scientific Python environment:

```sh
python pilot.py validate
python -m unittest discover -s . -p 'test_*.py'
python pilot.py run --output /tmp/round3-rerun
python pilot.py build-db \
  --base ../drxas_element_round2_v1/round2-candidates.db \
  --calculations /tmp/round3-rerun --output /tmp/round3-candidates.db
```

The pinned round-2 base SQLite has 480 path rows; the new separate candidate has 570. These differ from the combined XAScraft catalogue, which has 603 entries before ingestion. Successful ingestion of all 90 new review paths would give 693 entries, 159 absorber/scatterer pairs (previously 145), 76 absorbing elements, 229 review entries and the same 464 default-visible entries. `build-db` must open that exact quiescent base read-only and write a separate new candidate file; inherited rows and payload bytes must remain unchanged. General FeffDB consumers do not enforce the candidate review table, so this file is not a validated production replacement.

Final packaging must preserve exact source and calculation bytes and record the complete artifact inventory. Include source CIFs, factual transcription/provenance, code/tests and native evidence; exclude research journal PDFs or page images, search/cache directories, Python/pytest caches and native executables. The source license and original data notices remain separate from journal-article rights. A future publication SHA, app ingestion result or deployment claim must be recorded only after that action is explicitly performed and verified.
