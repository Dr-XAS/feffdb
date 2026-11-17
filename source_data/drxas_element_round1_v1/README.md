# Dr-XAS absorbing-element expansion, round 1

Local candidate source package, created 2026-11-16 (America/Chicago). No publication commit is assigned to this draft.

Nine ordered experimental crystal structures cover ten previously absent absorbing elements at the K edge: Na, Mg, Al, Si, P, S, Cl, K, Sc and Br. NaCl is calculated from both Na and Cl centers; Sc2O3 has two inequivalent Sc centers. The 11 site calculations produced 16 selected first-shell paths. Each calculation retains all generated single-scattering paths out to RPATH, while the candidate database includes the complete explicitly selected first coordination shell for every absorber site.

| Absorbing element | Structure / COD | Selected neighbor | Expected CN per site |
|---|---|---|---|
| Na, Cl | NaCl halite / 9003308 | Cl, Na respectively | 6 |
| Mg | MgO periclase / 1000053 | O | 6 |
| K | KCl sylvite / 9003112 | Cl | 6 |
| Br | NaBr / 9007464 | Na | 6 |
| Al | Al2O3 corundum / 9007498 | O | 6, split 3+3 |
| Si | alpha quartz / 9005017 | O | 4 |
| P | AlPO4 berlinite / 9010367 | O | 4 |
| S | ZnS sphalerite / 9000107 | Zn | 4 |
| Sc | Sc2O3 bixbyite / 4501808 | O | 6 at each of two sites |

## Source boundaries

Every raw CIF is unchanged COD-distributed data with its original notices. Adjacent provenance.json stores SHA256, article citation, verified DOI or original article URL, download URL, measurement conditions, source limitations, occupancy interpretation and symmetry-expanded site inspection. A missing DOI or measurement temperature/pressure remains unknown. The source-selection files record excluded KBr historical-unit records, mixed-composition sphalerite, the wrong-temperature quartz entry and historical Sc structures with assumed oxygen coordinates or an unresolved erratum.

Sc2O3 uses a modern 296(2) K single-crystal refinement (Cong et al., 2010, DOI 10.1021/cg1004962). Its original legacy multiplicity column is inconsistent with symmetry-expanded multiplicities; the unchanged CIF is retained, and the expanded site counts 8:24 determine absorber population weights 0.25:0.75. Each FEFF degeneracy is per absorber; an ensemble fit must apply normalized site populations explicitly. Adding both sites at full amplitude would double-count the population. Element-only background potentials are an approximation for inequivalent sites and are not an additional site-average calculation.

The ZnS paper measured the cubic cell edge at 25 C; its zincblende atomic coordinates are fixed by symmetry. Some source files omit an occupancy column: ordered full occupancy is the crystallographic parser interpretation, not a separately measured occupancy. This package does not change coordinates, normalize mixed occupancy, invent a textbook structure, or infer measurement conditions from synthesis.

## Recipe and review status

The pilot.py runner is adapted from the published Dr-XAS oxide pilot. It uses FEFF8L, K edges, isotropic/unpolarized EXAFS, S02=1, EXCHANGE=0, kmax=20 inverse angstrom, cluster radius=8 angstrom, RPATH=4 angstrom, NLEG=2, CRITERIA=0/0, and SCF 5.0 0 30 0.2 1. Raw files preserve the native core-hole and AFOLP defaults. Element potential stoichiometries come from the expanded unit cell, not sphere populations. Unsupported fixed recipe values are rejected rather than written only to metadata.

All six native modules must succeed. Every source and input is fingerprinted; package versions, executable hashes, native stdout/stderr, full path files and per-site technical validation are archived. Re-validation regenerates the FEFF input and all structural metadata independently from the pinned CIF. It compares the full first-shell distance multiset with FEFF native degeneracies (0.0002 angstrom tolerance), preserving close split shells, and checks finite seven-column arrays, increasing k, geometry, absorber/edge, native diagnostics and finite nonzero Larch chi(k).

Successful execution and early SCF termination do not prove convergence. **Every new path remains candidate/review.** Radius sensitivity and measured-reference validation have not been performed on this new round. The separate oxide sensitivity study applies only to its three oxide L3 references. A computed 0-20 k grid is not a measured usable range, a validation of soft-X-ray EXAFS, or an XANES calculation. No sample data are fitted here.

## SQLite and application use

round1-candidates.db is a separate evaluation copy of the pinned 358-path oxide candidate database. Every pre-existing row and BLOB is checked unchanged. The resulting candidate has 374 rows: 348 upstream, 10 earlier oxide candidates and 16 new first-shell candidates. New rows alone are listed in round1-candidates.manifest.json, and drxas_path_review carries their candidate status and complete evidence. The original upstream 348-row feffdb/feffdat.db and the earlier oxide pilot are unchanged. Generic upstream consumers do not enforce this review table: this is not a drop-in validated database release.

The XAScraft importer scripts/merge_feffdb_round1.py validates the database, manifest, exact full artifact inventory and complete shell evidence, then maps new rows to review. Existing paths and scientific metadata are preserved. A deterministic source ZIP pins this local draft by content hash; no old repository commit is falsely cited as containing these new calculations. A publication commit can be pinned after a separately authorized source release.

## Reproduce

Use an environment with xraylarch/FEFF8L, larixite, pymatgen, numpy and xraydb. Exact native and Python dependency fingerprints are in each run.json. No shared package was installed or changed.

```sh
python source_data/drxas_element_round1_v1/pilot.py validate
python -m unittest discover -s source_data/drxas_element_round1_v1 -p 'test_*.py'
python source_data/drxas_element_round1_v1/pilot.py run --output /tmp/element-round1-rerun
python source_data/drxas_element_round1_v1/pilot.py build-db \
  --base source_data/drxas_oxide_pilot_v1/pilot-candidates.db \
  --calculations /tmp/element-round1-rerun --output /tmp/round1-candidates.db
```

Existing output directories are never overwritten. artifact-sha256.json inventories all delivered files except itself and interpreter caches. Scientific sources retain their original attribution and COD notices; the existing FeffDB license is copied alongside, and FEFF binaries are not redistributed.
