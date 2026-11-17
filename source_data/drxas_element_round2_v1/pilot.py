#!/usr/bin/env python3
"""Reproducible, local-only absorbing-element reference expansion. See README.md for scope.

Requires larixite, xraylarch (including FEFF8L), pymatgen and numpy.
Does not import the XAScraft app or alter the bundled FeffDB database.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import lzma
import math
import os
import platform
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import warnings

ROOT = Path(__file__).resolve().parent
MODULES = ('rdinp', 'pot', 'xsph', 'pathfinder', 'genfmt', 'ff2x')
RECIPE = dict(cluster_radius_A=8.0, scf_radius_A=5.0, scf_iterations=30,
              scf_mixing=0.2, path_radius_A=4.0, max_legs=2, kmax_invA=20,
              s02=1.0, exchange=0, polarization=None, with_h=False, rng_seed=0)
DISTANCE_TOLERANCE_A = 0.0002  # coordinates 5 decimals; native reff 4 decimals


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def records(root=ROOT):
    return json.loads((root / 'references.json').read_text())['references']


def card_rows(text, name):
    active = False
    for i, line in enumerate(text.splitlines()):
        line = line.split('*', 1)[0].strip()
        if not line:
            continue
        if line.upper() == name:
            active = True
        elif active:
            if not re.match(r'^[+\-.0-9]', line):
                break
            yield i, line.split()


def complete_potentials(text, composition):
    """Explicit unit-cell ratios for every potential, never sphere atom counts."""
    from xraydb import atomic_number
    counts = {atomic_number(e): float(n) for e, n in composition.items()}
    rows = list(card_rows(text, 'POTENTIALS'))
    pots = {int(row[0]): int(row[1]) for _, row in rows}
    if not rows or set(pots) != set(range(len(rows))):
        raise ValueError('Missing or nonconsecutive potentials')
    if len(set(list(pots.values())[1:])) != len(rows) - 1:
        raise ValueError('Multiple background potentials for one element need explicit weights')
    if set(pots.values()) != set(counts):
        raise ValueError('POTENTIALS do not cover the CIF composition')
    if {int(row[3]) for _, row in card_rows(text, 'ATOMS')} != set(pots):
        raise ValueError('ATOMS/POTENTIALS mismatch')
    lines = text.splitlines()
    for i, row in rows:
        n = .01 if int(row[0]) == 0 and pots[0] in list(pots.values())[1:] else counts[int(row[1])]
        if len(row) < 5 or not math.isfinite(n) or n <= 0:
            raise ValueError('Invalid potential stoichiometry')
        lines[i] = '  ' + '  '.join([*row[:5], f'{n:.12g}', *row[6:]])
    return '\n'.join(lines) + '\n'


def source_identity(record):
    database = record.get('source_database', 'COD')
    identifier = record.get('source_record_id', record.get('cod_id'))
    if database not in ('COD', 'publication') or not isinstance(identifier, str) or not identifier.strip():
        raise ValueError('Invalid scientific source identity')
    return database, identifier


def load_structure(record, root=ROOT):
    from larixite.cif_cluster import CIF_Cluster
    from pymatgen.io.cif import CifParser
    path = root / 'references' / record['id'] / 'source.cif'
    if sha(path) != record['cif_sha256']:
        raise ValueError(f"CIF hash mismatch: {record['id']}")
    text = path.read_text()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        parser = CifParser.from_str(text, occupancy_tolerance=1.0)
        strict = parser.parse_structures(primitive=False, on_error='raise')[0]
        cluster = CIF_Cluster(ciftext=text, absorber=record['absorber'])
    # Mixed or partially occupied sites must not be sampled into an ordered model.
    if not strict.is_ordered or not cluster.struct.is_ordered:
        raise ValueError('Disordered/partial occupancy is excluded from this pilot')
    if strict.composition != cluster.struct.composition or len(strict) != len(cluster.struct):
        raise ValueError('Strict parser and generator disagree about the expanded cell')
    sites = [(i, site, mult, wyckoff) for i, (site, mult, wyckoff) in
             enumerate(cluster.unique_sites, 1) if site.specie.symbol == record['absorber']]
    if [i for i, *_ in sites] != record['absorber_site_indices']:
        raise ValueError('Absorber site inventory changed; explicit review needed')
    return text, cluster, sites, sorted(set(str(w.message) for w in caught))


def prepare(record, site_index, root=ROOT, recipe=None):
    from larixite.cif_cluster import cif2feffinp
    recipe = record.get('recipe', RECIPE) if recipe is None else recipe
    fixed = ('max_legs', 'kmax_invA', 's02', 'exchange', 'polarization', 'with_h', 'rng_seed')
    if set(recipe) - {'scf_mixing_iterations'} != set(RECIPE) or any(recipe[k] != RECIPE[k] for k in fixed):
        raise ValueError('Unsupported recipe: fixed scientific cards must match the runner')
    nmix = recipe.get('scf_mixing_iterations', 1)
    if not isinstance(nmix, int) or isinstance(nmix, bool) or not 1 <= nmix <= 30:
        raise ValueError('SCF mixing iterations must be an integer from 1 to 30')
    if type(recipe['scf_iterations']) is not int or not 0 < recipe['scf_mixing'] < 1 or not 1 < recipe['scf_iterations'] <= 100:
        raise ValueError('Invalid SCF iteration/mixing controls')
    cif, cluster, sites, parser_warnings = load_structure(record, root)
    _, site, mult, wyckoff = next(s for s in sites if s[0] == site_index)
    atoms = cluster.struct.get_sites_in_sphere(site.coords, recipe['cluster_radius_A'])
    if len(atoms) > 400:
        raise ValueError('Cluster exceeds the pilot 400-atom limit')
    text = cif2feffinp(cif, record['absorber'], edge=record['edge'],
                      absorber_site=site_index, cluster_size=recipe['cluster_radius_A'],
                      version8=True, with_h=False, rng_seed=0,
                      extra_titles=[f"Dr-XAS candidate {record['id']} {source_identity(record)[0]} {source_identity(record)[1]}"])
    text = re.sub(r'^\* FEFF8 input file generated .*$',
                  '* FEFF8 input generated by Dr-XAS element-round2-v1', text, flags=re.M)
    for card, value in [('EDGE', record['edge']), ('RPATH', recipe['path_radius_A']),
                        ('NLEG', 2), ('EXAFS', 20), ('S02', 1.0)]:
        text, n = re.subn(rf'^{card}\s+.*$', f'{card} {value}', text, flags=re.M)
        if n != 1:
            raise ValueError(f'Expected exactly one {card} card')
    text, n = re.subn(r'^\*SCF\s+.*$',
                     f"SCF {recipe['scf_radius_A']} 0 {recipe['scf_iterations']} {recipe['scf_mixing']} {nmix}",
                     text, flags=re.M)
    if n != 1:
        raise ValueError('SCF template changed')
    text = text.replace('NLEG 2\n', 'NLEG 2\nCRITERIA 0 0\n')
    composition = cluster.struct.composition.element_composition.reduced_composition.get_el_amt_dict()
    text = complete_potentials(text, composition)
    rows = list(card_rows(text, 'ATOMS'))
    extent = max(math.sqrt(sum(float(v)**2 for v in row[:3])) for _, row in rows)
    if extent < max(recipe['path_radius_A'], recipe['scf_radius_A']):
        raise ValueError('Actual cluster does not cover requested radii')
    neighbors = sorted(float(n.nn_distance) for n in cluster.struct.get_neighbors(site, recipe['path_radius_A'])
                       if n.specie.symbol == record['scatterer'])
    first = [r for r in neighbors if r < record['shell_cutoff_A']]
    if len(first) != record['expected_cn_by_site'][str(site_index)]:
        raise ValueError('CIF first-shell coordination differs from reviewed reference')
    metadata = dict(recipe=recipe, parser_warnings=parser_warnings,
                    cif_sha256=record['cif_sha256'], site_index=site_index,
                    fractional_coordinates=[float(v) for v in site.frac_coords],
                    site_multiplicity=int(mult), wyckoff=wyckoff,
                    site_weight=mult / sum(s[2] for s in sites),
                    space_group=cluster.struct.get_space_group_info()[0],
                    cell_composition=cluster.struct.composition.get_el_amt_dict(),
                    atom_count=len(rows), actual_cluster_radius_A=extent,
                    potential_rows=[row for _, row in card_rows(text, 'POTENTIALS')],
                    scatterer_distances_to_rpath_A=neighbors, first_shell_distances_A=first)
    return text, metadata


def run_reference(record, output, root=ROOT, recipe=None):
    from larch.xafs.feffrunner import find_exe
    executables = {}
    for module in MODULES:
        exe = find_exe(f'feff8l_{module}')
        if exe is None:
            raise RuntimeError(f'Missing FEFF8L {module}')
        executables[module] = Path(exe).resolve()
    packages = {p: importlib.metadata.version(p) for p in
                ('xraylarch', 'larixite', 'pymatgen', 'numpy', 'xraydb')}
    for index in record['absorber_site_indices']:
        target = output / record['id'] / f'site-{index}'
        target.mkdir(parents=True, exist_ok=False)
        text, metadata = prepare(record, index, root, recipe)
        (target / 'feff.inp').write_text(text)
        metadata.update(reference_id=record['id'], input_sha256=sha(target / 'feff.inp'),
                        packages=packages, executable_sha256={m: sha(e) for m, e in executables.items()},
                        engine='FEFF8L', review_status='candidate', experimental_validation=False,
                        execution_status='running', module_results=[],
                        python=sys.version.split()[0], platform=platform.platform(), machine=platform.machine())
        write_json(target / 'run.json', metadata)
        # Scratch files remain isolated; archive all raw paths and the diagnostics.
        with tempfile.TemporaryDirectory(prefix='feffdb-native-') as scratch:
            scratch = Path(scratch)
            shutil.copyfile(target / 'feff.inp', scratch / 'feff.inp')
            try:
                with (target / 'feff.log').open('wb') as log:
                    for module, exe in executables.items():
                        print(f"{record['id']} site-{index}: {module}", flush=True)
                        log.write(f'\n=== FEFF8L {module} ===\n'.encode()); log.flush()
                        started = time.monotonic()
                        result = subprocess.run([str(exe)], cwd=scratch, stdout=log, stderr=subprocess.STDOUT,
                                                timeout=1200, env={**os.environ, 'OMP_NUM_THREADS':'1', 'OPENBLAS_NUM_THREADS':'1'})
                        metadata['module_results'].append(dict(module=module, returncode=result.returncode,
                                                               elapsed_seconds=round(time.monotonic()-started, 3)))
                        write_json(target / 'run.json', metadata)
                        if result.returncode:
                            raise RuntimeError(f'{module} exit {result.returncode}; see {target}/feff.log')
                metadata['execution_status'] = 'complete'
            except BaseException as exc:
                metadata.update(execution_status='failed', error=str(exc))
                raise
            finally:
                write_json(target / 'run.json', metadata)
                archive = set(scratch.glob('feff[0-9][0-9][0-9][0-9].dat'))
                archive.update(scratch / name for name in ('files.dat', 'paths.dat', 'chi.dat', 'xmu.dat', 'misc.dat', 'log1.dat', 'log2.dat', 'log3.dat', 'log4.dat', 'log5.dat', 'log6.dat', 'pot.json', 'global.json'))
                for path in sorted(archive):
                    if path.is_file():
                        shutil.copyfile(path, target / path.name)
        report = validate_site(record, target, root)
        write_json(target / 'validation.json', report)
        print(f"{record['id']}: {len(report['selected_paths'])} first-shell paths, CN={report['coordination_sum']}", flush=True)


def match_coordination(expected, paths, tolerance=DISTANCE_TOLERANCE_A):
    """Compare the full distance multiset, preserving nearby split shells."""
    if not expected or not all(math.isfinite(v) and v > 0 for v in expected):
        raise ValueError('Invalid CIF neighbor distances')
    observed = []
    for p in paths:
        if not math.isfinite(p['reff']) or p['reff'] <= 0:
            raise ValueError('Invalid FEFF distance')
        n = float(p['degen'])
        if n <= 0 or not n.is_integer():
            raise ValueError('Nonintegral/invalid degeneracy for ordered isotropic SS path')
        observed.extend([p['reff']] * int(n))
    if len(observed) != len(expected):
        raise ValueError(f'Incomplete coordination: CIF {len(expected)}, FEFF {len(observed)}')
    differences = [abs(a-b) for a,b in zip(sorted(expected), sorted(observed))]
    if not differences or max(differences) > tolerance:
        raise ValueError(f'FEFF/CIF distance mismatch: {differences}')
    return max(differences)


def read_path(path):
    import numpy as np
    from larch.xafs.feffdat import FeffDatFile
    f = FeffDatFile(filename=str(path))
    arrays = np.array([f.k, f.real_phc, f.mag_feff, f.pha_feff, f.red_fact, f.lam, f.rep])
    if arrays.ndim != 2 or arrays.shape[0] != 7 or arrays.shape[1] < 20 or not np.isfinite(arrays).all():
        raise ValueError(f'Invalid FEFF arrays: {path}')
    if not np.all(np.diff(f.k) > 0) or f.k[0] > 0.01 or f.k[-1] < 19.99:
        raise ValueError(f'Incomplete/nonmonotonic k grid: {path}')
    if np.any(f.lam <= 0) or np.any(f.mag_feff < 0) or np.any(f.red_fact < 0):
        raise ValueError(f'Nonphysical path arrays: {path}')
    if f.nleg != 2 or len(f.geom) != 2:
        raise ValueError('Only single-scattering paths are in this pilot')
    geometry = [[str(g[0]), int(g[1]), int(g[2]), *[float(v) for v in g[3:]]] for g in f.geom]
    length = sum(math.dist(a[4:7], b[4:7]) for a, b in zip(geometry, geometry[1:]+geometry[:1])) / 2
    if not math.isfinite(length) or not math.isfinite(f.reff) or f.reff <= 0 or abs(length - f.reff) > DISTANCE_TOLERANCE_A:
        raise ValueError(f'Path geometry disagrees with reff: {path}')
    return dict(filename=path.name, sha256=sha(path), absorber=int(f.geom[0][1]),
                scatterer=int(f.geom[1][1]), nleg=int(f.nleg), edge=f.shell,
                reff=float(f.reff), degen=float(f.degen), geometry=geometry,
                kmin=float(f.k[0]), kmax=float(f.k[-1]), points=len(f.k))


def validate_site(record, directory, root=ROOT):
    from xraydb import atomic_number
    run = json.loads((directory / 'run.json').read_text())
    if (run.get('execution_status') != 'complete'
            or [(m['module'],m['returncode']) for m in run.get('module_results',[])] != [(m,0) for m in MODULES]):
        raise ValueError('FEFF run did not complete all six modules successfully')
    inp, expected = prepare(record, run['site_index'], root, run['recipe'])
    if any(run.get(key) != value for key, value in expected.items()):
        raise ValueError('Archived structural metadata differs from regenerated CIF evidence')
    if (directory / 'feff.inp').read_text() != inp or sha(directory/'feff.inp') != run['input_sha256']:
        raise ValueError('Input has changed or is not reproducible')
    paths = [read_path(p) for p in sorted(directory.glob('feff[0-9][0-9][0-9][0-9].dat'))]
    if not paths:
        raise ValueError('No FEFF path files')
    if any(p['absorber'] != atomic_number(record['absorber']) or p['edge'] != record['edge'] for p in paths):
        raise ValueError('Absorber or edge mismatch')
    selected = [p for p in paths if p['scatterer'] == atomic_number(record['scatterer']) and p['reff'] < record['shell_cutoff_A']]
    error = match_coordination(expected['first_shell_distances_A'], selected)
    log = (directory / 'feff.log').read_text(errors='replace')
    fatal = ('inconsistency in potential', 'results might be meaningless', 'ieee_invalid', 'ieee_divide_by_zero', 'ieee_overflow')
    if any(word in log.lower() for word in fatal):
        raise ValueError('Fatal FEFF diagnostic; inspect preserved log')
    diagnostics = [line.strip() for line in log.splitlines() if any(w in line.lower() for w in ('warn', 'converg', 'rmax', 'ratx', 'ieee'))]
    iterations = re.findall(r'SCF ITERATION NUMBER\s+(\d+)\s+OUT OF\s+(\d+)', log)
    if not iterations or 'Done with module 1: potentials.' not in log:
        raise ValueError('Missing native SCF execution evidence')
    if int(iterations[-1][0]) >= int(iterations[-1][1]):
        raise ValueError('SCF reached the iteration limit; review before inclusion')
    controls = json.loads((directory/'pot.json').read_text())
    recipe = run['recipe']
    if (controls['iunf'] != 0 or controls['ixc'] != recipe['exchange']
            or controls['nscmt'] != recipe['scf_iterations'] or controls['lfms1'] != 0
            or controls['nmix'] != recipe.get('scf_mixing_iterations', 1)
            or not math.isclose(controls['ca1'], recipe['scf_mixing'], abs_tol=1.e-7)
            or not math.isclose(controls['rfms1'], recipe['scf_radius_A'], abs_tol=1.e-6)
            or controls['ihole'] != {'K':1, 'L3':4}[record['edge']]):
        raise ValueError('Native potential controls differ from the recorded recipe')
    native_controls = {k:controls[k] for k in ('iunf','ixc','nscmt','lfms1','nmix','ca1','rfms1','ihole')}
    # Public Larch consumer smoke check: these are usable path inputs, not just
    # parsable headers. This is not a fit or experimental scientific validation.
    import numpy as np
    from larch.xafs import feffpath, path2chi
    for selected_path in selected:
        p = feffpath(str(directory / selected_path['filename']), s02=1, sigma2=0.003, e0=0, deltar=0)
        path2chi(p, k=np.linspace(0, 20, 401))
        if not np.isfinite(p.chi).all() or np.max(np.abs(p.chi)) <= 0:
            raise ValueError('Larch chi(k) evaluation failed')
    artifacts = {p.name:sha(p) for p in sorted(directory.iterdir()) if p.is_file() and p.name != 'validation.json'}
    previous = directory / 'validation.json'
    if previous.exists() and json.loads(previous.read_text())['artifacts'] != artifacts:
        raise ValueError('Archived calculation bytes changed; regenerate in a new directory')
    return dict(technical_validation='passed', review_status='candidate', experimental_validation=False,
                selected_paths=selected, all_paths=paths, coordination_sum=sum(p['degen'] for p in selected),
                cif_coordination=len(expected['first_shell_distances_A']), max_distance_error_A=error,
                distance_tolerance_A=DISTANCE_TOLERANCE_A, native_diagnostics=diagnostics,
                scf_iterations=int(iterations[-1][0]), scf_iteration_limit=int(iterations[-1][1]),
                scf_hit_iteration_limit=False, native_potential_controls=native_controls, larch_chi_evaluation='passed',
                artifacts=artifacts)


def _build_database(base, output, calculations, root=ROOT):
    """Append to a separate copy, keeping every upstream row/byte payload intact."""
    if output.exists() or base.resolve() == output.resolve():
        raise ValueError('Candidate output must be a new file, distinct from upstream')
    base_hash = sha(base)
    pinned_base = json.loads((root/'references.json').read_text())['upstream_base']['sha256']
    if base_hash != pinned_base or Path(str(base)+'-wal').exists() or Path(str(base)+'-shm').exists():
        raise ValueError('Use the pinned, quiescent upstream database without WAL/SHM')
    reports = [(r, index, calculations/r['id']/f'site-{index}') for r in records(root) for index in r['absorber_site_indices']]
    validations = [validate_site(r, d, root) for r, _, d in reports]
    with sqlite3.connect(base.resolve().as_uri() + '?mode=ro&immutable=1', uri=True) as src, sqlite3.connect(output) as conn:
        src.backup(conn)
        conn.execute('PRAGMA journal_mode=DELETE')
        conn.execute('PRAGMA foreign_keys=ON')
        conn.execute('CREATE TABLE IF NOT EXISTS drxas_path_review (feffdat_id INTEGER PRIMARY KEY REFERENCES feffdat(id), status TEXT NOT NULL CHECK(status="candidate"), reference_id TEXT NOT NULL, source_json TEXT NOT NULL, validation_json TEXT NOT NULL)')
        before = conn.execute('SELECT count(*) FROM feffdat').fetchone()[0]
        inserted = []
        for (r, index, directory), report in zip(reports, validations):
            ref_dir = root/'references'/r['id']
            run = json.loads((directory/'run.json').read_text())
            source = json.loads((ref_dir/'provenance.json').read_text())
            cif_id = conn.execute('INSERT INTO cif(ciftext,space_group,formula,compound,source_db,source_id) VALUES(?,?,?,?,?,?)',
                                 (lzma.compress((ref_dir/'source.cif').read_bytes()),run['space_group'],r['formula'],r['name'],*source_identity(r))).lastrowid
            pots = run['potential_rows']
            inp_id = conn.execute('INSERT INTO feffinp(inptext,absorber,edge,scatterers,natoms,cif_id) VALUES(?,?,?,?,?,?)',
                                 (lzma.compress((directory/'feff.inp').read_bytes()),int(pots[0][1]),r['edge'],json.dumps([int(p[1]) for p in pots[1:]]),run['atom_count'],cif_id)).lastrowid
            for p in report['selected_paths']:
                label = f"{r['id']}_site{index}_{p['filename'][:-4]}"
                article = f"DOI {r['doi']}" if r.get('doi') else source['primary_url']
                description = f"CANDIDATE: {r['name']}; {source_identity(r)[0]} {source_identity(r)[1]}; site {index}; {r['edge']}; first-shell {r['absorber']}-{r['scatterer']}; experimental validation pending. {article}"
                path_id = conn.execute('INSERT INTO feffdat(absorber,scatterer,nleg,reff,edge,degen,geometry,label,description,feffinp_id,feffdat) VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                                      (p['absorber'],p['scatterer'],p['nleg'],p['reff'],p['edge'],int(p['degen']),json.dumps(p['geometry']),label,description,inp_id,lzma.compress((directory/p['filename']).read_bytes()))).lastrowid
                conn.execute('INSERT INTO drxas_path_review VALUES(?,?,?,?,?)',
                             (path_id,'candidate',r['id'],json.dumps(source),json.dumps(dict(run=run,path=p,validation=report))))
                inserted.append(dict(id=path_id,reference=r['id'],filename=p['filename'],sha256=p['sha256']))
        conn.execute('INSERT INTO info(key,value) VALUES(?,?)',('drxas_element_round2',json.dumps(dict(status='candidate',base_sha256=base_hash,added=len(inserted)))))
        if conn.execute('PRAGMA integrity_check').fetchone()[0] != 'ok' or conn.execute('PRAGMA foreign_key_check').fetchall():
            raise ValueError('Candidate database integrity check failed')
        # Every pre-existing table and row must compare identically, including BLOBs.
        for (table,) in src.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"):
            if not re.fullmatch('[a-z_]+', table):
                raise ValueError('Unexpected base table name')
            original = src.execute(f'SELECT * FROM "{table}"').fetchall()
            updated = set(conn.execute(f'SELECT * FROM "{table}"').fetchall())
            if not all(row in updated for row in original):
                raise ValueError(f'Upstream rows changed: {table}')
        conn.commit()
    if sha(base) != base_hash:
        raise ValueError('Upstream database changed')
    result = dict(base_sha256=base_hash, candidate_sha256=sha(output), upstream_paths=before,
                  added_paths=len(inserted), total_paths=before+len(inserted), review_status='candidate', paths=inserted)
    write_json(output.with_suffix('.manifest.json'), result)
    return result


def build_database(base, output, calculations, root=ROOT):
    if output.exists() or output.with_suffix('.manifest.json').exists() or output.resolve() == base.resolve():
        raise ValueError('Candidate output and manifest must be new files')
    with tempfile.TemporaryDirectory(prefix='.feffdb-build-', dir=output.parent) as work:
        temp = Path(work)/output.name
        result = _build_database(base, temp, calculations, root)
        temp.replace(output)
        temp.with_suffix('.manifest.json').replace(output.with_suffix('.manifest.json'))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=('run','validate','build-db'))
    ap.add_argument('--root', type=Path, default=ROOT)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--reference')
    ap.add_argument('--base', type=Path)
    ap.add_argument('--calculations', type=Path)
    args = ap.parse_args()
    calcs = args.calculations or args.root/'calculations'
    if args.command == 'build-db':
        if args.base is None or args.output is None:
            ap.error('build-db requires --base and --output')
        print(json.dumps(build_database(args.base,args.output,calcs,args.root), indent=2))
        return
    chosen = [r for r in records(args.root) if args.reference is None or r['id'] == args.reference]
    if not chosen:
        ap.error('Unknown reference')
    for r in chosen:
        if args.command == 'run':
            run_reference(r, args.output or calcs, args.root)
        else:
            for index in r['absorber_site_indices']:
                result = validate_site(r, calcs/r['id']/f'site-{index}', args.root)
                print(json.dumps(dict(reference=r['id'],site=index,paths=len(result['selected_paths']),coordination=result['coordination_sum'],status=result['technical_validation'])))


if __name__ == '__main__':
    main()
