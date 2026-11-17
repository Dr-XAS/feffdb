#!/usr/bin/env python3
"""Bounded oxide FEFF sensitivity, never a promotion or experimental validation.

Run: python sensitivity.py run --workers 3
Resume: python sensitivity.py run --workers 3 (completed archives revalidated)
Reanalyse without native FEFF: python sensitivity.py analyse
All output stays beside this script; source_snapshot is an unchanged pilot copy.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import traceback

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/feffdb-round1-mpl')
ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'source_snapshot'
SOURCE_COMMIT = '547463406c812c2ea70333cb3861d80be8ceeb6c'
SETTINGS = {
    'baseline_scf5_cluster8': (5., 8.),
    'scf6_cluster8': (6., 8.),
    'scf5_cluster9': (5., 9.),
    'scf4_cluster8': (4., 8.),
    'scf5_cluster7': (5., 7.),
}
WINDOWS = {'3-8': (3., 8.), '8-12': (8., 12.), '12-16': (12., 16.),
           '16-20': (16., 20.), '3-12': (3., 12.)}
THRESHOLDS = {'complex_relative_rms': .05, 'amplitude_relative_rms': .03,
              'phase_weighted_rms_rad': .03}
DOCS = [
    'https://feff.phys.washington.edu/FEFF/Docs/feff8/feff84/feff84/SCF.html',
    'https://feff.phys.washington.edu/feff/Docs/feff8/feff8web/node49.html',
    'https://xraypy.github.io/xraylarch/xafs_feffpaths.html',
]


def load_pilot():
    spec = importlib.util.spec_from_file_location('oxide_pilot_snapshot', SOURCE / 'pilot.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    tmp.replace(path)


def verify_source():
    inventory = json.loads((SOURCE / 'artifact-sha256.json').read_text())
    bad = [name for name, digest in inventory.items() if sha(SOURCE / name) != digest]
    if bad:
        raise ValueError(f'Pinned source changed: {bad}')


def protocol():
    return dict(source_repository='https://github.com/Dr-XAS/feffdb', source_commit=SOURCE_COMMIT,
                source_inventory_sha256=sha(SOURCE / 'artifact-sha256.json'),
                source_script_sha256=sha(SOURCE / 'pilot.py'),
                settings={k:dict(scf_radius_A=s, cluster_radius_A=c) for k,(s,c) in SETTINGS.items()},
                windows_invA=WINDOWS, analysis_parameters=dict(s02=1, sigma2_A2=0, e0_eV=0,
                deltar_A=0, native_degeneracy=True, kstep_invA=.05, weight=2),
                review_thresholds=THRESHOLDS,
                threshold_reason='Predeclared numerical review triggers, not literature standards or experimental uncertainty. '
                    '5% complex / 3% amplitude changes flag material numerical dependence. A 0.03 rad phase change '
                    'would correspond roughly to 0.0025 A at k=6 if entirely represented as 2k*dR; this is not a fit bias estimate.',
                phase_mask='Keep points where baseline |complex path| >= 5% of its maximum in the window; '
                    'use circular phase differences weighted by squared k^2 baseline amplitude.',
                complex_convention='Larch complex path reconstructed as -chi_imag + 1j*chi. '
                    'Full shell is the complex sum with native degeneracies; sigma2=0 avoids damping high-k differences.',
                interpretation='Separate one-factor variations. Lower-radius probes are sensitivity bounds; '
                    'larger-radius comparisons test local numerical stability only, not global convergence.',
                experimental_validation=False, promote_candidates=False, status='candidate', primary_docs=DOCS)


def worker(reference_id, setting):
    pilot = load_pilot()
    record = next(r for r in pilot.records(SOURCE) if r['id'] == reference_id)
    scf, cluster = SETTINGS[setting]
    recipe = dict(pilot.RECIPE, scf_radius_A=scf, cluster_radius_A=cluster)
    output = ROOT / 'calculations' / setting
    target = output / reference_id / 'site-1'
    try:
        if target.exists():
            pilot.validate_site(record, target, SOURCE)
            run = json.loads((target / 'run.json').read_text())
            if run['recipe'] != recipe:
                raise ValueError('Existing run has different recipe')
        else:
            pilot.run_reference(record, output, SOURCE, recipe=recipe)
        return dict(reference=reference_id, setting=setting, status='complete')
    except Exception as exc:
        return dict(reference=reference_id, setting=setting, status='failed', error=str(exc),
                    traceback=traceback.format_exc())


def metrics(k, baseline, varied, limits):
    import numpy as np
    mask = (k >= limits[0]) & (k <= limits[1])
    weights = k[mask] ** 2
    a, b = baseline[mask], varied[mask]
    wa, wb = weights * a, weights * b
    denominator = float(np.linalg.norm(wa))
    if denominator <= 0:
        raise ValueError('Zero reference signal')
    phase_mask = np.abs(a) >= .05 * np.max(np.abs(a))
    delta = np.angle(b[phase_mask] * np.conj(a[phase_mask]))
    phase_weights = np.abs(wa[phase_mask]) ** 2
    phase_rms = float(np.sqrt(np.average(delta ** 2, weights=phase_weights)))
    chi_denominator = float(np.linalg.norm(wa.imag))
    result = dict(points=int(mask.sum()), phase_points=int(phase_mask.sum()),
                  complex_relative_rms=float(np.linalg.norm(wb-wa)/denominator),
                  amplitude_relative_rms=float(np.linalg.norm(np.abs(wb)-np.abs(wa))/denominator),
                  phase_weighted_rms_rad=phase_rms, phase_max_abs_rad=float(np.max(np.abs(delta))),
                  chi_relative_rms=float(np.linalg.norm(wb.imag-wa.imag)/chi_denominator),
                  baseline_weighted_rms=float(np.sqrt(np.mean(np.abs(wa)**2))))
    result['review_flags'] = [name for name,limit in THRESHOLDS.items() if result[name] > limit]
    return result


def read_complex(directory, selected, k):
    import numpy as np
    from larch.xafs import feffpath, path2chi
    result = {}
    for record in selected:
        key = f"O_R{record['reff']:.4f}_N{record['degen']:g}"
        if key in result:
            raise ValueError('Ambiguous split-path identity')
        path = feffpath(str(directory / record['filename']), degen=record['degen'],
                        s02=1, sigma2=0, e0=0, deltar=0)
        path2chi(path, k=k)
        z = -np.asarray(path.chi_imag) + 1j*np.asarray(path.chi)
        if not np.isfinite(z).all():
            raise ValueError('Nonfinite complex signal')
        result[key] = z
    result['complete_O_shell'] = sum(result.values())
    return result


def analyse():
    import numpy as np
    verify_source()
    pilot = load_pilot()
    k = np.linspace(0, 20, 401)
    report = dict(protocol=protocol(), results=[], baseline_reproduction=[], failures=[],
                  scientific_status='candidate', experimental_validation=False, auto_promoted=False)
    curve_root = ROOT / 'curves'; curve_root.mkdir(exist_ok=True)
    rows = []
    for reference in pilot.records(SOURCE):
        rid = reference['id']; baseline_dir = ROOT/'calculations'/'baseline_scf5_cluster8'/rid/'site-1'
        if not (baseline_dir/'validation.json').is_file():
            report['failures'].append(dict(reference=rid, setting='baseline_scf5_cluster8', reason='No validated run'))
            continue
        bv = pilot.validate_site(reference, baseline_dir, SOURCE)
        baseline = read_complex(baseline_dir, bv['selected_paths'], k)
        original = SOURCE/'calculations'/rid/'site-1'
        original_run = json.loads((original/'run.json').read_text())
        br = json.loads((baseline_dir/'run.json').read_text())
        reproduce_names = ['feff.inp']+[x['filename'] for x in bv['all_paths']]
        report['baseline_reproduction'].append(dict(reference=rid,
            compared_files=len(reproduce_names), all_input_and_raw_path_bytes_identical=all(
                sha(baseline_dir/name)==sha(original/name) for name in reproduce_names),
            executable_fingerprints_match_source=br['executable_sha256']==original_run['executable_sha256'],
            package_versions_match_source=br['packages']==original_run['packages'],
            differences=[name for name in reproduce_names if sha(baseline_dir/name)!=sha(original/name)]))
        for setting in SETTINGS:
            directory = ROOT/'calculations'/setting/rid/'site-1'
            try:
                v = pilot.validate_site(reference,directory,SOURCE)
                run = json.loads((directory/'run.json').read_text())
                varied = read_complex(directory,v['selected_paths'],k)
                if varied.keys()!=baseline.keys():raise ValueError('Path identity/split shell changed')
                if run['executable_sha256']!=br['executable_sha256']:raise ValueError('Executable changed')
            except Exception as exc:
                report['failures'].append(dict(reference=rid,setting=setting,reason=str(exc)));continue
            result = dict(reference=rid, setting=setting, recipe=run['recipe'], atom_count=run['atom_count'],
                          actual_cluster_radius_A=run['actual_cluster_radius_A'], scf_iterations=v['scf_iterations'],
                          native_diagnostics=v['native_diagnostics'], selected_path_count=len(v['selected_paths']),
                          coordination_sum=v['coordination_sum'], paths={})
            scf_sizes = {}
            for count, potential in re.findall(r'Doing FMS for a cluster of\s+(\d+) atoms around iph\s*=\s*(\d+)', (directory/'feff.log').read_text()):
                scf_sizes.setdefault(potential, set()).add(int(count))
            result['scf_cluster_atom_counts_by_potential'] = {key:sorted(value) for key,value in scf_sizes.items()}
            for key,z in varied.items():
                result['paths'][key] = {name:metrics(k,baseline[key],z,limits) for name,limits in WINDOWS.items()}
                matrix=np.column_stack((k, baseline[key].real, baseline[key].imag, z.real,z.imag,
                    np.abs(baseline[key]),np.abs(z),np.angle(z*np.conj(baseline[key]))))
                np.savetxt(curve_root/f'{rid}__{setting}__{key}.csv',matrix,delimiter=',',
                    header='k,baseline_real,baseline_chi,varied_real,varied_chi,baseline_amplitude,varied_amplitude,circular_phase_difference_rad')
                for window,m in result['paths'][key].items():
                    rows.append(dict(reference=rid, setting=setting, path=key, window=window,
                                     **{a:b for a,b in m.items() if a!='review_flags'},
                                     review_flags=';'.join(m['review_flags'])))
            report['results'].append(result)
    write_json(ROOT/'results.json',report)
    if rows:
        with (ROOT/'metrics.csv').open('w',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    write_summary(report)
    inventory={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*'))
        if p.is_file() and '__pycache__' not in p.parts and p not in
        (ROOT/'artifact-sha256.json', ROOT/'run-progress.json', ROOT/'driver.log')}
    write_json(ROOT/'artifact-sha256.json',inventory)
    return report


def write_summary(report):
    lines=['# Oxide pilot radius sensitivity', '',
           'Numerical sensitivity only. All three references remain candidate; no experimental validation or promotion.', '',
           'Baseline SCF=5 A, cluster=8 A. Vary one control at a time. RPATH=4 A and all other native settings retained.',
           'Signals use native degeneracy, S02=1, sigma2=0, E0=0, dR=0. Metrics use k^2 weighting.', '',
           'Review triggers (empirical, predeclared): complex NRMS >5%, amplitude NRMS >3%, amplitude-weighted circular phase RMS >0.03 rad.',
           'These are numerical review thresholds, not accuracy guarantees. Lower-radius probes are not convergence tests.', '',
           '| Reference | Setting | k window | Max per-path complex NRMS | Max amplitude NRMS | Max phase RMS (rad) | Shell complex NRMS |',
           '|---|---|---|---:|---:|---:|---:|']
    for result in report['results']:
        if result['setting'].startswith('baseline'):continue
        for window in WINDOWS:
            paths=[w[window] for key,w in result['paths'].items() if key!='complete_O_shell']
            shell=result['paths']['complete_O_shell'][window]
            lines.append(f"| {result['reference']} | {result['setting']} | {window} | "
                f"{max(x['complex_relative_rms'] for x in paths):.3%} | "
                f"{max(x['amplitude_relative_rms'] for x in paths):.3%} | "
                f"{max(x['phase_weighted_rms_rad'] for x in paths):.6f} | {shell['complex_relative_rms']:.3%} |")
    lines+=['','## Baseline repeatability','',json.dumps(report['baseline_reproduction'],indent=2),'',
            '## Failures','',json.dumps(report['failures'],indent=2),'',
            'The native radius-counter warning and underflow/denormal flags are preserved. Successful modules and stopping before the SCF iteration cap do not prove radius convergence.',
            '', 'Primary methodological references:', *['- '+url for url in DOCS], '']
    (ROOT/'README.md').write_text('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['run','analyse'])
    parser.add_argument('--workers',type=int,default=3);args=parser.parse_args()
    verify_source()
    if args.action=='run':
        write_json(ROOT/'protocol.json',protocol())
        records=load_pilot().records(SOURCE)
        tasks=[(r['id'],setting) for setting in SETTINGS for r in records]
        results=[]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures={pool.submit(worker,*task):task for task in tasks}
            for future in as_completed(futures):
                result=future.result();results.append(result);write_json(ROOT/'run-progress.json',results)
                print(json.dumps(result),flush=True)
    report=analyse()
    print(json.dumps(dict(completed=len(report['results']),failures=report['failures'])),flush=True)
    if report['failures']:sys.exit(1)


if __name__=='__main__':main()
