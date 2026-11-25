"""Verify the stored release without writing to it. --full also rebuilds the
ODE certificates and regressions in a temporary directory and compares them.
Giving an explicit --workers count selects that full replay for compatibility.
"""
import argparse
import hashlib
from fractions import Fraction
import json
from pathlib import Path
import os
import shutil
import tempfile
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def load_result(stem):
    return json.loads((HERE / f'results_{stem}.json').read_text())


def require(condition, description):
    if not condition:
        raise ArithmeticError(description)


def check_numerical_manifest():
    root = HERE.parent
    manifest_path = root / 'NUMERICAL_SHA256.txt'
    manifest = manifest_path.read_bytes()
    checksum = (root / 'NUMERICAL_SHA256.sha256').read_text(encoding='ascii')
    expected = f'{hashlib.sha256(manifest).hexdigest()}  {manifest_path.name}'
    require(checksum.splitlines() == [expected],
            'Numerical manifest differs from its detached checksum')
    entries = manifest.decode('utf-8').splitlines()
    for entry in entries:
        digest, relative_path = entry.split('  ', 1)
        path = root / relative_path
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                f'Numerical release file differs: {relative_path}')
    print(f'PASS: numerical manifest checksum and {len(entries)} numerical release files.', flush=True)


def check_exceptional_results():
    F = Fraction
    threshold = load_result('threshold')
    require(threshold['ell'] == 0 and F(threshold['frequency_radius']) == F(1, 1024),
            'Threshold normalization differs')
    require({case['sheet'] for case in threshold['cases']} == {'upper_spatial', 'lower_spatial'},
            'Expected both threshold continuations')
    for case in threshold['cases']:
        lo, hi = map(F, case['wronskian']['re'])
        require(F('-1.386') < lo <= hi < F('-0.621'), 'Threshold Wronskian bound failed')
    tests = load_result('lateral')['tested']
    wanted = {(imag, width) for imag in ('-0.25', '-0.5') for width in ('0', '1e-12')}
    require({(t['frequency_center_imag'], t['box_halfwidth']) for t in tests} == wanted,
            'Wrong set of axis cases')
    for test in tests:
        imag = test['frequency_center_imag']
        lower = F(test['min_lateral_modulus'][0])
        require(lower > F('3.22' if imag == '-0.25' else '1.04'), 'Lateral modulus bound failed')
        if imag == '-0.25' and test['box_halfwidth'] == '0':
            lo, hi = map(F, test['jump']['imag'])
            require(F('-5.962910178215331') < lo <= hi < F('-5.962910178215330'),
                    'Cut jump bound failed')
    frame = load_result('resonant_horizon')
    require(F(frame['chart_half_side_in_lambda']) == F(1, 1024), 'Wrong horizon box radius')
    require({(c['ell'], c['N']) for c in frame['cases']} ==
            {(ell, n) for ell in range(3) for n in (1, 2, 3)}, 'Unexpected horizon boxes')
    for case in frame['cases']:
        require(F(case['h_abs_lower']) > F('0.9775'), 'Horizon nonvanishing bound failed')
    accessory = load_result('resonant')
    require({F(c['epsilon']) for c in accessory['brackets']} ==
            {F(0), F('1e-7'), F('1e-5'), F('1e-3')}, 'Unexpected accessory parameters')
    for case in accessory['brackets']:
        mid, radius = F(case['midpoint_decimal']), F(case['radius'])
        require(radius == F('1e-28') and F(case['t_lower']) == mid-radius
                and F(case['t_upper']) == mid+radius, 'Accessory bracket endpoints differ')
        require(F(case['left']['complete_enclosure'][0]) > 0
                and F(case['right']['complete_enclosure'][1]) < 0, 'Accessory signs unresolved')
    print('Threshold, lateral, resonant horizon and accessory inequalities passed.', flush=True)


def invoke(arguments, cwd=HERE.parent):
    env = {**os.environ, 'PYTHONDONTWRITEBYTECODE':'1'}
    # The historical regression programs include assertions. Ensure subprocesses
    # retain these checks even if the caller exported PYTHONOPTIMIZE.
    env.pop('PYTHONOPTIMIZE', None)
    subprocess.run([sys.executable, '-B', *arguments], cwd=cwd, env=env, check=True)


def without_times(item):
    if isinstance(item, dict):
        return {key:without_times(value) for key,value in item.items()
                if key not in {'elapsed_seconds','seconds','refined_seconds'}}
    if isinstance(item, list):
        return [without_times(value) for value in item]
    return item


def full_replay(workers):
    with tempfile.TemporaryDirectory(prefix='che-certificate-') as tmp:
        fresh = Path(tmp)
        dest = fresh/'cert'
        dest.mkdir()
        (dest/'results').mkdir()
        for source in HERE.glob('*.py'):
            shutil.copy2(source, dest/source.name)
        invoke(['-m','cert.elementary'], fresh)
        for module in ('cert.reproduce','cert.neighborhoods'):
            invoke(['-m',module,'--workers',str(workers)], fresh)
        for module in ('cert.lateral','cert.threshold','cert.resonant_horizon','cert.resonant'):
            invoke(['-m',module], fresh)
        invoke(['cert/generic.py'], fresh)
        for name in ('test_transport','test_quadrature','test_isolate',
                     'test_neighborhoods','test_lateral'):
            invoke(['-m','cert.'+name], fresh)
        expected = {p.relative_to(HERE) for p in HERE.rglob('*.json')}
        regenerated = {p.relative_to(dest) for p in dest.rglob('*.json')}
        require(regenerated == expected,
                f'Certificate output file set differs: {regenerated ^ expected}')
        for path in sorted(expected):
            a = json.loads((HERE/path).read_text())
            b = json.loads((dest/path).read_text())
            require(without_times(a)==without_times(b),f'Recomputed certificate differs: {path}')
        print(f'PASS: {len(expected)} complete certificate JSON files reproduce '
              'with exact agreement after wall times are removed.', flush=True)
        comparison = fresh/'schultze.json'
        invoke(['research/check_schultze.py','--output',str(comparison)])
        saved = HERE.parent/'research'/'schultze_comparison.json'
        require(json.loads(comparison.read_text()) == json.loads(saved.read_text()),
                'Full Schultze recurrence replay differs')
        print('PASS: full Schultze recurrence replay reproduces the saved record.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--full', '--recompute', action='store_true')
    group.add_argument('--check-saved', action='store_true',
                       help='check saved data and tables (default)')
    parser.add_argument('--workers', type=int, default=None,
                        help='full replay worker count (default one for --full)')
    args = parser.parse_args()
    if args.workers is not None and (args.workers < 1 or args.check_saved):
        parser.error('--workers must be positive and selects full recomputation')
    check_numerical_manifest()
    for module in ('cert.reproduce','cert.neighborhoods'):
        invoke(['-m',module,'--check-saved'])
    check_exceptional_results()
    manuscript = HERE.parent / 'reports' / 'Confluent_Heun_Generic_Progress.tex'
    if manuscript.exists():
        invoke(['-m','cert.check_tables'])
    else:
        print('SKIP: printed tables, since no manuscript source is present in reports/.', flush=True)
    invoke(['analysis/qnm_spectral_release.py'])
    invoke(['research/check_accessory_strip.py'])
    invoke(['research/check_schultze_phase.py'])
    invoke(['research/check_schultze.py','--check-saved'])
    if args.full or args.workers is not None:
        full_replay(args.workers or 1)


if __name__ == '__main__':
    main()
