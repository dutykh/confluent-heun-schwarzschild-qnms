"""Regenerate the 68 finite-matrix frequencies, including the N=640 run.

Default: check the three archived raw grid records, finite-depth continued
fractions and contour comparisons against the paper's historical summary.
--recompute independently seeds from the numerical spectrum at N=40.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import time
import mpmath as mp

HERE = Path(__file__).resolve().parent
BASE_GRIDS = (40,60,80,100,160,240,320)


def read(stem):
    return json.loads((HERE/(stem+'.json')).read_text())


def mode(row):
    return row['ell'],row['provisional_overtone']


def by_grid(rows):
    return {(r['ell'],r['provisional_overtone'],r['N']):r for r in rows}


def grid_schedule():
    for ell in range(3):
        for grid in BASE_GRIDS:
            for overtone in range(3):
                yield ell,overtone,grid
    for grid in (400,480):
        for overtone in (1,2):
            yield 0,overtone,grid
    yield 0,2,640


def track(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[mode(r)].append(r)
    latest = []
    with mp.workdps(90):
        for key in sorted(groups):
            prev = None
            for row in sorted(groups[key],key=lambda r:r['N']):
                now = mp.mpc(*row['omega_strings'])
                row['previous_grid_difference'] = None if prev is None else mp.nstr(abs(now-prev),20)
                prev = now
            latest.append(max(groups[key],key=lambda r:r['N']))
    return latest


def without_timings(row):
    return {k:v for k,v in row.items() if k != 'elapsed_seconds'}


def check_saved(spectral_path=None):
    raw = []
    for name in ('highprecision','highprecision_extra','highprecision_640'):
        raw.extend(read('qnm_spectral_'+name)['results'])
    archived = (json.loads(spectral_path.read_text()) if spectral_path
                else read('qnm_spectral_results'))
    records = by_grid(archived['results'])
    if (len(raw) != len(records) or len(archived['results']) != len(records)
            or set(by_grid(raw)) != set(grid_schedule())
            or set(records) != set(by_grid(raw))):
        raise ValueError('Wrong spectral grid or duplicate record')
    actual = track(raw)
    for key,r in by_grid(raw).items():
        if without_timings(r) != without_timings(records[key]):
            raise ValueError('Grid record differs: '+str(key))
    if list(map(without_timings,actual)) != list(map(without_timings,archived['latest_by_branch'])):
        raise ValueError('Released terminal spectral rows do not match grid refinements')
    for key,sidecar in [('precision_checks','precision_check'),('residual_diagnostics','diagnostics')]:
        if archived[key] != read('qnm_spectral_'+sidecar)['results']:
            raise ValueError('Spectral diagnostic sidecar differs: '+key)
    cf = {(r['ell'],r['overtone']):r for r in read('qnm_leaver_results')['modes']}
    double = {(r['ell'],r['n']):r for r in read('qnm_contour_results')['modes']}
    multi = {(r['ell'],r['n']):r for r in read('qnm_contour_mp_results')['modes']}
    saved = {(r['ell'],r['n']):r for r in read('qnm_comparison_results')}
    if set(cf) != set(saved) or set(cf) != set(double) or len(saved) != 9:
        raise ValueError('Mode coverage differs across numerical methods')
    with mp.workdps(70):
        for r in actual:
            key = mode(r)
            record = saved[key]
            f = cf[key]['omega']
            ref = mp.mpc(f['real'],f['imag'])
            base = mp.mpc(*double[key]['runs'][-2]['omega'])
            h = multi[key]['root'] if key in multi else None
            contour = mp.mpc(h['real'],h['imag']) if h else base
            spec = mp.mpc(*r['omega_strings'])
            if (record['spectral_degree'] != r['N'] or
                    record['spectral_last_grid_change'] != r['previous_grid_difference'] or
                    record['contour_multiprecision'] != bool(h) or
                    mp.mpc(record['contour']['real'],record['contour']['imag']) != contour or
                    str(abs(base-ref)) != record['double_contour_cf_difference'] or
                    str(abs(contour-ref)) != record['contour_cf_difference'] or
                    str(abs(spec-ref)) != record['spectral_cf_difference']):
                raise ValueError(f'Cross-method comparison differs at {key}')
    print('PASS: 68 spectral grid roots, nine terminal grid changes and nine historical '
          'CF/spectral/contour comparisons agree with their raw records.')


def recompute(output):
    from qnm_spectral import metadata
    from qnm_spectral_chebcoef import solve
    from qnm_spectral_refine import refine, collocation_diagnostics
    seeds = {ell:solve(ell,BASE_GRIDS[0])[:3].copy() for ell in range(3)}
    data = metadata()
    started = time.monotonic()
    rows = []
    for ell,n,grid in grid_schedule():
        t0 = time.monotonic()
        r = refine(ell,grid,seeds[ell][n],50)
        seeds[ell][n] = complex(*map(float,r['omega_strings']))
        r.update(provisional_overtone=n,elapsed_seconds=time.monotonic()-t0)
        rows.append(r)
        print(ell,n,grid,r['omega_strings'],flush=True)
    latest = track(rows)
    registry = by_grid(rows)
    repeats = []
    for ell,grid,n in ((0,160,2),(2,100,2)):
        reference = registry[ell,n,grid]['omega_strings']
        r = refine(ell,grid,complex(*map(float,reference)),70,tol_digits=50)
        r.update(provisional_overtone=n,difference_from_50_digit_run_printed_values=
                 mp.nstr(abs(mp.mpc(*r['omega_strings'])-mp.mpc(*reference)),20))
        repeats.append(r)
    diagnostics = []
    for ell in range(3):
        for n in range(3):
            target = registry[ell,n,320]['omega_strings']
            r = collocation_diagnostics(ell,320,target,60)
            r['provisional_overtone'] = n
            diagnostics.append(r)
    data.update(grids=list(BASE_GRIDS),extra_grids={'ell0 n1':[400,480],'ell0 n2':[400,480,640]},
                results=rows,latest_by_branch=latest,precision_checks=repeats,
                residual_diagnostics=diagnostics)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data,indent=2)+'\n')
    print('Regenerated 68 finite matrix roots in',round(time.monotonic()-started,2),'seconds.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check-result',type=Path,metavar='SPECTRAL_JSON',
                       help='compare regenerated grid histories and diagnostics with the archive, excluding times')
    group.add_argument('--recompute',type=Path,metavar='OUTPUT_JSON',
                        help='run all finite-matrix eigensolves and write this file')
    args = parser.parse_args()
    if args.recompute:
        recompute(args.recompute)
    else:
        check_saved(args.check_result)
