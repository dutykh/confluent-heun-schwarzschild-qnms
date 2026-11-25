"""Compare independently computed frequency records; no ground-truth claim."""
import argparse
import json
from pathlib import Path
import mpmath as mp

BASE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--spectral',type=Path, default=BASE/'qnm_spectral_results.json',
                    help='complete spectral grid/terminal-frequency record')
parser.add_argument('--output-dir',type=Path, default=BASE,
                    help='directory for the comparison JSON and historical TeX fragment')
args=parser.parse_args()
mp.mp.dps=70
def read(name):return json.loads((BASE/name).read_text())
def z(d):return mp.mpc(d['real'],d['imag'])
def sci(t):
    s=mp.nstr(t,3,min_fixed=0,max_fixed=0)
    if 'e' in s:
        a,b=s.split('e');return a+r'\times10^{'+str(int(b))+'}'
    return s

cf={(d['ell'],d['overtone']):d for d in read('qnm_leaver_results.json')['modes']}
co={(d['ell'],d['n']):d for d in read('qnm_contour_results.json')['modes']}
hi={(d['ell'],d['n']):d for d in read('qnm_contour_mp_results.json')['modes']}
spectral=json.loads(args.spectral.read_text())
if 'latest_by_branch' not in spectral:
    parser.error('This is a partial spectral grid. Generate the full record with qnm_spectral_release.py --recompute SPECTRAL_JSON.')
sp={(d['ell'],d['provisional_overtone']):d for d in spectral['latest_by_branch']}
degrees={(ell,n):320 for ell in range(3) for n in range(3)}
degrees.update({(0,1):480,(0,2):640})
if len(spectral['latest_by_branch']) != 9 or {k:r['N'] for k,r in sp.items()} != degrees:
    parser.error('The historical comparison requires all nine branches, including the monopole grids at degrees 480 and 640.')
rows=[]; records=[]
for key in sorted(cf):
    ref=z(cf[key]['omega'])
    baseline=mp.mpc(*co[key]['runs'][-2]['omega'])
    c=z(hi[key]['root']) if key in hi else baseline
    s=mp.mpc(*sp[key]['omega_strings'])
    dc=abs(c-ref);ds=abs(s-ref)
    label=str(key[0])+(r'^{\dagger}' if key in hi else '')
    rows.append(f"${label}$ & {key[1]} & ${float(c.real):.12f}$ & ${float(-c.imag):.12f}$ & ${sci(dc)}$ & ${sci(ds)}$ "+r'\\')
    records.append(dict(ell=key[0],n=key[1],contour={'real':str(c.real),'imag':str(c.imag)},
        contour_multiprecision=key in hi,contour_cf_difference=str(dc),spectral_cf_difference=str(ds),
        double_contour_cf_difference=str(abs(baseline-ref)),spectral_degree=sp[key]['N'],
        spectral_last_grid_change=sp[key]['previous_grid_difference']))
args.output_dir.mkdir(parents=True, exist_ok=True)
(args.output_dir/'qnm_comparison_results.json').write_text(json.dumps(records,indent=2))
section=r'''\subsection{Results and observed limitations}
The table reports $\Omega=2M\omega$. Differences are absolute complex-frequency differences from the independently computed finite-depth continued-fraction value. They are comparisons, not certified errors. A dagger marks the three arbitrary-precision contour calculations; other contour entries use the final double-precision run. The frequency columns are rounded to twelve decimal places.
\begin{center}\small
\begin{tabular}{rrrrrr}\toprule
$\ell$ & $n$ & $\operatorname{Re}\Omega$ & $-\operatorname{Im}\Omega$ & Contour--CF & Spectral--CF\\\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule\end{tabular}\end{center}

The spectral degrees are $320$ except for $(\ell,n)=(0,1)$ and $(0,2)$, which use $480$ and $640$. For $(0,2)$ its last grid change is $4.62\times10^{-9}$, despite the smaller cross-method difference; neither number is a proved error bound. A small residual cannot resolve this discrepancy: an earlier degree-$320$ polynomial had a continuum ODE residual of order $10^{-44}$ while its frequency still moved by order $10^{-7}$ under refinement. The endpoint-continuation and non-self-adjoint conditioning issues are substantive.

The double-precision contour $(0,2)$ root differed from the continued fraction by about $7.0\times10^{-11}$ and changed by $2.0\times10^{-10}$ when the matching radius was varied. Its infinity-tail-only determinant bound was already of order $10^{-52}$, so that bound could not explain the observed numerical error. Repeating the calculation with a 40-digit Taylor propagator, then 55-digit checks, reduced the cross-method difference to $1.48\times10^{-28}$. The continued-fraction reference itself has a double-depth Newton correction estimate of $1.51\times10^{-28}$; the agreement is therefore reference-limited at that scale. It is not a 28-digit certification.

The selected arbitrary-precision contour roots use $R=140$, inverse-power order $65$, Taylor order $48$, and step fraction $0.20$. Root searches began from rounded seeds. For $(0,2)$, changing the matching radius $5\to5.4$, Taylor order $48\to56$, and step fraction $0.20\to0.18$ at 55 digits gives $|A_{\rm in}|=3.57\times10^{-29}$ at the same root. This supplies an additional contour/propagation check. Separate function and derivative enclosures for the numerical Taylor transport are still required for certification.

As a normalization check independent of root location, real-frequency scattering at $\Omega=0.5,1$ for $\ell=0,1,2$ satisfies $|A_{\rm in}|^2-|A_{\rm out}|^2=1$ with absolute discrepancies between $4.4\times10^{-16}$ and $3.7\times10^{-12}$ in the double-precision calculation. The largest discrepancy occurs when subtracting large scattering amplitudes. The exceptional-horizon evaluator also passed 53 numerical checks, including 36 exact resonance cases through $j=32$, generic and near-resonant comparisons, and static Legendre checks. The largest exceptional-identity discrepancy was $1.95\times10^{-93}$ at elevated working precision; these are consistency checks, not interval proofs.

The archived JSON records retain precision, truncation, grid, contour, residual, and comparison data. No speed advantage is claimed: elapsed times from different arithmetic settings and implementations are not an equitable performance benchmark.
'''
(args.output_dir/'qnm_numerical_section.tex').write_text(section)
print(json.dumps(records,indent=2))
