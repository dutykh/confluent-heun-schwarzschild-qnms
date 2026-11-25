#!/usr/bin/env python3
"""Independent scalar Schwarzschild QNM Chebyshev-collocation benchmark.

Run, from the containing project directory:
  OPENBLAS_NUM_THREADS=1 python3 analysis/qnm_spectral.py --output /tmp/che-spectral.json
With no grid or precision overrides this regenerates the complete release
benchmark through degree 640. Specifying an override requests a custom grid.
Requires numpy, scipy, sympy, mpmath. Companion local modules:
qnm_spectral_chebcoef.py and qnm_spectral_refine.py.

Omega=2M*omega, exp(-i omega t), z=r/(2M), y=(z-1)/z.
u=e^(i Omega z)*(z-1)^(-i Omega)*z^(2 i Omega)*v(y).
Writing lambda=i Omega, h=y(1-y)^2, k=1-2(1-y)^2 gives
 h v''+(h'+2 lambda k)v'
 +[4 lambda(1-y)-4 lambda^2(2-y)-ell(ell+1)-1+y]v=0.
Endpoint rows are the equation's limiting Robin relations.

The default uses exactly assembled Chebyshev-coefficient matrices, equivalent
to Lobatto collocation: the residual has degree <=N+1 and T_(N+1)=T_(N-1)
on those nodes. This avoids severe nodal D^2 conditioning. Arbitrary-precision
Newton refines the finite-dimensional roots, not the continuum eigenproblem.
No Leaver frequencies or externally tabulated values are used as input.

CAUTION: in damped frequencies the incoming relative solution at infinity is
exp[-2i Omega/(1-y)]*(1-y)^(2i Omega), exponentially flat on this real interval.
Smooth endpoint behavior alone does not prove the outgoing analytic continuation.
Only independently checked stable branches should be interpreted as QNMs.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import sympy as sp
from qnm_spectral_chebcoef import coefficient_matrices,solve
from qnm_spectral_refine import refine, logdet_derivative


def symbolic_checks():
    y,O,L=sp.symbols('y O L')
    h=y*(1-y)**2
    phase=1/(1-y)-sp.log(y)-sp.log(1-y)
    k=1-2*(1-y)**2
    V=y*(L*(1-y)**2+(1-y)**3)
    ids=[sp.simplify(h*sp.diff(phase,y)-k),
         sp.simplify((sp.I*O*h*sp.diff(k,y)+O**2*(1-k**2)-V)/h
                     -(4*sp.I*O*(1-y)+4*O**2*(2-y)-L-1+y))]
    assert ids==[0,0]
    # Verify exact coefficient representation on symbolic basis polynomials.
    x,lam=sp.symbols('x lam')
    A0,A1,A2=coefficient_matrices(2,8)
    for j in range(6):
        v=sp.chebyshevt(j,x)
        R=(1-x-x*x+x**3)/2*sp.diff(v,x,2)+((-1-2*x+3*x*x)/2+2*lam*(1+2*x-x*x))*sp.diff(v,x)+(2*lam*(1-x)-2*lam**2*(3-x)-sp.Rational(13,2)+x/2)*v
        matrix=sum((sp.Rational(float(A0[k,j]))+lam*sp.Rational(float(A1[k,j]))-lam**2*sp.Rational(float(A2[k,j])))*sp.chebyshevt(k,x) for k in range(9))
        assert sp.expand(R-matrix)==0
    return {'operator_identities':[str(x) for x in ids],
            'exact_symbolic_Chebyshev_basis_checks':6,
            'coefficient_assembly':'integer arithmetic, final exact integer/32 entries'}


def algorithm_checks():
    import mpmath as mp
    mp.mp.dps=60
    mats=coefficient_matrices(1,8)
    A=[[[mp.mpf(float(x)) for x in row] for row in aa] for aa in mats]
    lam=mp.mpc('0.27','0.41')
    ld,_=logdet_derivative(A,lam)
    Q=mp.matrix([[A[0][i][j]+lam*A[1][i][j]-lam**2*A[2][i][j] for j in range(9)] for i in range(9)])
    Qp=mp.matrix([[A[1][i][j]-2*lam*A[2][i][j] for j in range(9)] for i in range(9)])
    prod=Q**-1*Qp
    independent=mp.fsum(prod[j,j] for j in range(9))
    error=abs(ld-independent)
    assert error<mp.mpf('1e-50')
    return {'Hessenberg_logdet_derivative_vs_inverse_trace_absolute_error':mp.nstr(error,12),
            'small_matrix_dimension':9,'working_decimal_digits':60}


def metadata():
    return {'field':'massless scalar Schwarzschild',
       'conventions':{'Omega':'2M omega','time':'exp(-i omega t)','z':'r/(2M)',
                      'y':'(z-1)/z','lambda':'i Omega'},
       'method':'Chebyshev-coefficient Lobatto collocation and arbitrary-precision finite-matrix root refinement',
       'limitations':[
         'No interval error enclosures; grid differences are empirical diagnostics.',
         'Smoothness at compactified infinity does not remove exponentially flat incoming contamination in the damped half-plane.',
         'Outgoing scattering continuation must be checked independently.',
         'Labels track initial low-damping positive-frequency candidates; they are provisional, not a completeness theorem.',
         'High working precision controls algebraic root error, not spectral truncation error.',
         'Raw nodal double-precision D^2 formulation loses accuracy much earlier than coefficient form.'
       ],'symbolic_checks':symbolic_checks(),'algorithm_checks':algorithm_checks()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ells',nargs='+',type=int)
    parser.add_argument('--grids',nargs='+',type=int)
    parser.add_argument('--overtones',type=int)
    parser.add_argument('--dps',type=int)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    if all(value is None for value in (args.ells,args.grids,args.overtones,args.dps)):
        from qnm_spectral_release import recompute
        recompute(args.output or here/'qnm_spectral_results.json')
        return
    # Preserve customizable base-grid studies and keep them separate from the
    # released complete benchmark, which is consumed by qnm_compare.py.
    args.ells = args.ells if args.ells is not None else [0,1,2]
    args.grids = args.grids if args.grids is not None else [40,60,80,100,160,240,320]
    args.overtones = args.overtones if args.overtones is not None else 3
    args.dps = args.dps if args.dps is not None else 50
    args.output = args.output or here/'qnm_spectral_custom.json'
    args.output.parent.mkdir(parents=True,exist_ok=True)
    out=metadata()
    out.update(grids=args.grids,working_decimal_digits=args.dps,results=[],initial_double_spectra={})
    for ell in args.ells:
        candidates=solve(ell,args.grids[0])
        out['initial_double_spectra'][str(ell)]=[[float(v.real),float(v.imag)] for v in candidates]
        seeds=candidates[:args.overtones]
        previous=[None]*args.overtones
        for N in args.grids:
            for overtone in range(args.overtones):
                started=time.time()
                row=refine(ell,N,seeds[overtone],args.dps)
                current=complex(*map(float,row['omega_strings']))
                row.update(provisional_overtone=overtone,
                           previous_grid_difference=None if previous[overtone] is None else abs(current-previous[overtone]),
                           elapsed_seconds=time.time()-started)
                previous[overtone]=current
                seeds[overtone]=current
                out['results'].append(row)
                Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
                print(f'ell={ell} n~={overtone} N={N}: {current}',flush=True)
    print(Path(args.output).resolve())

if __name__=='__main__':main()
