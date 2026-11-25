#!/usr/bin/env python3
"""Independent Chebyshev-collocation experiment for scalar Schwarzschild QNMs.

Conventions: z=r/(2M), Omega=2M*omega, exp(-i*omega*t),
y=(z-1)/z in [0,1]. Factor u=exp(i*Omega*z)*(z-1)^(-i*Omega)
*z^(2*i*Omega)*v(y). Polynomial eigenvalue in lambda=i*Omega:
  h v''+(h'+2 lambda k)v'+[lambda k'-4 lambda^2(2-y)-L-1+y]v=0,
h=y(1-y)^2, k=1-2(1-y)^2, L=ell(ell+1).

CAVEAT: finite/smooth endpoint residuals alone do not uniquely impose outgoing
analytic continuation at infinity in Im(Omega)<0, because the incoming solution
relative to the outgoing factor is exponentially flat on this compactified real
interval. Thus stable collocation roots are independent numerical consistency
checks, not a proof of the intended continued scattering boundary condition.
No Leaver data or externally tabulated frequencies enter this script.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eig
from scipy.fft import dct
from numpy.polynomial.chebyshev import chebval, chebder


def cheb_grid(n):
    j=np.arange(n+1)
    y=(1-np.cos(np.pi*j/n))/2
    weights=(-1.)**j
    weights[[0,-1]]*=0.5
    delta=y[:,None]-y[None,:]
    np.fill_diagonal(delta,1.)
    D=weights[None,:]/weights[:,None]/delta
    np.fill_diagonal(D,0.)
    np.fill_diagonal(D,-np.sum(D,axis=1))
    return y,D


def matrices(ell,n):
    y,D=cheb_grid(n)
    h=y*(1-y)**2
    hp=(1-y)*(1-3*y)
    k=1-2*(1-y)**2
    L=ell*(ell+1)
    A0=h[:,None]*(D@D)+hp[:,None]*D-np.diag(L+1-y)
    A1=2*k[:,None]*D+np.diag(4*(1-y))
    A2=np.diag(4*(2-y))
    return y,D,A0,A1,A2


def diagnostics(lam,vec,ell,n,A0,A1,A2):
    v=vec[:n+1].astype(complex)
    v/=np.max(np.abs(v))
    p0=A0@v
    p1=lam*(A1@v)
    p2=lam**2*(A2@v)
    residual=np.linalg.norm(p0+p1-p2,np.inf)/(
        np.linalg.norm(p0,np.inf)+np.linalg.norm(p1,np.inf)+np.linalg.norm(p2,np.inf))
    # DCT sees descending Chebyshev x; convert to ascending x=2*y-1.
    coeff=dct(v,type=1)/n*((-1.)**np.arange(n+1))
    coeff[[0,-1]]*=.5
    tail=np.max(np.abs(coeff[-min(10,n//4):]))/np.max(np.abs(coeff))
    xx=np.cos(np.pi*(np.arange(n+13)+.5)/(n+13))
    yy=(xx+1)/2
    vv=chebval(xx,coeff)
    vp=2*chebval(xx,chebder(coeff))
    vpp=4*chebval(xx,chebder(coeff,2))
    t0=yy*(1-yy)**2*vpp
    t1=((1-yy)*(1-3*yy)+2*lam*(1-2*(1-yy)**2))*vp
    t2=(4*lam*(1-yy)-4*lam**2*(2-yy)-ell*(ell+1)-1+yy)*vv
    off=np.max(np.abs(t0+t1+t2))/(np.max(np.abs(t0))+np.max(np.abs(t1))+np.max(np.abs(t2)))
    # Endpoints are the differential equation limits, not independently imposed BCs.
    v0=chebval(-1,coeff); vp0=2*chebval(-1,chebder(coeff))
    v1=chebval(1,coeff); vp1=2*chebval(1,chebder(coeff))
    bc0=(1-2*lam)*vp0+(4*lam-8*lam**2-ell*(ell+1)-1)*v0
    bc1=2*lam*vp1+(-4*lam**2-ell*(ell+1))*v1
    return {'quadratic_relative_residual':float(residual),
            'off_grid_relative_residual':float(off),
            'last_ten_chebyshev_tail_ratio':float(tail),
            'endpoint_absolute_residuals':[float(abs(bc0)),float(abs(bc1))]}


def solve(ell,n,omega_real_max=2.,omega_damping_max=2.):
    y,D,A0,A1,A2=matrices(ell,n)
    size=n+1
    Z=np.zeros((size,size)); I=np.eye(size)
    M=np.block([[Z,I],[A0,A1]])
    B=np.block([[I,Z],[Z,A2]])
    vals,vecs=eig(M,B,check_finite=False)
    good=[]
    for j,lam in enumerate(vals):
        om=-1j*lam
        if np.isfinite(om) and 1e-5<om.real<omega_real_max and -omega_damping_max<om.imag<-1e-5:
            row={'omega':[float(om.real),float(om.imag)],
                 'lambda':[float(lam.real),float(lam.imag)]}
            row.update(diagnostics(lam,vecs[:,j],ell,n,A0,A1,A2))
            good.append(row)
    good.sort(key=lambda row:-row['omega'][1])
    return good


def symbolic_check():
    import sympy as s
    y,O,L=s.symbols('y O L')
    z=1/(1-y)
    h=y*(1-y)**2
    q=1/(1-y)-s.log(y)-s.log(1-y)
    k=1-2*(1-y)**2
    check1=s.simplify(h*s.diff(q,y)-k)
    potential=y*(L*(1-y)**2+(1-y)**3)
    c0=s.simplify((s.I*O*h*s.diff(k,y)+O**2*(1-k**2)-potential)/h)
    expected=4*s.I*O*(1-y)+4*O**2*(2-y)-L-1+y
    check2=s.simplify(c0-expected)
    check3=s.simplify(s.diff(h,y)-(1-y)*(1-3*y))
    assert all(x==0 for x in (check1,check2,check3))
    yy,D=cheb_grid(24)
    differentiation_error=float(np.max(np.abs(D@(yy**3)-3*yy**2)))
    assert differentiation_error<1e-11
    return {'symbolic_operator_identities':[str(x) for x in (check1,check2,check3)],
            'grid_cubic_derivative_max_error':differentiation_error}


def match_tracks(runs,grids,threshold):
    # Identify branches converged on the first three grids, independently of any
    # tabulated frequencies. Do not hide subsequent loss of accuracy in double precision.
    base=grids[0]
    tracks=[]
    for candidate in runs[str(base)]:
        om0=complex(*candidate['omega'])
        seq=[]
        for n in grids:
            choices=runs[str(n)]
            if not choices: break
            nearest=min(choices,key=lambda row:abs(complex(*row['omega'])-om0))
            seq.append({'N':n,**nearest,'distance_from_first':float(abs(complex(*nearest['omega'])-om0))})
        if len(seq)!=len(grids): continue
        first_three=[complex(*row['omega']) for row in seq[:3]]
        if max(abs(om-first_three[0]) for om in first_three)<threshold:
            pair_deltas=[abs(complex(*seq[j]['omega'])-complex(*seq[j-1]['omega'])) for j in range(1,len(seq))]
            best_index=int(np.argmin(pair_deltas))+1
            tracks.append({'branch_sorted_by_damping':None,
                           'provisional_overtone_label':None,
                           'grid_sequence':seq,
                           'successive_grid_differences':[float(x) for x in pair_deltas],
                           'closest_pair_grids':[grids[best_index-1],grids[best_index]],
                           'closest_pair_difference':float(pair_deltas[best_index-1]),
                           'representative_omega':seq[best_index]['omega']})
    tracks.sort(key=lambda tr:-tr['representative_omega'][1])
    for j,tr in enumerate(tracks):
        tr['branch_sorted_by_damping']=j
        tr['provisional_overtone_label']=j
    return tracks


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--grids',nargs='+',type=int,default=[40,60,80,100,120])
    p.add_argument('--ells',nargs='+',type=int,default=[0,1,2])
    p.add_argument('--stability-tolerance',type=float,default=1e-5)
    p.add_argument('--output',default='analysis/qnm_spectral_nodal_results.json')
    args=p.parse_args()
    checks=symbolic_check()
    results={'conventions':{'field':'massless scalar','z':'r/(2M)','Omega':'2M omega',
        'time':'exp(-i omega t)','compact_coordinate':'y=(z-1)/z',
        'factor':'u=exp(i Omega z)*(z-1)^(-i Omega)*z^(2 i Omega)*v(y)',
        'eigenparameter':'lambda=i Omega'},
        'limitations':[
          'Double-precision generalized quadratic Chebyshev collocation; no rigorous error enclosures.',
          'Endpoint smoothness alone does not eliminate exponentially flat incoming contamination at infinity for Im(Omega)<0.',
          'Tracks selected by convergence on first three grids, without using Leaver or tabulated frequencies.',
          'Overtone labels are provisional damping order among retained positive-frequency branches, not a completeness theorem.',
          'Closest grid-pair discrepancy is an empirical diagnostic, not an error estimate; higher N may lose accuracy.'
        ],'symbolic_and_grid_checks':checks,'grids':args.grids,'ell':{}}
    for ell in args.ells:
        runs={}
        for n in args.grids:
            runs[str(n)]=solve(ell,n)
            print(f'ell={ell} N={n}: {len(runs[str(n)])} filtered candidates',flush=True)
        tracks=match_tracks(runs,args.grids,args.stability_tolerance)
        results['ell'][str(ell)]={'all_filtered_candidates':runs,'stable_tracks':tracks}
        for tr in tracks[:5]:
            print(' track',ell,tr['provisional_overtone_label'],tr['representative_omega'],
                  'pair difference',tr['closest_pair_difference'],flush=True)
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(results,indent=2)+'\n')
    print(dest.resolve(),flush=True)

if __name__=='__main__':main()
