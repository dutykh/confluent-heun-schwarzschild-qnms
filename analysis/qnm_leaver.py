#!/usr/bin/env python3
"""Scalar Schwarzschild QNMs from independently derived Leaver recurrence.

Conventions: e^{-i omega t}, z=r/(2M), Omega=2M*omega, x=z+log(z-1).
The dimensionless radial equation is u_xx+[Omega**2-U_l]u=0,
U_l=(1-1/z)[l(l+1)/z**2+1/z**3].  Substitute
  t=(z-1)/z,
  u=exp(i*Omega*z)*t**(-i*Omega)*(1-t)**(-i*Omega)*sum(a_n*t**n).
The resulting series equation is
  t(1-t)^2 S'' + [(3-4iO)t^2+(8iO-4)t+1-2iO]S'
  +[(1-4iO-4O^2)t+8O^2+4iO-l(l+1)-1]S=0.
Thus alpha_n a_{n+1}+beta_n a_n+gamma_n a_{n-1}=0, a_{-1}=0.
The minimal-solution continued fraction implements the outgoing condition.
We use overtone-index inversions and a_{N+1}=0 at finite depth N.
Depth convergence is empirical, not an interval-certified truncation bound.
No contour-determinant outputs are read or used as reference values.

Dependency: mpmath (symbolic self-check additionally uses sympy).
Run: python3 analysis/qnm_leaver.py --output analysis/qnm_leaver_results.json
"""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
import mpmath as mp


def coefficients(n, omega, ell):
    n=mp.mpf(n)
    alpha=(n+1)*(n+1-2j*omega)
    beta=-2*n*n+(8j*omega-2)*n+8*omega*omega+4j*omega-ell*(ell+1)-1
    gamma=(n-2j*omega)**2
    return alpha,beta,gamma


def cf_with_derivative(omega, ell, inversion, depth):
    """nth-inverted CF and exact derivative, zero ratio beyond index depth."""
    if depth <= inversion:
        raise ValueError('depth must exceed inversion')
    I=mp.j
    # Coefficients and derivatives evaluated through low-degree polynomials.
    bc=8*omega*omega+4*I*omega-ell*(ell+1)-1
    bl=8*I*omega-2
    bdc=16*omega+4*I
    def coeff(n):
        n=mp.mpf(n)
        a=(n+1)*(n+1-2*I*omega)
        b=-2*n*n+bl*n+bc
        g=(n-2*I*omega)**2
        da=-2*I*(n+1)
        db=8*I*n+bdc
        dg=-4*I*(n-2*I*omega)
        return a,b,g,da,db,dg
    # Left finite fraction, eliminating a_0,...,a_{inversion-1}.
    ap,left,g,apd,leftd,gd=coeff(0)
    for n in range(1,inversion+1):
        a,b,g,ad,bd,gd=coeff(n)
        numer=ap*g; numerd=apd*g+ap*gd
        leftd=bd-numerd/left+numer*leftd/(left*left)
        left=b-numer/left
        ap,apd=a,ad
    # Right fraction T_k=alpha_k gamma_{k+1}/(beta_{k+1}-T_{k+1}).
    tail=mp.mpc(0); taild=mp.mpc(0)
    a,bn,gn,ad,bdn,gdn=coeff(depth)
    for n in range(depth-1,inversion-1,-1):
        a,b,g,ad,bd,gd=coeff(n)
        den=bn-tail; dend=bdn-taild
        numer=a*gn; numerd=ad*gn+a*gdn
        taild=numerd/den-numer*dend/(den*den)
        tail=numer/den
        bn,gn,bdn,gdn=b,g,bd,gd
    return left-tail,leftd-taild


def newton_root(seed,ell,inversion,depth,dps=60):
    with mp.workdps(dps):
        w=mp.mpc(seed)
        tolerance=mp.mpf(10)**(-(dps-12))
        for it in range(30):
            F,dF=cf_with_derivative(w,ell,inversion,depth)
            step=F/dF
            w-=step
            if abs(step)<tolerance:
                F,dF=cf_with_derivative(w,ell,inversion,depth)
                return +w,+F,+dF,it+1
        raise ArithmeticError(f'Newton failed ell={ell} n={inversion} N={depth}')


def scomplex(z,digits=48):
    return {'real':mp.nstr(z.real,digits),'imag':mp.nstr(z.imag,digits)}


def symbolic_self_check():
    import sympy as s
    t,w,L=s.symbols('t Omega L')
    bx=t*(1-t)**2
    r=s.I*w/(1-t)**2-s.I*w/t+s.I*w/(1-t)
    U=t*(L*(1-t)**2+(1-t)**3)
    B=s.cancel((2*r*bx*bx+bx*s.diff(bx,t))/bx)
    C=s.cancel((bx*bx*(s.diff(r,t)+r*r)+bx*s.diff(bx,t)*r+w*w-U)/bx)
    expected_B=(3-4*s.I*w)*t*t+(8*s.I*w-4)*t+1-2*s.I*w
    expected_C=(1-4*s.I*w-4*w*w)*t+8*w*w+4*s.I*w-L-1
    assert s.simplify(B-expected_B)==0
    assert s.simplify(C-expected_C)==0
    n=s.symbols('n',integer=True,positive=True)
    ar=(n+1)*n+(1-2*s.I*w)*(n+1)
    br=-2*n*(n-1)+(8*s.I*w-4)*n+8*w*w+4*s.I*w-L-1
    gr=(n-1)*(n-2)+(3-4*s.I*w)*(n-1)+1-4*s.I*w-4*w*w
    assert s.simplify(ar-(n+1)*(n+1-2*s.I*w))==0
    assert s.simplify(br-(-2*n*n+(8*s.I*w-2)*n+8*w*w+4*s.I*w-L-1))==0
    assert s.simplify(gr-(n-2*s.I*w)**2)==0
    # Validate analytically differentiated CF against numerical complex derivative.
    with mp.workdps(40):
        w=mp.mpc('.53','-.61')
        F,dF=cf_with_derivative(w,1,1,80)
        dFn=mp.diff(lambda z:cf_with_derivative(z,1,1,80)[0],w)
        assert abs(dF-dFn)<mp.mpf('1e-35')
    return 'Exact gauge equation and all recurrence coefficients verified symbolically; CF derivative independently checked.'


# Coarse search seeds only; values are not used as target frequencies or tests.
SEEDS={
  (0,0):('.22','-.21'),(0,1):('.17','-.70'),(0,2):('.15','-1.20'),
  (1,0):('.59','-.20'),(1,1):('.53','-.61'),(1,2):('.46','-1.08'),
  (2,0):('.97','-.19'),(2,1):('.93','-.59'),(2,2):('.86','-1.02')}


def campaign(output, max_depth=25600, dps=60, target='1e-19'):
    mp.mp.dps=dps
    results={
        'method':'scalar Schwarzschild Leaver three-term recurrence; overtone-inverted continued fraction',
        'conventions':{'time':'exp(-i omega t)','frequency':'Omega=2M omega','radial':'z=r/(2M)>1','tortoise':'x=r_*/(2M)=z+log(z-1)',
                       'finite_tail':'a_(N+1)=0; depth N includes beta_N',
                       'inversion':'inversion index equals requested overtone n'},
        'recurrence':{'alpha_n':'(n+1)(n+1-2i Omega)',
                      'beta_n':'-2n^2+(8i Omega-2)n+8Omega^2+4iOmega-l(l+1)-1',
                      'gamma_n':'(n-2i Omega)^2'},
        'derivation':'Direct exact substitution of stated gauge into radial ODE; symbolic self-check included.',
        'self_check':symbolic_self_check(),
        'working_dps':dps,'maximum_depth':max_depth,'target_successive_depth_difference':target,
        'caveat':'Depth convergence and residuals are empirical; no certified bound on the infinite-CF truncation error. Coarse seeds are not ground truth.',
        'modes':[]}
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    for (ell,n),(sr,si) in SEEDS.items():
        start=time.time(); prev=None; w=mp.mpc(sr,si); records=[]; depth=50
        while depth<=max_depth:
            w,F,dF,iterations=newton_root(w,ell,n,depth,dps)
            diff=abs(w-prev) if prev is not None else None
            record={'depth':depth,'omega':scomplex(w), 'absolute_cf_residual':mp.nstr(abs(F),12),
                    'absolute_cf_derivative':mp.nstr(abs(dF),12),'newton_iterations':iterations,
                    'successive_root_difference':None if diff is None else mp.nstr(diff,12)}
            records.append(record)
            print(f'l={ell} n={n} N={depth} Omega={mp.nstr(w,28)} delta={None if diff is None else mp.nstr(diff,5)}',flush=True)
            # Require a sizeable depth even where convergence is exceptionally quick.
            if depth>=800 and diff is not None and diff<mp.mpf(target): break
            prev=+w; depth*=2
        final_depth=records[-1]['depth']
        # Precision check at the SAME tail: separates arithmetic from truncation.
        wh,Fh,dFh,ith=newton_root(w,ell,n,final_depth,dps+25)
        precision_difference=abs(wh-w)
        # Unrooted validation at twice the final depth quantifies residual drift.
        F2,dF2=cf_with_derivative(wh,ell,n,2*final_depth)
        mode={'ell':ell,'overtone':n,'coarse_seed':[sr,si], 'depth_study':records,
              'final_depth':final_depth,'omega':scomplex(wh,65),'Momega':scomplex(wh/2,min(55,dps-5)),
              'precision_check':{'dps_low':dps,'dps_high':dps+25,'absolute_root_difference':mp.nstr(precision_difference,12),'cf_residual_high_precision':mp.nstr(abs(Fh),12)},
              'double_depth_evaluation':{'depth':2*final_depth,'absolute_cf_residual':mp.nstr(abs(F2),12),'newton_frequency_correction_estimate':mp.nstr(abs(F2/dF2),12)},
              'reached_successive_depth_target':records[-1]['successive_root_difference'] is not None and mp.mpf(records[-1]['successive_root_difference'])<mp.mpf(target),
              'elapsed_seconds':round(time.time()-start,3)}
        results['modes'].append(mode)
        output.write_text(json.dumps(results,indent=2)+'\n')
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',default='analysis/qnm_leaver_results.json')
    p.add_argument('--max-depth',type=int,default=25600)
    p.add_argument('--dps',type=int,default=60)
    p.add_argument('--target',default='1e-19')
    args=p.parse_args()
    campaign(args.output,args.max_depth,args.dps,args.target)
