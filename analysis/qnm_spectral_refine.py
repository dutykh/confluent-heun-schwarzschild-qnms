"""High-precision root refinement of exact Chebyshev collocation matrices.

The matrix is upper Hessenberg. No-pivot elimination and differentiated pivots
compute det'(Q)/det(Q), avoiding costly full eigensolves. Starting guesses are
from the independently computed collocation spectrum, never Leaver data.
All coefficient matrices are dyadic rationals exactly represented as doubles
at the tested resolutions, checked before conversion to arbitrary precision.
"""
import time,json,sys
from pathlib import Path
import mpmath as mp
import numpy as np
from qnm_spectral_chebcoef import coefficient_matrices,solve

def logdet_derivative(A,lam):
    A0,A1,A2=A
    n=len(A0)
    H=[[A0[i][j]+lam*A1[i][j]-lam**2*A2[i][j] for j in range(n)] for i in range(n)]
    Hp=[[A1[i][j]-2*lam*A2[i][j] for j in range(n)] for i in range(n)]
    out=mp.mpc(0)
    minpivot=mp.inf
    for k in range(n):
        p=H[k][k]; dp=Hp[k][k]
        minpivot=min(minpivot,abs(p))
        out+=dp/p
        if k==n-1: break
        fac=H[k+1][k]/p
        facp=(Hp[k+1][k]-fac*dp)/p
        for j in range(k+1,n):
            Hp[k+1][j]-=facp*H[k][j]+fac*Hp[k][j]
            H[k+1][j]-=fac*H[k][j]
    return out,minpivot

def refine(ell,n,omega,dps=50,tol_digits=30):
    mp.mp.dps=dps
    Af=coefficient_matrices(ell,n)
    for a in Af:
        assert np.all(a*32==np.round(a*32)), 'non-dyadic data: exact assembly must be reviewed'
        assert np.all(np.tril(a,-2)==0), 'matrix is not Hessenberg'
    A=[[[mp.mpf(float(x)) for x in row] for row in aa] for aa in Af]
    lam=mp.j*mp.mpc(float(omega.real),float(omega.imag))
    history=[]
    for it in range(12):
        ld,minp=logdet_derivative(A,lam)
        step=1/ld
        lam-=step
        history.append(mp.nstr(abs(step),12))
        if abs(step)<mp.mpf(10)**(-tol_digits): break
    else: raise ArithmeticError('Newton did not converge')
    om=-mp.j*lam
    return {'ell':ell,'N':n,'working_decimal_digits':dps,
            'omega_strings':[mp.nstr(om.real,dps-5),mp.nstr(om.imag,dps-5)],
            'newton_absolute_steps':history,
            'final_step':float(abs(step))}

if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser()
    pa.add_argument('--ells',nargs='+',type=int,default=[0,1,2])
    pa.add_argument('--grids',nargs='+',type=int,default=[40,80,160])
    pa.add_argument('--overtones',nargs='+',type=int,default=[0,1,2])
    pa.add_argument('--dps',type=int,default=50)
    pa.add_argument('--output',default='analysis/qnm_spectral_highprecision.json')
    ar=pa.parse_args()
    result={'method':'exact dyadic Chebyshev-coefficient Lobatto collocation, arbitrary-precision determinant Newton','results':[]}
    for ell in ar.ells:
        seeds=solve(ell,ar.grids[0])
        for n in ar.grids:
            for overtone in ar.overtones:
                start=time.time()
                rr=refine(ell,n,seeds[overtone],ar.dps)
                seeds[overtone]=complex(float(rr['omega_strings'][0]),float(rr['omega_strings'][1]))
                rr['provisional_overtone']=overtone
                rr['elapsed_seconds']=time.time()-start
                result['results'].append(rr)
                Path(ar.output).write_text(json.dumps(result,indent=2)+'\n')
                print(ell,n,overtone,rr['omega_strings'],'seconds',rr['elapsed_seconds'],flush=True)


def collocation_diagnostics(ell,n,omega_strings,dps=60):
    """Residual/tail checks for finite collocation, not scattering-error estimates."""
    mp.mp.dps=dps
    lam=mp.j*mp.mpc(*omega_strings)
    Af=coefficient_matrices(ell,n)
    A=[[[mp.mpf(float(x)) for x in row] for row in aa] for aa in Af]
    Q=[[A[0][i][j]+lam*A[1][i][j]-lam**2*A[2][i][j] for j in range(n+1)] for i in range(n+1)]
    H=[row[:] for row in Q]
    for k in range(n):
        fac=H[k+1][k]/H[k][k]
        for j in range(k+1,n+1):H[k+1][j]-=fac*H[k][j]
    v=[mp.mpc(0)]*(n+1);v[n]=1
    for k in range(n-1,-1,-1):
        v[k]=-mp.fsum(H[k][j]*v[j] for j in range(k+1,n+1))/H[k][k]
    scale=max(map(abs,v));v=[vv/scale for vv in v]
    res=[mp.fsum(Q[i][j]*v[j] for j in range(n+1)) for i in range(n+1)]
    prod=[[mp.fsum(A[k][i][j]*v[j] for j in range(n+1)) for i in range(n+1)] for k in range(3)]
    denom=max(map(abs,prod[0]))+abs(lam)*max(map(abs,prod[1]))+abs(lam)**2*max(map(abs,prod[2]))
    last=((n+1)/mp.mpf(2)-lam)**2*v[n]
    vp0=2*mp.fsum(((-1)**(k-1))*k*k*v[k] for k in range(1,n+1))
    vp1=2*mp.fsum(k*k*v[k] for k in range(1,n+1))
    v0=mp.fsum((-1)**k*v[k] for k in range(n+1));v1=mp.fsum(v)
    b0=(1-2*lam)*vp0+(4*lam-8*lam**2-ell*(ell+1)-1)*v0
    b1=2*lam*vp1+(-4*lam**2-ell*(ell+1))*v1
    return {'ell':ell,'N':n,'omega_strings':omega_strings,
            'normalization':'max absolute Chebyshev coefficient = 1',
            'finite_matrix_relative_residual':mp.nstr(max(map(abs,res))/denom,12),
            'last_ten_coefficient_tail_ratio':mp.nstr(max(map(abs,v[-10:])),12),
            'continuous_ODE_residual_sup_bound':mp.nstr(mp.fsum(map(abs,res))+2*abs(last),12),
            'endpoint_absolute_residuals':[mp.nstr(abs(b0),12),mp.nstr(abs(b1),12)],
            'meaning':'ODE bound follows from exact residual = r_(N+1)*(T_(N+1)-T_(N-1)) plus finite-matrix residual; it is not a QNM error bound.'}
