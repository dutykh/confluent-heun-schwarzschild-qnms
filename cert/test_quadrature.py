"""Numerical regression checks for the analytic interval quadrature enclosure.

The enclosure proof is in quadrature.py. High-precision mpmath quadrature here
is an independent regression oracle, not part of the certificate or its proof.
Run: python -m cert.test_quadrature
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp

from .amplitude import polynomial
from .elementary import backend_info, endpoints
from .quadrature import defect_integrals


def point(raw):
    return mp.make_mpf(raw)


def contains(box, value):
    return (point(box.real._mpi_[0]) <= value.real <= point(box.real._mpi_[1])
            and point(box.imag._mpi_[0]) <= value.imag <= point(box.imag._mpi_[1]))


def real_upper(x):
    return point(x._mpi_[1])


def serialize(value):
    if hasattr(value,'_mpi_'):
        return {'lower_exact':str(endpoints(value)[0]),'upper_exact':str(endpoints(value)[1])}
    if hasattr(value,'_mpci_'):
        return {'real':serialize(value.real),'imag':serialize(value.imag)}
    if isinstance(value,dict):
        return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [serialize(v) for v in value]
    if isinstance(value,Fraction):
        return str(value)
    if isinstance(value,(mp.mpf,mp.mpc)):
        return str(value)
    return value


def residual(omega,ell,N):
    """Point-arithmetic reference recurrence, independent from interval code."""
    aa=[mp.mpc(1)]
    for k in range(N):
        previous=aa[k-1] if k else mp.mpc(0)
        aa.append(((k*(k+1)-ell*(ell+1))*aa[k]-k*k*previous)/(2j*omega*(k+1)))
    c=(N*(N+1)-ell*(ell+1))*aa[N]-N*N*(aa[N-1] if N else 0)
    d=-(N+1)**2*aa[N]
    return c,d


def reference(omega,c,d,N,a,R,T):
    zR=a+1j*R
    scale=abs(c*zR**(-N-2))+abs(d*zR**(-N-3))
    def f(t,exponential=True):
        zz=zR+1j*t
        out=(c*zz**(-N-2)+d*zz**(-N-3))/scale
        if exponential:
            out*=mp.exp(2j*omega*(1j*t+mp.log((zz-1)/(zR-1))))
        return out
    I0=1j*scale*mp.quad(lambda t:f(t,False),[0,1,4,16,64,mp.inf])
    Ie=1j*scale*mp.quad(f,[0,1,4,16,64,mp.inf])
    local=1j*scale*mp.quad(f,[0,T/4,T/2,T])
    tail=1j*scale*mp.quad(f,[T,T+1,T+4,T+16,T+64,mp.inf])
    return I0,Ie,local,tail


def campaign():
    mp.mp.dps=100
    mp.iv.dps=65
    results={'backend':backend_info(),'interval_dps':mp.iv.dps,
             'reference_dps':mp.mp.dps,
             'qualification':'References are numerical regression checks. The rigorous enclosure uses analytic Taylor-Cauchy and infinite-tail bounds with outward arithmetic.',
             'point_cases':[]}
    cases=[('generic',None,2,('.22','-.21'),'1.5'),
           ('scalar_l0',0,11,('.22','-.21'),'1.5'),
           ('scalar_l2',2,37,('.97','-.19'),'1.5'),
           ('scalar_overtone',0,47,('.15','-1.20'),'5')]
    for label,ell,N,ww,anchor in cases:
        w=mp.iv.mpc(*ww)
        if ell is None:
            c,d=mp.iv.mpc('.3','.5'),mp.iv.mpc('-.7','.2')
            cp,dp=mp.mpc('.3','.5'),mp.mpc('-.7','.2')
        else:
            _,_,c,d=polynomial(w,ell,N,mp.iv.mpc(anchor,20))
            cp,dp=residual(mp.mpc(*ww),ell,N)
        st=time.perf_counter()
        i0,ie,diag=defect_integrals(w,c,d,N,anchor,20,rtol='1e-30',return_diagnostics=True)
        elapsed=time.perf_counter()-st
        T=mp.mpf(diag['truncation_endpoint'].numerator)/diag['truncation_endpoint'].denominator
        r0,re,rl,rt=reference(mp.mpc(*ww),cp,dp,N,mp.mpf(anchor),mp.mpf(20),T)
        side=mp.iv.mpf([-1,1])*diag['local_remainder']
        local_box=diag['finite_taylor_sum']+mp.iv.mpc(side,side)
        checks={'I0_contains_reference':contains(i0,r0),
                'Iexp_contains_reference':contains(ie,re),
                'finite_integral_contains_reference':contains(local_box,rl),
                'reference_tail_below_bound':abs(rt)<=real_upper(diag['tail_remainder']),
                'analytic_error_below_target':real_upper(diag['analytic_error'])<=point(diag['target']._mpi_[0])}
        assert all(checks.values()),(label,checks)
        row={'label':label,'ell':ell,'N':N,'omega':ww,'anchor':anchor,'radius':20,
             'seconds':elapsed,'I0':i0,'Iexp':ie,'diagnostics':diag,'checks':checks,
             'reference_Iexp':str(re),'reference_tail_over_bound':str(abs(rt)/real_upper(diag['tail_remainder']))}
        results['point_cases'].append(row)
        print(label,'PASS',elapsed,'seconds',diag['panels'],'panels',flush=True)
    # A nontrivial frequency rectangle, with c,d also interval-dependent.
    N,ell=11,0
    wr=('0.2199999999','0.2200000001');wi=('-0.2100000001','-0.2099999999')
    w=mp.iv.mpc(mp.iv.mpf(wr),mp.iv.mpf(wi))
    _,_,c,d=polynomial(w,ell,N,mp.iv.mpc('1.5',20))
    i0,ie,diag=defect_integrals(w,c,d,N,'1.5',20,rtol='1e-30',return_diagnostics=True)
    checks=[]
    for x,y in [(wr[0],wi[0]),(wr[0],wi[1]),(wr[1],wi[0]),(wr[1],wi[1]),('.22','-.21')]:
        wp=mp.mpc(x,y);cp,dp=residual(wp,ell,N)
        r0,re,_,_=reference(wp,cp,dp,N,mp.mpf('1.5'),mp.mpf(20),mp.mpf(128))
        check={'omega':(x,y),'I0_contains_reference':contains(i0,r0),'Iexp_contains_reference':contains(ie,re)}
        assert check['I0_contains_reference'] and check['Iexp_contains_reference'],check
        checks.append(check)
    results['frequency_box']={'real':wr,'imag':wi,'I0':i0,'Iexp':ie,'checks':checks,'diagnostics':diag}
    zero=defect_integrals(mp.iv.mpc('.22','-.21'),mp.iv.mpc(0),mp.iv.mpc(0),0,'1.5',20)
    assert zero[0]==0 and zero[1]==0
    try:
        defect_integrals(mp.iv.mpc('.01','-1'),mp.iv.mpc(1),mp.iv.mpc(1),1,'1.5',20)
    except ValueError:
        results['invalid_damping_rejected']=True
    else:
        raise AssertionError('Damping condition was not checked')
    results['zero_case_exact']=True
    out=Path(__file__).with_name('quadrature_results.json')
    out.write_text(json.dumps(serialize(results),indent=2)+'\n')
    print('All quadrature regression checks passed:',out,flush=True)
    return results


if __name__=='__main__':
    campaign()
