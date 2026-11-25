#!/usr/bin/env python3
"""Explicit free-tail defect correction for the incoming Schwarzschild amplitude.

A_corr=A_N+(a I0+b Iexp)/(2iOmega), a=h_R-hx_R/(2iOmega),
b=hx_R/(2iOmega). Iexp is a rapidly decaying contour quadrature.
Remaining analytic tail bound: (|a|+|b|)expm1(kV_R)F_R/(2|Omega|).
Transport, floating-point, and quadrature errors are NOT certified.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import mpmath as mp
import qnm_contour_mp as contour
import qnm_tail_optimization as bounds
import qnm_forward_coefficient as fwd


def residual_data(omega,ell,c,R,N):
    coeff=contour.asymptotic_coefficients(omega,ell,N)
    C=(N*(N+1)-ell*(ell+1))*coeff[N]-N*N*coeff[N-1]
    D=-(N+1)**2*coeff[N]
    z=mp.mpf(c)+mp.j*mp.mpf(R)
    return C,D,z


def defect_integrals(omega,ell,c,R,N,breaks=(0,1,4,12,mp.inf),method='tanh-sinh',check_i0=False):
    """Normalized quadrature retains relative accuracy for tiny defects."""
    omega=mp.mpc(omega);c=mp.mpf(c);R=mp.mpf(R)
    C,D,zR=residual_data(omega,ell,c,R,N)
    I0=C*zR**(-N-1)/(N+1)+D*zR**(-N-2)/(N+2)
    dom=bounds.Domain(omega);_,_,_,mu=dom.phase(c,R)
    scale=1/(mu+(N+2)*R/(c*c+R*R))
    T0=C*zR**(-N-2);T1=D*zR**(-N-3);normalizer=abs(T0)+abs(T1)
    count=0
    def normalized(v,oscillatory=True):
        nonlocal count
        count+=1
        if mp.isinf(v): return mp.mpc(0)
        u=scale*v;zz=zR+mp.j*u;ratio=zR/zz
        defect=(T0/normalizer)*ratio**(N+2)+(T1/normalizer)*ratio**(N+3)
        if oscillatory:
            delta_x=mp.j*u+mp.log1p(mp.j*u/(zR-1))
            defect*=mp.exp(2*mp.j*omega*delta_x)
        return defect
    if not normalizer:
        return I0,mp.mpc(0),{'evaluations':0,'scale':scale,'I0_quadrature_error':0}
    Iexp=mp.j*scale*normalizer*mp.quad(normalized,list(breaks),method=method)
    check=None
    if check_i0:
        I0quad=mp.j*scale*normalizer*mp.quad(lambda v:normalized(v,False),list(breaks),method=method)
        check=abs(I0quad-I0)/max(abs(I0),mp.mpf('1e-200'))
    return I0,Iexp,{'evaluations':count,'scale':scale,'normalizer':normalizer,
                     'I0_quadrature_relative_error':check,'quadrature_method':method,
                     'breakpoints':[str(x) for x in breaks]}


def corrected_coefficient(omega,ell,settings,diagnostics=False,breaks=(0,1,4,12,mp.inf)):
    omega=mp.mpc(omega)
    AN,diag=fwd.forward_coefficient(omega,ell,diagnostics=True,**settings)
    h,hx=diag['h_R'],diag['hx_R']
    a=h-hx/(2*mp.j*omega);b=hx/(2*mp.j*omega)
    I0,Iexp,qd=defect_integrals(omega,ell,settings['anchor'],settings['radius'],settings['asympt_order'],breaks)
    correction=(a*I0+b*Iexp)/(2*mp.j*omega)
    A=AN+correction
    if diagnostics:
        rem=(abs(a)+abs(b))*mp.expm1(diag['tail_kernel_bound']*diag['tail_potential_bound'])*diag['residual_integral_bound']/(2*abs(omega))
        return A,{'uncorrected_A_N':AN,'correction':correction,'a':a,'b':b,'I0':I0,'Iexp':Iexp,
                  'remaining_tail_bound_conditional':rem,'quadrature':qd,**diag}
    return A


def root_search(ell,seed,settings,domain,dps=40):
    calls=0;trace=[]
    with mp.workdps(dps):
        initial=mp.mpc(*seed)
        def fn(w):
            nonlocal calls
            calls+=1
            if abs(w-domain.center)>domain.radius:raise ValueError('Root left chosen frequency disk')
            A=corrected_coefficient(w,ell,settings)
            trace.append({'omega':w,'abs_A_corrected':abs(A)})
            return A
        w=mp.findroot(fn,(initial,initial+mp.mpc('.001','.001')),solver='secant',tol=mp.mpf('1e-29'),maxsteps=20,verify=True)
    return w,calls,trace


def campaign(output='analysis/qnm_free_correction_results.json',dps=40):
    mp.mp.dps=dps
    prior=json.loads(Path(__file__).with_name('qnm_forward_coefficient_results.json').read_text())
    data={'frequency':'Omega=2Momega','working_dps':dps,'method':'explicit free-tail defect correction',
          'qualification':'Analytic remainder conditions on exact forward state and exact integrals. Numerical transport and quadrature errors are separately checked but uncertified.','modes':[]}
    for previous in prior['modes']:
        ell,n,c=previous['ell'],previous['n'],previous['anchor'];seed=previous['seed']
        dom=bounds.Domain(mp.mpc(*seed),'.01');wr=mp.mpc(previous['CF_frequency']['real'],previous['CF_frequency']['imag'])
        rec={'ell':ell,'n':n,'CF_frequency':wr,'radius_diagnostics':[]}
        # Reference is an independent higher-resolution backward evaluation at
        # the same CF frequency; zero is NOT assumed to be exact at finite-depth CF.
        with mp.workdps(55):
            refcfg={'anchor':c,'radius':140,'asympt_order':65,'taylor_order':56,'step_fraction':'.18','horizon_terms':180,'scaled_horizon':True}
            refA,refd=contour.incoming_amplitude(wr,ell,diagnostics=True,**refcfg)
        rec['higher_resolution_reference']={'A':refA,'settings':refcfg,'dps':55,'qualification':'Numerical reference, not certified ground truth.'}
        for priorrun in previous['radius_diagnostics']:
            cfg=priorrun['settings'];R,N=cfg['radius'],cfg['asympt_order']
            st=time.perf_counter();A,diag=corrected_coefficient(wr,ell,cfg,diagnostics=True);elapsed=time.perf_counter()-st
            with mp.workdps(55):
                I0h,Ieh,qh=defect_integrals(wr,ell,c,R,N,breaks=(0,mp.mpf('.5'),2,8,24,mp.inf),check_i0=(R==20))
            qchange=abs(diag['b'])*abs(Ieh-diag['Iexp'])/(2*abs(wr))
            unerr=abs(diag['uncorrected_A_N']-refA);correrr=abs(A-refA)
            rr={'settings':cfg,'A_corrected_at_CF':A,'abs_A_corrected_at_CF':abs(A),
                'uncorrected_discrepancy_from_reference':unerr,'corrected_discrepancy_from_reference':correrr,
                'observed_discrepancy_reduction':unerr/correrr if correrr else None,
                'wall_seconds':elapsed,'quadrature_refinement':{'dps':55,'I0':I0h,'Iexp':Ieh,
                  'Iexp_change':abs(Ieh-diag['Iexp']),'induced_coefficient_change_estimate':qchange,**qh},**diag}
            rec['radius_diagnostics'].append(rr)
            print('free correction',ell,n,'R,N',R,N,'uncorrected discrepancy',mp.nstr(unerr,8),'corrected',mp.nstr(correrr,8),'bound',mp.nstr(diag['remaining_tail_bound_conditional'],8),'seconds',elapsed,flush=True)
        data['modes'].append(rec);bounds.dump(output,data)
        # Two economical root examples, always independently rounded seeds.
        chosenR=40 if (ell,n)==(0,0) else 20 if (ell,n)==(0,2) else None
        if chosenR is not None:
            cfg=next(x['settings'] for x in rec['radius_diagnostics'] if x['settings']['radius']==chosenR)
            st=time.perf_counter();w,calls,trace=root_search(ell,seed,cfg,dom,dps);elapsed=time.perf_counter()-st
            A,diag=corrected_coefficient(w,ell,cfg,diagnostics=True)
            with mp.workdps(55):
                Ah,dh=corrected_coefficient(w,ell,cfg,diagnostics=True,breaks=(0,mp.mpf('.5'),2,8,24,mp.inf))
            rec['corrected_root_search']={'settings':cfg,'seed':seed,'root':w,'calls':calls,'trace':trace,
                'wall_seconds':elapsed,'abs_corrected_coefficient':abs(A),'difference_from_CF':abs(w-wr),
                'root_diagnostics':diag,'55_digit_check':{'A':Ah,'change':abs(Ah-A)},
                'same_endpoint_uncorrected_run':next(x for x in previous['root_searches'] if x['settings']['radius']==chosenR)}
            print('corrected ROOT',ell,n,'R',chosenR,'seconds',elapsed,'CFdiff',mp.nstr(abs(w-wr),9),'bound',mp.nstr(diag['remaining_tail_bound_conditional'],9),flush=True)
        bounds.dump(output,data)
    return data


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',default='analysis/qnm_free_correction_results.json')
    ap.add_argument('--dps',default=40,type=int)
    args=ap.parse_args();campaign(args.output,args.dps)
