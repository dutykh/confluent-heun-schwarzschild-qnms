#!/usr/bin/env python3
"""Sharper analytic outgoing-tail bounds and fixed-domain cost experiments.

No earlier deliverable is modified. Uses qnm_contour_mp solely as an unchanged
contour evaluator. Bounds are mathematical exact-arithmetic expressions evaluated
with mpmath, not outward-rounded certified numerics. They exclude Taylor and
floating-point propagation errors. Domain settings are selected once on a fixed
frequency disk and kept fixed throughout each independent secant root search.
"""
from __future__ import annotations
import argparse
import functools
import json
import time
from pathlib import Path
import mpmath as mp
import qnm_contour_mp as contour


def power_tail(c,R,p):
    """Integral_R^inf (c*c+s*s)^(-p/2) ds, p>1."""
    c,R=mp.mpf(c),mp.mpf(R)
    x=c*c/(c*c+R*R)
    return c**(1-p)*mp.betainc(mp.mpf(p-1)/2,mp.mpf('.5'),0,x)/2


def weighted_power_upper(c,R,p,mu):
    """Rigorous convex-tangent upper bound on exponential-weighted power tail."""
    c,R,mu=map(mp.mpf,(c,R,mu))
    r=mp.sqrt(c*c+R*R); k=R/(r*r); a=mu/k
    # expint at large integer order/argument can use cancellation-prone
    # recurrences. Compute the positive scaled quantity with guarded precision.
    extra=max(30,int(mp.ceil(a/mp.log(10)))+20)
    with mp.workdps(mp.mp.dps+extra):
        guarded=r**(-p)*mp.exp(a)*mp.expint(p,a)/k
        if guarded<=0:
            raise ArithmeticError('Lost positivity in exponential integral')
    return min(power_tail(c,R,p),+guarded)


class Domain:
    def __init__(self,center,radius=0):
        self.center=mp.mpc(center); self.radius=mp.mpf(radius)
        self.p=self.center.real-self.radius
        self.q=-self.center.imag+self.radius
        self.omin=abs(self.center)-self.radius
        if self.p<=0 or self.center.imag+self.radius>=0 or self.omin<=0:
            raise ValueError('This study uses disks strictly in Re>0, Im<0')
    def phase(self,c,R):
        c,R=mp.mpf(c),mp.mpf(R); a=c-1
        rho0=1/(2*mp.sqrt(c*a))
        rhot=rho0 if R<=mp.sqrt(c*a) else R/(R*R+c*a)
        ratio=self.q/self.p
        def K(rho):
            if ratio*rho>=1: raise ValueError('Frequency disk fails phase cone condition')
            kap=(ratio+rho)/(1-ratio*rho)
            S=mp.sqrt(1+kap*kap)
            return S/(1+S)/self.omin
        O=mp.sqrt(self.p*self.p+self.q*self.q)
        mustar=2*(self.p-(O-self.p)/(2*a))
        sstar=a*(self.p+O)/self.q
        ss=max(R,sstar)
        mut=2*(self.p+(self.p*a-self.q*ss)/(a*a+ss*ss))
        return K(rho0),K(rhot),mustar,mut


def coefficient_data(domain,ell,max_order):
    if not domain.radius:
        return contour.asymptotic_coefficients(domain.center,ell,max_order),False
    L=ell*(ell+1); a=[mp.mpf(1)]
    for k in range(max_order):
        prev=a[k-1] if k else 0
        a.append((abs(k*(k+1)-L)*a[k]+k*k*prev)/(2*domain.omin*(k+1)))
    return a,True


def tail_bound(domain,ell,c,R,N,coeffs=None,level='improved',refine=True):
    """Uniform outgoing state norm bound Em+Ev/omin on the specified disk."""
    c,R=mp.mpf(c),mp.mpf(R); L=ell*(ell+1)
    aa,majorized=coefficient_data(domain,ell,N) if coeffs is None else coeffs
    if majorized:
        C=abs(N*(N+1)-L)*aa[N]+N*N*(aa[N-1] if N else 0)
        DD=(N+1)**2*aa[N]
        cc=dd=None
    else:
        cc=(N*(N+1)-L)*aa[N]-N*N*(aa[N-1] if N else 0)
        dd=-(N+1)**2*aa[N]
        C,DD=abs(cc),abs(dd)
    oldF=C*R**(-N-1)/(N+1)+DD*R**(-N-2)/(N+2)
    Ftriangle=C*power_tail(c,R,N+2)+DD*power_tail(c,R,N+3)
    F=Ftriangle
    cancelF=None; cancelC=None
    if not majorized:
        # c_N z+d_N=i c_N s+(c c_N+d_N). Convexity on 1/s∈[0,1/R].
        cancelC=max(abs(cc),abs(mp.j*cc+(c*cc+dd)/R))
        cancelF=cancelC*(c*c+R*R)**(-mp.mpf(N+1)/2)/(N+1)
        F=min(Ftriangle,cancelF)
    if level=='legacy': F=oldF
    elif level=='geometry': pass
    elif level!='improved': raise ValueError(level)
    V0=mp.pi*L/(2*c)+1/c**2
    VR=L*power_tail(c,R,2)+power_tail(c,R,3); Vseg=V0-VR
    Ks,Kt,mus,mut=domain.phase(c,R)
    if level in ('legacy','geometry'): Ks=Kt=1/domain.omin
    E0=Kt*F*mp.exp(Kt*VR)
    E1=F*mp.exp(Kt*VR)
    Fw=Vw=None
    if level=='improved' and refine:
        Fw=C*weighted_power_upper(c,R,N+2,mut)+DD*weighted_power_upper(c,R,N+3,mut)
        if cancelC is not None:
            Fw=min(Fw,cancelC*weighted_power_upper(c,R,N+2,mut))
        Fw=min(F,Fw)
        Vw=L*weighted_power_upper(c,R,2,mut)+weighted_power_upper(c,R,3,mut)
        E0=min(Kt*F,(F+Fw)/(2*domain.omin))*mp.exp(Kt*VR)
        E1=min(E1,Fw+E0*Vw)
    A=E0+Ks*E1
    Em=A*mp.exp(Ks*Vseg)
    Ev=E1+A/Ks*mp.expm1(Ks*Vseg)
    if level=='improved' and refine:
        # Uniform bound on the exponentially weighted segment potential integral.
        Vweighted=min(Vseg,(L/c**2+1/c**3)*(-mp.expm1(-mus*R))/mus)
        Ev=min(Ev,E1*mp.exp(-mus*R)+Em*Vweighted)
    return {'state_norm':Em+Ev/domain.omin,'matching_value':Em,'matching_derivative':Ev,
            'terminal_value':E0,'terminal_derivative':E1,'residual_integral_bound':F,
            'legacy_residual_bound':oldF,'geometry_triangle_residual':Ftriangle,
            'combined_residual_bound':cancelF,'weighted_residual_bound':Fw,
            'segment_kernel_bound':Ks,'tail_kernel_bound':Kt,
            'phase_decay_segment':mus,'phase_decay_tail':mut,
            'uniform_disk':bool(domain.radius),'level':level,'weighted_refinement':refine}


def optimize(domain,ell,c,target,level='improved'):
    """First feasible radius in five-unit grid; best order at each radius."""
    records=[]
    for R in range(20,156,5):
        nmax=max(80,int(mp.ceil(4*domain.omin*mp.sqrt(mp.mpf(c)**2+R*R)))+20)
        aa=coefficient_data(domain,ell,nmax)
        vals=[]
        for N in range(6,nmax+1):
            bd=tail_bound(domain,ell,c,R,N,aa,level,refine=False)
            vals.append((bd['state_norm'],N,bd))
        best=min(vals,key=lambda a:a[0])
        records.append({'radius':R,'best_order':best[1],'best_state_bound':best[0]})
        feasible=[x for x in vals if x[0]<=target]
        if feasible:
            # Minimize terminal polynomial cost subject to the common target.
            chosen=min(feasible,key=lambda a:a[1])
            return {'radius':R,'order':chosen[1],'best_order_at_radius':best[1],
                    'uniform_bound':tail_bound(domain,ell,c,R,chosen[1],aa,level,refine=True),
                    'search':records}
    raise RuntimeError('No feasible radius <=155 for target')



def global_residual_bounds(domain,ell,c,R,N):
    """Single full-ray Volterra estimate for the C1 joined approximation."""
    c=mp.mpf(c)
    local=tail_bound(domain,ell,c,R,N,level='improved',refine=False)
    F=local['residual_integral_bound']; k=local['segment_kernel_bound']
    V0=mp.pi*ell*(ell+1)/(2*c)+1/c**2
    H=F*mp.exp(k*V0)
    return {'matching_value':k*H,'matching_derivative':H,
            'state_norm':k*H+H/domain.omin,'resolvent_residual_integral':H,
            'full_ray_kernel':k,'residual_integral_bound':F}


def direct_amplitude_bounds(omega,ell,c,R,N):
    """Pointwise horizon factors are computed numerically, not certified."""
    domain=Domain(omega); g=global_residual_bounds(domain,ell,c,R,N)
    hh,hx,_=contour.horizon_scaled_state(omega,ell,str(c),180,48,'.20')
    H=g['resolvent_residual_integral'];k=g['full_ray_kernel']
    direct=(abs(hh)+k*abs(hx))*H/(2*abs(omega))
    by_state=(abs(hh)*g['matching_derivative']+
              abs(2*mp.j*omega*hh-hx)*g['matching_value'])/(2*abs(omega))
    return {'global_state':g,'direct_kernel_identity_amplitude':direct,
            'global_state_amplitude':by_state,'best_of_global_forms':min(direct,by_state),
            'abs_horizon_scaled':abs(hh),'abs_horizon_scaled_derivative':abs(hx),
            'state_to_amplitude_norm_amplification':
                max(abs(hh)*abs(omega),abs(2*mp.j*omega*hh-hx))/(2*abs(omega))}


def amplitude_bound(state_bound,omega,ell,c):
    # Horizon values are numerical here, so this is a diagnostic evaluation of
    # the analytic linear functional, not a uniform disk amplitude certificate.
    hh,hx,_=contour.horizon_scaled_state(omega,ell,str(c),180,48,'.20')
    return (abs(hh)*state_bound['matching_derivative']+
            abs(2*mp.j*omega*hh-hx)*state_bound['matching_value'])/(2*abs(omega))


def dump(path,data):
    Path(path).write_text(json.dumps(contour.serialize(data),indent=2)+'\n')


def campaign(output='analysis/qnm_tail_optimization_results.json',dps=40,target='1e-25'):
    mp.mp.dps=dps; target=mp.mpf(target)
    cf=json.loads(Path(__file__).with_name('qnm_leaver_results.json').read_text())
    result={'target_uniform_outgoing_state_tail_norm':target,
            'state_norm':'|delta m|+|delta m_x|/omega_min at the matching anchor',
            'working_dps':dps,'certified_arithmetic':False,
            'qualification':'Exact-arithmetic tail estimates only; no Taylor or roundoff certification and no frequency-error certificate.',
            'optimization_grid':'R=20,25,...; N>=6 scanned through at least four |Omega_min| sqrt(c^2+R^2)+20',
            'modes':[]}
    cases=[(0,0,('.22','-.21'),'1.5'),(2,0,('.97','-.19'),'1.5'),(0,2,('.15','-1.20'),'5')]
    for ell,n,seed,c in cases:
        start=time.perf_counter(); dom=Domain(mp.mpc(*seed),'.01')
        ref=next(x for x in cf['modes'] if x['ell']==ell and x['overtone']==n)
        wr=mp.mpc(ref['omega']['real'],ref['omega']['imag']); point=Domain(wr)
        record={'ell':ell,'n':n,'frequency_disk_center':dom.center,'frequency_disk_radius':dom.radius,'anchor':c}
        record['selection_legacy']=optimize(dom,ell,c,target,'legacy')
        print('selected legacy',ell,n,record['selection_legacy']['radius'],record['selection_legacy']['order'],flush=True)
        record['selection_improved']=optimize(dom,ell,c,target,'improved')
        sel=record['selection_improved']; R,N=sel['radius'],sel['order']
        print('selected improved',ell,n,R,N,'norm=',mp.nstr(sel['uniform_bound']['state_norm'],8),flush=True)
        # Compare sharpening levels at exactly the same selected endpoint settings.
        record['point_bound_comparison']={lvl:tail_bound(point,ell,c,R,N,level=lvl,refine=True) for lvl in ('legacy','geometry','improved')}
        common={'anchor':c,'taylor_order':48,'step_fraction':'.20','horizon_terms':180,'scaled_horizon':True}
        evals={}
        for tag,r,nn in [('previous_fixed',140,65),('legacy_adaptive',record['selection_legacy']['radius'],record['selection_legacy']['order']),('adaptive_fixed_domain',R,N)]:
            settings={**common,'radius':r,'asympt_order':nn}
            timings=[]
            for _ in range(2):
                st=time.perf_counter(); A,diag=contour.incoming_amplitude(wr,ell,diagnostics=True,**settings)
                timings.append(time.perf_counter()-st)
            bd=tail_bound(point,ell,c,r,nn,level='improved',refine=True)
            evals[tag]={'settings':settings,'A_at_CF':A,'abs_A_at_CF':abs(A),'wall_seconds':timings,
                        'outgoing_steps':diag['steps'],'horizon_steps':diag['horizon_steps'],
                        'improved_tail_bound':bd,'amplitude_tail_bound_diagnostic':amplitude_bound(bd,wr,ell,c)}
            print('evaluation',ell,n,tag,'steps',diag['steps'],'seconds',timings,'|A|',mp.nstr(abs(A),8),flush=True)
        record['evaluations']=evals
        record['evaluation_speedup']=sum(evals['previous_fixed']['wall_seconds'])/sum(evals['adaptive_fixed_domain']['wall_seconds'])
        record['amplitude_difference_adaptive_minus_fixed']=abs(evals['previous_fixed']['A_at_CF']-evals['adaptive_fixed_domain']['A_at_CF'])
        result['modes'].append(record);dump(output,result)
        # A genuinely independent rounded-seed root search with endpoints fixed on disk.
        rootruns={}
        for tag in ('previous_fixed','legacy_adaptive','adaptive_fixed_domain'):
            settings=evals[tag]['settings'];st=time.perf_counter()
            root,calls=contour.solve_root(ell,seed,settings,dps)
            elapsed=time.perf_counter()-st
            A,diag=contour.incoming_amplitude(root,ell,diagnostics=True,**settings)
            rootruns[tag]={'root':root,'calls':calls,'wall_seconds':elapsed,'abs_A':abs(A),
                           'distance_from_disk_center':abs(root-dom.center),'difference_from_cf':abs(root-wr)}
            print('root completed',ell,n,tag,'seconds',elapsed,'CFdiff',mp.nstr(abs(root-wr),8),flush=True)
            record['root_searches']=rootruns;dump(output,result)
        record['root_search_speedup']=rootruns['previous_fixed']['wall_seconds']/rootruns['adaptive_fixed_domain']['wall_seconds']
        record['difference_between_contour_roots']=abs(rootruns['previous_fixed']['root']-rootruns['adaptive_fixed_domain']['root'])
        record['elapsed_seconds']=time.perf_counter()-start
        dump(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',default='analysis/qnm_tail_optimization_results.json')
    p.add_argument('--dps',type=int,default=40)
    p.add_argument('--target',default='1e-25')
    args=p.parse_args();campaign(args.output,args.dps,args.target)
