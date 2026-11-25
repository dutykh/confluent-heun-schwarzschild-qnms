#!/usr/bin/env python3
"""Incoming scalar Schwarzschild coefficient by forward horizon transport only.

No outgoing backward ODE is solved. With h=exp(iOmega*x)u_H, propagate (h,h_x)
up z=c+i s using the unchanged rational Taylor engine with Omega -> -Omega.
At z_R, A_N=[(2iOmega h-h_x)P_N+h f P_N']/2iOmega.
Exact-arithmetic identity: A_in-A_N=integral_R^infinity h r_N dz/(2iOmega).
The reported tail bound is conditional on exact h,h_x at R. The computed states
have uncertified Taylor/roundoff errors, so no certified QNM claim is made.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import mpmath as mp
import qnm_contour_mp as contour
import qnm_tail_optimization as bounds


def forward_horizon(omega,ell,anchor='1.5',radius=30,taylor_order=48,
                    step_fraction='.20',horizon_terms=180):
    omega=mp.mpc(omega);c=mp.mpf(anchor);R=mp.mpf(radius)
    h,hx,hsteps=contour.horizon_scaled_state(omega,ell,c,horizon_terms,taylor_order,step_fraction)
    s=mp.mpf(0);steps=0
    while s<R:
        z=c+mp.j*s
        ds=min(mp.mpf(step_fraction)*min(abs(z),abs(z-1)),1/abs(omega),R-s)
        h,hx=contour.taylor_step_outgoing(z,h,hx,mp.j*ds,-omega,ell,taylor_order)
        s+=ds;steps+=1
    return h,hx,steps,hsteps


def forward_coefficient(omega,ell,anchor='1.5',radius=30,asympt_order=60,
                        taylor_order=48,step_fraction='.20',horizon_terms=180,
                        diagnostics=False):
    omega=mp.mpc(omega);c=mp.mpf(anchor);R=mp.mpf(radius)
    h,hx,steps,hsteps=forward_horizon(omega,ell,c,R,taylor_order,step_fraction,horizon_terms)
    P,Px=contour.terminal_state(omega,ell,c,R,asympt_order)
    A=((2*mp.j*omega*h-hx)*P+h*Px)/(2*mp.j*omega)
    if diagnostics:
        domain=bounds.Domain(omega)
        b=bounds.tail_bound(domain,ell,c,R,asympt_order,level='improved',refine=False)
        F=b['residual_integral_bound'];kt=b['tail_kernel_bound']
        VR=ell*(ell+1)*bounds.power_tail(c,R,2)+bounds.power_tail(c,R,3)
        hs=(abs(h)+kt*abs(hx))*mp.exp(kt*VR)
        tail=F*hs/(2*abs(omega))
        return A,{'h_R':h,'hx_R':hx,'forward_steps':steps,'horizon_steps':hsteps,
                  'residual_integral_bound':F,'tail_kernel_bound':kt,
                  'tail_potential_bound':VR,'h_tail_sup_bound_conditional':hs,
                  'amplitude_tail_bound_conditional':tail}
    return A



def coefficient_tail_comparison(omega,ell,c,R,N,h,hx):
    """Several valid pointwise coefficient bounds at IDENTICAL endpoints."""
    omega=mp.mpc(omega);c=mp.mpf(c);R=mp.mpf(R)
    dom=bounds.Domain(omega)
    bb=bounds.tail_bound(dom,ell,c,R,N,level='improved',refine=False)
    F=bb['residual_integral_bound'];k=bb['tail_kernel_bound'];mu=bb['phase_decay_tail']
    VR=ell*(ell+1)*bounds.power_tail(c,R,2)+bounds.power_tail(c,R,3)
    aa=contour.asymptotic_coefficients(omega,ell,N)
    C=(N*(N+1)-ell*(ell+1))*aa[N]-N*N*aa[N-1]
    D=-(N+1)**2*aa[N]
    Fw=abs(C)*bounds.weighted_power_upper(c,R,N+2,mu)+abs(D)*bounds.weighted_power_upper(c,R,N+3,mu)
    combined=max(abs(C),abs(mp.j*C+(c*C+D)/R))
    Fw=min(F,Fw,combined*bounds.weighted_power_upper(c,R,N+2,mu))
    common=F*mp.exp(k*VR)/(2*abs(omega))
    direct=common*(abs(h)+k*abs(hx))
    state=common*(abs(h)+k*abs(2*mp.j*omega*h-hx))
    a=h-hx/(2*mp.j*omega);b=hx/(2*mp.j*omega)
    free=common*(abs(a)+abs(b))
    weighted=(abs(a)*F+abs(b)*Fw+(abs(a)+abs(b))*mp.expm1(k*VR)*F)/(2*abs(omega))
    return {'conventional_terminal_state_projection':state,
            'direct_forward_sup_bound':direct,'free_components_sup_bound':free,
            'weighted_free_components_bound':weighted,'best_bound':min(state,direct,free,weighted),
            'state_over_direct_ratio':state/direct,'state_over_weighted_ratio':state/weighted,
            'weighted_residual_over_unweighted':Fw/F,
            'abs_free_constant_component':abs(a),'abs_free_exponential_component':abs(b),
            'conditional_on_exact_horizon_state_at_R':True}


def best_order_on_disk(domain,ell,c,R):
    """Choose a fixed N minimizing the uniform residual bound in a scan."""
    nmax=max(80,int(mp.ceil(4*domain.omin*mp.sqrt(mp.mpf(c)**2+R*R)))+20)
    aa=bounds.coefficient_data(domain,ell,nmax)
    records=[]
    for N in range(6,nmax+1):
        bd=bounds.tail_bound(domain,ell,c,R,N,aa,level='improved',refine=False)
        records.append((bd['residual_integral_bound'],N))
    val,N=min(records)
    return N,{'scan_min_order':6,'scan_max_order':nmax,'uniform_residual_bound':val}


def root_search(ell,seed,settings,domain,dps=40):
    calls=0;trace=[]
    with mp.workdps(dps):
        initial=mp.mpc(*seed)
        def fn(w):
            nonlocal calls
            calls+=1
            if abs(w-domain.center)>domain.radius:
                raise ValueError('Root iterate left chosen frequency disk')
            A=forward_coefficient(w,ell,**settings)
            trace.append({'omega':w,'abs_A_N':abs(A)})
            return A
        w=mp.findroot(fn,(initial,initial+mp.mpc('.001','.001')),solver='secant',
                      tol=mp.mpf('1e-29'),maxsteps=20,verify=True)
    return w,calls,trace


def campaign(output='analysis/qnm_forward_coefficient_results.json',dps=40,target='1e-30'):
    mp.mp.dps=dps;target=mp.mpf(target)
    directory=Path(__file__).parent
    cf=json.loads((directory/'qnm_leaver_results.json').read_text())
    old=json.loads((directory/'qnm_tail_optimization_results.json').read_text())
    result={'frequency':'Omega=2Momega','method':'forward horizon transport only; finite residual identity for incoming coefficient',
            'working_dps':dps,'taylor_order':48,'step_fraction':'.20','horizon_terms':180,
            'root_seed_policy':'Rounded seed and shifted rounded seed; no CF frequency used to seed or select roots.',
            'conditional_amplitude_tail_target':target,
            'qualification':'Tail bounds condition on exact transported horizon states. Numerical transport, arithmetic, and frequency errors are not certified.',
            'modes':[]}
    cases=[(0,0,('.22','-.21'),'1.5'),(2,0,('.97','-.19'),'1.5'),(0,2,('.15','-1.20'),'5')]
    for ell,n,seed,c in cases:
        start=time.perf_counter();dom=bounds.Domain(mp.mpc(*seed),'.01')
        ref=next(x for x in cf['modes'] if x['ell']==ell and x['overtone']==n)
        wr=mp.mpc(ref['omega']['real'],ref['omega']['imag'])
        rec={'ell':ell,'n':n,'anchor':c,'seed':seed,'disk_radius':'.01','CF_frequency':wr,'radius_diagnostics':[],'root_searches':[]}
        settings_by_radius={}
        common={'anchor':c,'taylor_order':48,'step_fraction':'.20','horizon_terms':180}
        # Fixed-R diagnostics requested, without using CF values in root searches.
        for R in (20,30,40):
            N,order_info=best_order_on_disk(dom,ell,c,R)
            cfg={**common,'radius':R,'asympt_order':N};settings_by_radius[R]=(cfg,order_info)
            timings=[]
            for _ in range(2):
                st=time.perf_counter();A,diag=forward_coefficient(wr,ell,diagnostics=True,**cfg)
                timings.append(time.perf_counter()-st)
            rec['radius_diagnostics'].append({'settings':cfg,'order_selection':order_info,'A_at_CF':A,'abs_A_at_CF':abs(A),'wall_seconds':timings,**diag})
            print('forward diagnostic',ell,n,'R,N',R,N,'|A|',mp.nstr(abs(A),8),'tail',mp.nstr(diag['amplitude_tail_bound_conditional'],8),'steps',diag['forward_steps'],flush=True)
        result['modes'].append(rec);bounds.dump(output,result)
        # Independent root searches, fixed settings on each disk. Increase R only
        # between completed searches based on the root's a posteriori tail bound.
        for R in (20,30,40,50,60,80):
            if R not in settings_by_radius:
                N,order_info=best_order_on_disk(dom,ell,c,R)
                settings_by_radius[R]=({**common,'radius':R,'asympt_order':N},order_info)
            cfg,order_info=settings_by_radius[R]
            st=time.perf_counter();w,calls,trace=root_search(ell,seed,cfg,dom,dps);elapsed=time.perf_counter()-st
            A,diag=forward_coefficient(w,ell,diagnostics=True,**cfg)
            rr={'settings':cfg,'order_selection':order_info,'root':w,'calls':calls,'trace':trace,
                'wall_seconds':elapsed,'abs_A_N':abs(A),'difference_from_CF':abs(w-wr),**diag}
            rec['root_searches'].append(rr)
            print('forward root',ell,n,'R',R,'seconds',elapsed,'CFdiff',mp.nstr(abs(w-wr),9),'tail',mp.nstr(diag['amplitude_tail_bound_conditional'],9),flush=True)
            bounds.dump(output,result)
            if diag['amplitude_tail_bound_conditional']<=target:
                rec['selected_root_run_index']=len(rec['root_searches'])-1
                break
        if 'selected_root_run_index' not in rec:
            rec['target_reached']=False
        else:
            rec['target_reached']=True
            chosen=rec['root_searches'][rec['selected_root_run_index']]
            oldrec=next(x for x in old['modes'] if x['ell']==ell and x['n']==n)
            rec['backward_comparison']={tag:{'settings':oldrec['evaluations'][tag]['settings'],
                'root_wall_seconds':v['wall_seconds'],'root_difference_from_CF':v['difference_from_cf'],
                'root_amplitude_tail_bound_diagnostic':v['amplitude_tail_bound_at_root_diagnostic']}
                for tag,v in oldrec['root_searches'].items()}
            with mp.workdps(dps+15):
                Ap,Dp=forward_coefficient(chosen['root'],ell,diagnostics=True,**chosen['settings'])
            rec['precision_check']={'dps':dps+15,'A_N':Ap,'amplitude_change':abs(Ap-forward_coefficient(chosen['root'],ell,**chosen['settings']))}
            # Fixed-root changed order/step check; tail N remains the same.
            with mp.workdps(dps+15):
                fine={**chosen['settings'],'taylor_order':56,'step_fraction':'.18'}
                Af,Df=forward_coefficient(chosen['root'],ell,diagnostics=True,**fine)
            rec['transport_refinement_check']={'dps':dps+15,'settings':fine,'A_N':Af,'abs_A_N':abs(Af),**Df}
        rec['elapsed_seconds']=time.perf_counter()-start
        bounds.dump(output,result)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',default='analysis/qnm_forward_coefficient_results.json')
    ap.add_argument('--dps',default=40,type=int)
    ap.add_argument('--target',default='1e-30')
    args=ap.parse_args();campaign(args.output,args.dps,args.target)
