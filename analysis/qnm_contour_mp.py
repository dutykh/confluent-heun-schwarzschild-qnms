#!/usr/bin/env python3
"""High-precision Taylor diagnostic for scalar Schwarzschild contour matching.

Independent implementation: no SciPy integrator or qnm_contour.py import is used.
Omega=2M*omega, exp(-i omega t). Propagate (m,v=m_x) from z=c+iR
vertically to c using local Taylor series for the rational first-order system.
The outgoing inverse-power polynomial is finite, not presumed convergent.
Analytic bounds concern exact terminal data and propagation, not floating-point
arithmetic or Taylor truncation. Order/step/precision comparisons are empirical.
Dependencies: mpmath; standard library.
"""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
import mpmath as mp


def asymptotic_coefficients(omega,ell,order):
    L=ell*(ell+1); a=[mp.mpc(1)]
    for k in range(order):
        prev=a[k-1] if k else 0
        a.append(((k*(k+1)-L)*a[k]-k*k*prev)/(2*mp.j*omega*(k+1)))
    return a


def terminal_state(omega,ell,anchor,radius,order):
    a=asymptotic_coefficients(omega,ell,order)
    z=anchor+mp.j*radius
    P=mp.fsum(a[n]*z**(-n) for n in range(order+1))
    Pz=mp.fsum(-n*a[n]*z**(-n-1) for n in range(1,order+1))
    return P,(z-1)/z*Pz


def analytic_tail_bounds(omega,ell,anchor,radius,order):
    """Analytic absolute terminal/matching error bounds; evaluated with mp."""
    p=omega.real; q=max(mp.mpf(0),-omega.imag)
    if not (p>0 and q*q<=4*p*p*anchor*(anchor-1)):
        raise ValueError('Vertical ray fails sufficient monotone phase condition')
    L=ell*(ell+1); a=asymptotic_coefficients(omega,ell,order); N=order
    c=(N*(N+1)-L)*a[N]-N*N*(a[N-1] if N else 0)
    d=-(N+1)**2*a[N]
    F=abs(c)*radius**(-N-1)/(N+1)+abs(d)*radius**(-N-2)/(N+2)
    V0=mp.pi*L/(2*anchor)+1/anchor**2
    sq=mp.sqrt(anchor*anchor+radius*radius)
    VR=L/anchor*mp.atan(anchor/radius)+1/(sq*(sq+radius))
    E0=F/abs(omega)*mp.exp(VR/abs(omega)); E1=abs(omega)*E0
    bseg=(V0-VR)/abs(omega); AA=E0+E1/abs(omega)
    return {'terminal_value':E0,'terminal_derivative':E1,
            'matching_value':AA*mp.exp(bseg),
            'matching_derivative':E1+abs(omega)*AA*mp.expm1(bseg),
            'residual_integral':F,'total_potential_integral_bound':V0,
            'phase_condition_margin':4*p*p*anchor*(anchor-1)-q*q}


def taylor_step_outgoing(z,m,v,h,omega,ell,order):
    """Taylor coefficients of m_z=b v, v_z=V m-2iOmega b v."""
    L=ell*(ell+1)
    iz=1/z; izm=1/(z-1)
    b=[z*izm]
    b.extend((-1)**k*izm**(k+1) for k in range(1,order))
    vv=[(-1)**k*((k+1)*L*iz**(k+2)+(k+1)*(k+2)*iz**(k+3)/2)
        for k in range(order)]
    mc=[m]; vc=[v]
    for n in range(order):
        bv=mp.fsum(b[k]*vc[n-k] for k in range(n+1))
        vm=mp.fsum(vv[k]*mc[n-k] for k in range(n+1))
        mc.append(bv/(n+1)); vc.append((vm-2*mp.j*omega*bv)/(n+1))
    mn=mc[-1]; vn=vc[-1]
    for j in range(order-1,-1,-1):
        mn=mn*h+mc[j]; vn=vn*h+vc[j]
    return mn,vn


def propagate_outgoing(omega,ell,anchor,radius,asympt_order=65,
                       taylor_order=40,step_fraction=mp.mpf('.25')):
    anchor=mp.mpf(anchor); radius=mp.mpf(radius); step_fraction=mp.mpf(step_fraction)
    m,v=terminal_state(omega,ell,anchor,radius,asympt_order)
    s=radius; steps=0
    while s>0:
        z=anchor+mp.j*s
        distance=min(step_fraction*min(abs(z),abs(z-1)),1/abs(omega),s)
        h=-mp.j*distance
        m,v=taylor_step_outgoing(z,m,v,h,omega,ell,taylor_order)
        s-=distance; steps+=1
    return m,v,steps


def horizon_state(omega,ell,anchor=mp.mpf('1.5'),terms=180):
    """Convergent horizon series, restricted here to |1-anchor|<1."""
    anchor=mp.mpf(anchor)
    if not 1<anchor<2:
        raise ValueError('This diagnostic uses direct horizon series only at 1<c<2')
    w=1-anchor; L=ell*(ell+1); lam=2*mp.j*omega
    h0=mp.mpc(1); hprev=mp.mpc(0); H=h0; Hz=mp.mpc(0)
    for n in range(1,terms+1):
        den=n*(n-lam)
        if not den: raise ValueError('Horizon resonance')
        hn=((n*(n-1)+4*omega**2-L)*h0-(lam*(n-1)+4*omega**2)*hprev)/den
        H+=hn*w**n; Hz-=n*hn*w**(n-1)
        hprev,h0=h0,hn
    f=(anchor-1)/anchor
    g=mp.exp(-2*mp.j*omega+mp.j*omega*anchor)*anchor*(anchor-1)**(-mp.j*omega)
    rg=mp.j*omega+1/anchor-mp.j*omega/(anchor-1)
    return g*H,g*f*(Hz+rg*H)



def horizon_scaled_state(omega,ell,anchor='1.5',terms=180,
                         taylor_order=48,step_fraction='.20'):
    """h=exp(iOmega*x)u_H and h_x; propagate with Omega -> -Omega."""
    z0=mp.mpf('1.5'); anchor=mp.mpf(anchor); w=1-z0
    L=ell*(ell+1); lam=2*mp.j*omega
    h0=mp.mpc(1); hprev=mp.mpc(0); H=h0; Hz=mp.mpc(0)
    for n in range(1,terms+1):
        den=n*(n-lam)
        if not den: raise ValueError('Horizon resonance')
        hn=((n*(n-1)+4*omega**2-L)*h0-(lam*(n-1)+4*omega**2)*hprev)/den
        H+=hn*w**n; Hz-=n*hn*w**(n-1)
        hprev,h0=h0,hn
    G=z0*mp.exp(2*mp.j*omega*(z0-1)); f=(z0-1)/z0
    h=G*H; hx=f*G*(Hz+(2*mp.j*omega+1/z0)*H)
    z=z0; steps=0
    while abs(anchor-z)>0:
        remaining=abs(anchor-z)
        distance=min(mp.mpf(step_fraction)*min(abs(z),abs(z-1)),1/abs(omega),remaining)
        step=(anchor-z)/remaining*distance
        h,hx=taylor_step_outgoing(z,h,hx,step,-omega,ell,taylor_order)
        z+=step; steps+=1
        if distance==remaining: break
    return h,hx,steps


def incoming_amplitude(omega,ell,anchor='1.5',radius=140,asympt_order=65,
                       taylor_order=40,step_fraction='.25',horizon_terms=180,
                       diagnostics=False,scaled_horizon=False):
    omega=mp.mpc(omega); anchor=mp.mpf(anchor); radius=mp.mpf(radius)
    m,mx,steps=propagate_outgoing(omega,ell,anchor,radius,asympt_order,taylor_order,step_fraction)
    x=anchor+mp.log(anchor-1); phase=mp.exp(mp.j*omega*x)
    if scaled_horizon:
        hh,hhx,hsteps=horizon_scaled_state(omega,ell,anchor,horizon_terms,taylor_order,step_fraction)
        uH=hh/phase; uHx=(hhx-mp.j*omega*hh)/phase
        D=2*mp.j*omega*hh*m+hh*mx-hhx*m
    else:
        uH,uHx=horizon_state(omega,ell,anchor,horizon_terms)
        hsteps=0
        D=phase*(uH*mx+(mp.j*omega*uH-uHx)*m)
    A=D/(2*mp.j*omega)
    if diagnostics:
        bounds=analytic_tail_bounds(omega,ell,anchor,radius,asympt_order)
        bounds['determinant_tail_only']=abs(phase)*(abs(uH)*bounds['matching_derivative']+abs(mp.j*omega*uH-uHx)*bounds['matching_value'])
        bounds['amplitude_tail_only']=bounds['determinant_tail_only']/(2*abs(omega))
        return A,{'tail_bounds':bounds,'steps':steps,'horizon_steps':hsteps,'m':m,'mx':mx,'uH':uH,'uHx':uHx,'determinant':D}
    return A


def complex_record(z,digits=35):
    return {'real':mp.nstr(z.real,digits),'imag':mp.nstr(z.imag,digits)}


def serialize(value):
    if isinstance(value,dict): return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,list): return [serialize(v) for v in value]
    if isinstance(value,mp.mpc): return complex_record(value)
    if isinstance(value,mp.mpf): return mp.nstr(value,18)
    return value


def solve_root(ell,seed,settings,dps=40):
    """Secant root from explicit rounded seeds; no continued-fraction data used."""
    with mp.workdps(dps):
        initial=mp.mpc(seed[0],seed[1]); calls=0
        def f(w):
            nonlocal calls
            calls+=1
            val=incoming_amplitude(w,ell,**settings)
            print(f'root l={ell} eval={calls} Omega={mp.nstr(w,22)} |A|={mp.nstr(abs(val),7)}',flush=True)
            return val
        root=mp.findroot(f,(initial,initial+mp.mpc('.001','.001')),
                         solver='secant',tol=mp.mpf('1e-29'),maxsteps=20,verify=True)
        return +root,calls


def campaign(output='analysis/qnm_contour_mp_results.json',dps=40):
    mp.mp.dps=dps
    cf=json.loads(Path(__file__).with_name('qnm_leaver_results.json').read_text())
    result={'conventions':'scalar Schwarzschild, Omega=2Momega, exp(-iomega t)',
            'arithmetic_dps':dps,'method':'independent mpmath local Taylor propagation on vertical ray',
            'certified':False,'qualification':'Analytic terminal-tail bounds only; Taylor truncation and floating-point errors are empirically assessed. Frequencies found from rounded seeds independently of CF.',
            'modes':[]}
    output=Path(output)
    for ell,seed in [(0,('.22','-.21')),(2,('.97','-.19'))]:
        start=time.time()
        ref=next(v for v in cf['modes'] if v['ell']==ell and v['overtone']==0)
        wr=mp.mpc(ref['omega']['real'],ref['omega']['imag'])
        mode={'ell':ell,'n':0,'reference_cf':complex_record(wr),'reference_diagnostics':[]}
        configs=[dict(radius=100,asympt_order=50,taylor_order=32,step_fraction='.25'),
                 dict(radius=140,asympt_order=65,taylor_order=40,step_fraction='.25'),
                 dict(radius=140,asympt_order=65,taylor_order=48,step_fraction='.20')]
        for cfg in configs:
            a,diag=incoming_amplitude(wr,ell,diagnostics=True,**cfg)
            mode['reference_diagnostics'].append({'settings':cfg,'A':a,'abs_A':abs(a),**diag})
            print('diagnostic',ell,cfg,'|A|=',mp.nstr(abs(a),12),'tail=',mp.nstr(diag['tail_bounds']['amplitude_tail_only'],12),'steps=',diag['steps'],flush=True)
            result['modes']=[*result['modes'],mode] if not result['modes'] or result['modes'][-1] is not mode else result['modes']
            output.write_text(json.dumps(serialize(result),indent=2)+'\n')
        settings=configs[-1]
        root,calls=solve_root(ell,seed,settings,dps)
        ar,diagr=incoming_amplitude(root,ell,diagnostics=True,**settings)
        mode.update({'root':root,'root_seed':seed,'root_calls':calls,'root_settings':settings,
                     'root_amplitude':ar,'root_amplitude_abs':abs(ar),
                     'difference_from_cf':abs(root-wr),'root_diagnostics':diagr})
        # Arithmetic check at a fixed root, without importing any double-precision solver.
        with mp.workdps(dps+15):
            ahi,dhi=incoming_amplitude(root,ell,diagnostics=True,**settings)
        mode['precision_check']={'dps':dps+15,'A':ahi,'amplitude_change':abs(ahi-ar)}
        mode['elapsed_seconds']=round(time.time()-start,3)
        output.write_text(json.dumps(serialize(result),indent=2)+'\n')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--campaign',action='store_true')
    ap.add_argument('--output',default='analysis/qnm_contour_mp_results.json')
    ap.add_argument('--dps',default=40,type=int)
    args=ap.parse_args()
    mp.mp.dps=args.dps
    if args.campaign: campaign(args.output,args.dps)
    else:
        a,d=incoming_amplitude(mp.mpc('.2209098781608393717509230123','-.2097914341737619175634781337'),0,diagnostics=True)
        print(json.dumps(serialize({'A':a,**d}),indent=2))
