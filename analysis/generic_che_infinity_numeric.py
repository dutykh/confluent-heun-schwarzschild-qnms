#!/usr/bin/env python3
"""Generic CHE exponential connection amplitude on specified recessive charts.

Only mpmath and standard-library routines are used for the infinity transport.
Finite inverse-power terminal data are transported by arbitrary-precision local
Taylor series along a radial ray, a circle of radius>1, and an inward radial
segment. Formal-series convergence is not assumed. Analytic terminal-tail bounds
are reported separately from empirical ODE/precision checks.

The path's starting log z uses angle(c) in (0,2pi], assigning 2pi on the positive real ray; the circular arc reaches
pi/4 without crossing z=0 or1. At z*=.5+.5i, log z and log(z-1) are principal.
Changing this continuation path can change the connection data.
"""
from __future__ import annotations
import json,time,argparse
from pathlib import Path
import mpmath as mp
from generic_che_path_verify import local_coefficients,evaluate_series,ordinary_basis


def dictionary(p):
    alpha,beta,gamma,delta,eta=map(mp.mpc,p)
    A=delta+alpha*(beta+gamma+2)/2
    B=eta+beta/2+(gamma-alpha)*(beta+1)/2
    if alpha==0:raise ValueError('alpha=0 is outside this infinity chart')
    rho=-A/alpha;u=beta+1+2*rho;v=gamma+1
    C=rho*(rho+beta);D=A+B+rho*(gamma+1)
    return alpha,beta,gamma,A,B,rho,u,v,C,D


def asymptotic_coefficients(p,N):
    alpha,_,_,_,_,_,u,v,C,D=dictionary(p)
    aa=[mp.mpc(1)]
    for n in range(N):
        am=aa[n-1] if n else 0
        aa.append(((n*(n+1-u-v+alpha)+C+D)*aa[n]-((n-1)*(n-u)+C)*am)/(alpha*(n+1)))
    return aa


def polynomial_state(aa,z):
    return (mp.fsum(aa[n]*z**(-n) for n in range(len(aa))),
            mp.fsum(-n*aa[n]*z**(-n-1) for n in range(1,len(aa))))


def radial_tail_bounds(p,R,N,ray=None):
    alpha,_,_,_,_,_,u,v,C,D=dictionary(p)
    aa=asymptotic_coefficients(p,N)
    EN=(N*(N+1-u-v+alpha)+C+D)*aa[N]-((N-1)*(N-u)+C)*(aa[N-1] if N else 0)
    FN=-(N*(N+1-u)+C)*aa[N]
    c=-mp.conj(alpha)/abs(alpha) if ray is None else mp.mpc(ray)/abs(mp.mpc(ray))
    sigma=-mp.re(alpha*c)-abs(u)/R-abs(v)/(R-1)
    if sigma<=0:raise ValueError('Terminal ray monotonicity bound is not positive')
    F=R/(R-1)*(abs(EN)*R**(-N-1)/(N+1)+abs(FN)*R**(-N-2)/(N+2))
    Bstar=abs(C+D)/R+abs(D)*(-mp.log1p(-1/R)-1/R)
    error=F/sigma*mp.exp(Bstar/sigma)
    return {'sigma_on_terminal_tail':sigma,'residual_integral_majorant':F,
            'potential_tail_majorant':Bstar,'terminal_m_error_bound':error,
            'terminal_mz_error_bound':sigma*error,
            'qualification':'Analytic terminal truncation only; subsequent floating arithmetic and Taylor transport are not interval certified.'}


def taylor_step(z,m,w,h,pars,order):
    alpha,u,v,C,D=pars[:5]
    extra=pars[5] if len(pars)>5 else 0
    invz=1/z;invm=1/(z-1);az=invz;am=invm
    pp=[];qq=[]
    for k in range(order):
        pp.append(u*az+v*am+(alpha if k==0 else 0))
        qq.append(C*(k+1)*az*invz+D*(am-az)+extra*(k+1)*am*invm)
        az*=-invz;am*=-invm
    mc=[m];wc=[w]
    for n in range(order):
        mc.append(wc[n]/(n+1))
        wc.append(-mp.fsum(qq[k]*mc[n-k]+pp[k]*wc[n-k] for k in range(n+1))/(n+1))
    mn=mc[-1];wn=wc[-1]
    for k in range(order-1,-1,-1):mn=mn*h+mc[k];wn=wn*h+wc[k]
    return mn,wn


def line_transport(z,target,state,pars,order,step_fraction):
    m,w=state;steps=0
    while z!=target:
        remainder=target-z;distance=abs(remainder)
        length=min(distance,step_fraction*min(abs(z),abs(z-1)),1/abs(pars[0]))
        h=remainder/distance*length
        m,w=taylor_step(z,m,w,h,pars,order)
        z=target if length==distance else z+h
        steps+=1
    return (m,w),steps


def infinity_state(p,R=64,N=32,taylor_order=48,step_fraction='.2',circle_radius=2,ray=None):
    data=dictionary(p);alpha,rho,u,v,C,D=data[0],data[5],data[6],data[7],data[8],data[9]
    pars=(alpha,u,v,C,D)
    R=mp.mpf(R);radius=mp.mpf(circle_radius);step_fraction=mp.mpf(step_fraction)
    if ray is None:c=-mp.conj(alpha)/abs(alpha)
    else:c=mp.mpc(ray)/abs(mp.mpc(ray))
    if mp.re(alpha*c)>=0:raise ValueError('Ray does not select the recessive algebraic branch')
    if radius<=1 or R<=radius:raise ValueError('Need R>circle radius>1')
    theta=mp.arg(c)
    if theta<=0:theta+=2*mp.pi
    z=c*R;state=polynomial_state(asymptotic_coefficients(p,N),z)
    state,steps=line_transport(z,c*radius,state,pars,taylor_order,step_fraction)
    z=c*radius;target_theta=mp.pi/4
    # Chords approximate the stated circle by a homotopic polygon; no singularity
    # lies between each short chord and arc, and Taylor steps remain local.
    count=max(1,int(mp.ceil(abs(theta-target_theta)/mp.mpf('.12'))))
    for j in range(1,count+1):
        angle=theta+(target_theta-theta)*j/count
        target=radius*mp.exp(mp.j*angle)
        state,ns=line_transport(z,target,state,pars,taylor_order,step_fraction)
        z=target;steps+=ns
    match=mp.mpc('.5','.5')
    state,ns=line_transport(z,match,state,pars,taylor_order,step_fraction);steps+=ns
    return state,{'steps':steps,'ray':c,'starting_log_angle':theta,'matching_point':match,
                  'circle_radius':radius,'terminal_tail_bounds':radial_tail_bounds(p,R,N,ray=c)}


def algebraic_state_from_m(p,z,mstate):
    rho=dictionary(p)[5];m,w=mstate
    power=mp.exp(rho*mp.log(z))
    return power*m,power*(w+rho*m/z)


def coefficient_from_states(p,z,hstate,mstate):
    alpha,beta,gamma,*_=dictionary(p)
    y,yp=algebraic_state_from_m(p,z,mstate)
    h,hp=hstate
    return mp.exp(alpha*z+(beta+1)*mp.log(z)+(gamma+1)*mp.log(z-1))*(h*yp-hp*y)/alpha


def coefficient(p,diagnostics=False,**settings):
    state,diag=infinity_state(p,**settings)
    z=diag['matching_point']
    h=evaluate_series(local_coefficients(p,800),z)
    value=coefficient_from_states(p,z,h,state)
    if diagnostics:
        data=dictionary(p);pars=(data[0],data[6],data[7],data[8],data[9])
        other=[]
        for target in [mp.mpc('.4','.3'),mp.mpc('.6','.2')]:
            ss,ns=line_transport(z,target,state,pars,settings.get('taylor_order',48),mp.mpf(settings.get('step_fraction','.2')))
            hh=evaluate_series(local_coefficients(p,800),target)
            vv=coefficient_from_states(p,target,hh,ss)
            other.append({'z':target,'coefficient':vv,'absolute_difference':abs(vv-value),'extra_steps':ns})
        diag.update({'m':state[0],'mz':state[1],'yalg':algebraic_state_from_m(p,z,state),
                     'matching_point_checks':other})
        return value,diag
    return value


def record(value):
    if isinstance(value,dict):return {k:record(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [record(v) for v in value]
    if isinstance(value,mp.mpc):return {'real':mp.nstr(value.real,45),'imag':mp.nstr(value.imag,45)}
    if isinstance(value,mp.mpf):return mp.nstr(value,25)
    return value


def asymptotic_residual_check(p,N=6):
    data=dictionary(p);alpha,u,v,C,D=data[0],data[6],data[7],data[8],data[9]
    aa=asymptotic_coefficients(p,N);z=mp.mpc('2.3','1.7')
    m,w=polynomial_state(aa,z)
    mpp=mp.fsum(n*(n+1)*aa[n]*z**(-n-2) for n in range(1,N+1))
    actual=mpp+(alpha+u/z+v/(z-1))*w+(C/z**2+D/(z*(z-1)))*m
    E=(N*(N+1-u-v+alpha)+C+D)*aa[N]-((N-1)*(N-u)+C)*aa[N-1]
    F=-(N*(N+1-u)+C)*aa[N]
    exact=(E*z+F)/(z**(N+2)*(z-1))
    return abs(actual-exact)


def generic_campaign(output='analysis/generic_che_infinity_numeric.json'):
    mp.mp.dps=75
    examples=[('complex_A',[( '-1.2','.3'),('.3','.2'),('-.4','.1'),('.2','-.15'),('.6','.05')]),
              ('complex_B',[('.8','.6'),('-.2','.15'),('.4','-.25'),('.35','.2'),('-.1','.4')])]
    configs=[dict(R=24,N=16,taylor_order=36,step_fraction='.22'),
             dict(R=40,N=24,taylor_order=44,step_fraction='.20'),
             dict(R=64,N=32,taylor_order=52,step_fraction='.18'),
             dict(R=96,N=48,taylor_order=64,step_fraction='.16')]
    out={'method':'finite inverse-power terminal polynomial plus arbitrary-precision Taylor continuation',
         'working_decimal_digits':75,'certified':False,
         'qualification':'Terminal-tail analytic bounds are separate from empirical transport/roundoff checks; connection path and branches explicitly fixed.',
         'examples':[]}
    for name,raw in examples:
        p=tuple(mp.mpc(*pair) for pair in raw);start=time.time()
        case={'name':name,'parameters':p,'asymptotic_residual_identity_error':asymptotic_residual_check(p),'runs':[]}
        previous=None
        for cfg in configs:
            val,diag=coefficient(p,diagnostics=True,**cfg)
            case['runs'].append({'settings':cfg,'coefficient':val,'change_from_previous':None if previous is None else abs(val-previous),**diag})
            previous=val
            print(name,cfg,'B=',mp.nstr(val,30),'delta=',mp.nstr(case['runs'][-1]['change_from_previous'],8),flush=True)
        # Same terminal data, finer transport and a different homotopic circular path.
        refined={**configs[-1],'taylor_order':76,'step_fraction':'.12','circle_radius':3}
        val,diag=coefficient(p,diagnostics=True,**refined)
        case['transport_path_check']={'settings':refined,'coefficient':val,'difference_from_last':abs(val-previous),**diag}
        case['elapsed_seconds']=time.time()-start
        out['examples'].append(case)
        Path(output).write_text(json.dumps(record(out),indent=2)+'\n')
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign',action='store_true')
    parser.add_argument('--output',default='analysis/generic_che_infinity_numeric.json')
    args=parser.parse_args()
    if args.campaign:generic_campaign(args.output)
    else:
        mp.mp.dps=60
        p=tuple(map(mp.mpc,['.7','-.7','.4','.2','.6']))
        value,diag=coefficient(p,diagnostics=True,R=32,N=20,taylor_order=44)
        print(json.dumps(record({'coefficient':value,**diag}),indent=2))
