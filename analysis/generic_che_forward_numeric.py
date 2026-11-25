#!/usr/bin/env python3
"""Forward transport of the regular CHE solution in its bounded infinity gauge.

q=exp(alpha*z)*z^(u-rho)*(z-1)^v*h0 satisfies
q''-p*q'+(V-p')q=0. No outgoing/recessive solution is numerically propagated.
The same finite asymptotic polynomial PN yields the coefficient approximant
B_N=((p*PN+PN')q-PN*q')/alpha at the far endpoint. Its exact analytic tail error
is bounded separately; numerical state uncertainty is not interval certified.
"""
import json,time
from pathlib import Path
import mpmath as mp
from generic_che_infinity_numeric import dictionary,asymptotic_coefficients,polynomial_state,line_transport,radial_tail_bounds,record
from generic_che_path_verify import local_coefficients,evaluate_series


def forward_coefficient(p,R=64,N=32,taylor_order=52,step_fraction='.18',circle_radius=2,ray=None):
    data=dictionary(p)
    alpha,beta,gamma,A,B,rho,u,v,C,D=data
    R=mp.mpf(R);radius=mp.mpf(circle_radius);fraction=mp.mpf(step_fraction)
    c=-mp.conj(alpha)/abs(alpha) if ray is None else mp.mpc(ray)/abs(mp.mpc(ray))
    if mp.re(alpha*c)>=0:raise ValueError('Ray must be recessive for the algebraic branch')
    theta=mp.arg(c)
    if theta<=0:theta+=2*mp.pi
    z=mp.mpc('.5','.5')
    h,hp=evaluate_series(local_coefficients(p,800),z)
    factor=mp.exp(alpha*z+(u-rho)*mp.log(z)+v*mp.log(z-1))
    state=(factor*h,factor*(hp+(alpha+(u-rho)/z+v/(z-1))*h))
    # q''+pq*q'+Vq*q=0, Vq has an additional inverse square at z=1.
    qpars=(-alpha,-u,-v,C+u,D,v)
    target=radius*mp.exp(mp.j*mp.pi/4)
    state,steps=line_transport(z,target,state,qpars,taylor_order,fraction);z=target
    count=max(1,int(mp.ceil(abs(theta-mp.pi/4)/mp.mpf('.12'))))
    for j in range(1,count+1):
        angle=mp.pi/4+(theta-mp.pi/4)*j/count
        target=radius*mp.exp(mp.j*angle)
        state,ns=line_transport(z,target,state,qpars,taylor_order,fraction);steps+=ns;z=target
    target=c*R
    state,ns=line_transport(z,target,state,qpars,taylor_order,fraction);steps+=ns;z=target
    q,qp=state
    PN,PNp=polynomial_state(asymptotic_coefficients(p,N),z)
    pz=alpha+u/z+v/(z-1)
    BN=((pz*PN+PNp)*q-PN*qp)/alpha
    tail=radial_tail_bounds(p,R,N,ray=c)
    sigma=tail['sigma_on_terminal_tail']
    W=abs(C+u)/R+abs(D)*mp.log(R/(R-1))+abs(v)/(R-1)
    qmajorant=(abs(q)+abs(qp)/sigma)*mp.exp(W/sigma)
    error=tail['residual_integral_majorant']/abs(alpha)*qmajorant
    return BN,{'q_at_R':q,'qprime_at_R':qp,'ray':c,'starting_log_angle':theta,
               'steps':steps,'sigma':sigma,'transformed_potential_tail_majorant':W,
               'q_tail_majorant_from_endpoint':qmajorant,
               'coefficient_analytic_tail_bound':error,
               'qualification':'Tail bound assumes exact q(R),qprime(R); computed states and Taylor propagation are not interval certified.'}


def precision_transport_check(p,baseline=None):
    config=dict(R=128,N=64,taylor_order=80,step_fraction='.14')
    with mp.workdps(75):
        if baseline is None:baseline,_=forward_coefficient(p,**config)
    with mp.workdps(95):
        higher={**config,'taylor_order':92,'step_fraction':'.12','circle_radius':3}
        value,_=forward_coefficient(p,**higher)
        change=abs(value-baseline)
    return {'baseline_dps':75,'refined_dps':95,'refined_settings':higher,
            'coefficient_change':change,'qualification':'Combined precision, Taylor-order/step, and homotopic-path test at fixed terminal R,N.'}


def main():
    mp.mp.dps=75
    source=json.loads(Path('analysis/generic_che_infinity_numeric.json').read_text())
    out={'method':'forward regular-solution transport in bounded q gauge; no infinity solution propagation',
         'working_decimal_digits':75,'cases':[],'qualification':'Forward and backward approximants share PN and agree algebraically with exact transport. This checks independent equations/directions and normalization, not an independent asymptotic model.'}
    for case in source['examples']:
        p=tuple(mp.mpc(x['real'],x['imag']) for x in case['parameters'])
        row={'name':case['name'],'parameters':p,'runs':[]}
        for back in [case['runs'][1],case['runs'][-1]]:
            cfg=back['settings'];bv=mp.mpc(back['coefficient']['real'],back['coefficient']['imag'])
            start=time.time();fv,diag=forward_coefficient(p,**cfg)
            rr={'settings':cfg,'forward_coefficient':fv,'backward_coefficient':bv,
                'difference':abs(fv-bv),'elapsed_seconds':time.time()-start,**diag}
            row['runs'].append(rr)
            print(case['name'],cfg,'diff',mp.nstr(abs(fv-bv),10),'tail',mp.nstr(diag['coefficient_analytic_tail_bound'],10),flush=True)
        finer=dict(R=128,N=64,taylor_order=80,step_fraction='.14')
        fvhi,dh=forward_coefficient(p,**finer)
        row['forward_truncation_refinement']={'settings':finer,'coefficient':fvhi,'difference_from_R96_N48':abs(fvhi-fv),
            'difference_divided_by_previous_tail_bound':abs(fvhi-fv)/diag['coefficient_analytic_tail_bound'],**dh}
        print(case['name'],'forward refinement',mp.nstr(fvhi,38),'delta',mp.nstr(abs(fvhi-fv),10),'newtail',mp.nstr(dh['coefficient_analytic_tail_bound'],10),flush=True)
        row['precision_transport_check']=precision_transport_check(p,fvhi)
        out['cases'].append(row)
        Path('analysis/generic_che_forward_numeric.json').write_text(json.dumps(record(out),indent=2)+'\n')
    return out

if __name__=='__main__':main()
