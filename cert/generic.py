#!/usr/bin/env python3
"""Rigorous finite enclosures of two generic CHE connection coefficients.

Arithmetic: mpmath.iv rational operations plus MPFI transcendental enclosures.
This computes a single forward adjoint coefficient; it does not integrate an
outgoing solution or infer numerical error from agreement at two precisions.
See generic_notes.md for the explicit analytic majorants and branch paths.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath import iv
from elementary import cabs, cexp, clog, rexp, rlog, endpoints, backend_info


CASES = {
    'complex_A': [('-1.2','.3'),('.3','.2'),('-.4','.1'),('.2','-.15'),('.6','.05')],
    'complex_B': [('.8','.6'),('-.2','.15'),('.4','-.25'),('.35','.2'),('-.1','.4')],
}


def rat(x):
    x=F(x)
    return iv.mpf(x.numerator)/x.denominator


def point(z):
    return iv.mpc(rat(z[0]),rat(z[1]))


def upper_abs(z):
    return cabs(z).b


def box(e):
    e=e.b
    side=iv.mpf([-e,e])
    return iv.mpc(side,side)


def dictionary(raw):
    alpha,beta,gamma,delta,eta=[point(x) for x in raw]
    A=delta+alpha*(beta+gamma+2)/2
    B=eta+beta/2+(gamma-alpha)*(beta+1)/2
    rho=-A/alpha
    u=beta+1+2*rho
    v=gamma+1
    C=rho*(rho+beta)
    D=A+B+rho*(gamma+1)
    return alpha,beta,gamma,A,B,rho,u,v,C,D


def seed(data,center=(F(1,64),F(1,64)),degree=220):
    alpha,beta,gamma,A,B,rho,u,v,C,D=data
    z=point(center)
    c0=(1+beta.real).a
    if not (c0>0):
        raise ArithmeticError("The seed denominator lower bound is nonpositive")
    K=upper_abs(beta+gamma+1-alpha)
    S=max((1+K+upper_abs(B)/c0).b,
          (upper_abs(alpha)+upper_abs(A-alpha)/c0).b)
    M=(S+1).b
    rr=(M*upper_abs(z)).b
    if not (rr<rat('0.5')):
        raise ArithmeticError("The seed series ratio is too large")
    aa=[iv.mpc(1)]
    previous=iv.mpc(0)
    for n in range(degree):
        new=((n*(n+beta+gamma+1-alpha)+B)*aa[-1]
             +(alpha*(n-1)+A)*previous)/((n+1)*(n+beta+1))
        previous=aa[-1]
        aa.append(new)
    h=aa[-1]
    hp=iv.mpc(0)
    for coefficient in reversed(aa[:-1]):
        hp=hp*z+h
        h=h*z+coefficient
    err=(rr**(degree+1)/(1-rr)).b
    derr=(M*rr**degree*((degree+1)-degree*rr)/(1-rr)**2).b
    h+=box(err)
    hp+=box(derr)
    gauge=cexp(alpha*z+(u-rho)*clog(z)+v*clog(z-1))
    q=gauge*h
    qp=gauge*(hp+(alpha+(u-rho)/z+v/(z-1))*h)
    return q,qp,dict(degree=degree,center=center,coefficient_majorant=M,
                     ratio=rr,h_tail=err,hprime_tail=derr)


def shifted(poly,z):
    return [sum((poly[k]*comb(k,j)*z**(k-j) for k in range(j,len(poly))),iv.mpc(0))
            for j in range(len(poly))]


def ode_polynomials(data):
    alpha,_,_,_,_,_,u,v,C,D=data
    return ([iv.mpc(0),iv.mpc(0),iv.mpc(1),iv.mpc(-2),iv.mpc(1)],
            [iv.mpc(0),-u,-alpha+2*u+v,2*alpha-u-v,-alpha],
            [C+u,-2*(C+u)-D,C+u+D+v])


def disk_radius(center):
    a,b=center
    # These exact rational infinity norms are lower bounds on the Euclidean
    # distances to the singular points. Keeping d<=1 improves Cauchy constants.
    l0=max(abs(a),abs(b))
    l1=max(abs(a-1),abs(b))
    return min(l0,l1,F(2))/2,l0,l1


def taylor_step(data,polys,center,state,step,degree):
    z=point(center)
    hh=point(step)
    q,qp=state
    d,l0,l1=disk_radius(center)
    di=rat(d)
    r0=rat(l0-d)
    r1=rat(l1-d)
    alpha,_,_,_,_,_,u,v,C,D=data
    P=(upper_abs(alpha)+upper_abs(u)/r0+upper_abs(v)/r1).b
    Q=(upper_abs(C+u)/r0**2+upper_abs(D)/(r0*r1)+upper_abs(v)/r1**2).b
    K=max(iv.mpf(1),(di*P+di*di*Q).b)
    state_bound=max(upper_abs(q),(di*upper_abs(qp)).b)
    # Exact rational l1 step norm is an upper bound for |step|.
    rr=rat((abs(step[0])+abs(step[1]))/d)
    if not (rr<=rat('0.25')):
        raise ArithmeticError("The Taylor step ratio exceeds its proved limit")
    rem=(state_bound*rexp(K)*rr**(degree+1)/(1-rr)).b
    Ap,Bp,Cp=[shifted(poly,z) for poly in polys]
    aa=[q,qp]
    for n in range(degree):
        rhs=iv.mpc(0)
        for j in range(1,min(4,n)+1):
            rhs+=Ap[j]*(n-j+2)*(n-j+1)*aa[n-j+2]
        for j in range(min(4,n)+1):
            rhs+=Bp[j]*(n-j+1)*aa[n-j+1]
        for j in range(min(2,n)+1):
            rhs+=Cp[j]*aa[n-j]
        aa.append(-rhs/(Ap[0]*(n+2)*(n+1)))
    value=aa[degree]
    derivative=(degree+1)*aa[degree+1]
    for n in range(degree-1,-1,-1):
        value=value*hh+aa[n]
        derivative=derivative*hh+(n+1)*aa[n+1]
    return (value+box(rem),derivative+box(rem/di)),rem,K


def transport(data,state,center,targets,degree):
    polys=ode_polynomials(data)
    count=0
    max_remainder=iv.mpf(0)
    max_exponent=iv.mpf(0)
    for target in targets:
        target=tuple(map(F,target))
        while center!=target:
            d,_,_=disk_radius(center)
            step=(target[0]-center[0],target[1]-center[1])
            while abs(step[0])+abs(step[1])>d/4:
                step=(step[0]/2,step[1]/2)
            state,err,K=taylor_step(data,polys,center,state,step,degree)
            max_remainder=max(max_remainder,err)
            max_exponent=max(max_exponent,K)
            center=(center[0]+step[0],center[1]+step[1])
            count+=1
            if count%100==0:
                print('  Taylor steps',count,'center',str(center[0])[:20],str(center[1])[:20],flush=True)
    return state,dict(steps=count,degree=degree,max_local_q_tail=max_remainder,
                      max_Cauchy_exponent=max_exponent)


def coefficient(data,state,R,ray,N):
    if R<=1 or ray not in (-1,1) or N<1:
        raise ValueError('Need R>1, ray in {-1,+1}, and N>=1')
    alpha,_,_,_,_,_,u,v,C,D=data
    q,qp=state
    z=iv.mpc(R*ray)
    aa=[iv.mpc(1)]
    for n in range(N):
        am=aa[n-1] if n else 0
        aa.append(((n*(n+1-u-v+alpha)+C+D)*aa[n]
                  -((n-1)*(n-u)+C)*am)/(alpha*(n+1)))
    PN=aa[-1]
    PNds=iv.mpc(0)
    for a in reversed(aa[:-1]):
        PNds=PNds/z+PN
        PN=PN/z+a
    PNp=-PNds/z**2
    pz=alpha+u/z+v/(z-1)
    finite=((pz*PN+PNp)*q-PN*qp)/alpha
    E=(N*(N+1-u-v+alpha)+C+D)*aa[N]-((N-1)*(N-u)+C)*aa[N-1]
    G=-(N*(N+1-u)+C)*aa[N]
    sigma=(-alpha.real*ray-upper_abs(u)/R-upper_abs(v)/(R-1)).a
    if not (sigma>0):
        raise ArithmeticError("The terminal ray fails strict damping")
    residual=(rat(F(R,R-1))*(upper_abs(E)*rat(F(1,(N+1)*R**(N+1)))
                            +upper_abs(G)*rat(F(1,(N+2)*R**(N+2))))).b
    W=(upper_abs(C+u)/R+upper_abs(D)*rlog(rat(F(R,R-1)))+upper_abs(v)/(R-1)).b
    qmajorant=((upper_abs(q)+upper_abs(qp)/sigma)*rexp(W/sigma)).b
    tail=(residual*qmajorant/cabs(alpha).a).b
    return finite+box(tail),dict(finite_part=finite,tail_bound=tail,sigma=sigma,
             residual_integral_bound=residual,W_bound=W,q_tail_majorant=qmajorant,
             q_at_R=q,qprime_at_R=qp,ray=ray,R=R,N=N)


def exact_interval(x):
    a,b=endpoints(x)
    return {'lower':str(a),'upper':str(b)}


def decimal_interval(x,digits=60):
    # Decimal display is explicitly not used as the certification boundary.
    a,b=endpoints(x)
    with mp.workdps(digits):
        return [mp.nstr(mp.mpf(t.numerator)/t.denominator,digits) for t in (a,b)]


def outward_fixed_interval(x,places=40):
    """Decimal endpoint strings rounded outwards by exact integer arithmetic."""
    a,b=endpoints(x)
    scale=10**places
    lo=(a*scale).__floor__()
    hi=(b*scale).__ceil__()
    def fixed(n):
        prefix='-' if n<0 else ''
        n=abs(n)
        return prefix+str(n//scale)+'.'+str(n%scale).zfill(places)
    return [fixed(lo),fixed(hi)]


def serialize(x):
    if isinstance(x,dict):return {k:serialize(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serialize(v) for v in x]
    if isinstance(x,F):return str(x)
    if hasattr(x,'_mpi_'):
        return {**exact_interval(x),'decimal_display':decimal_interval(x)}
    if hasattr(x,'_mpci_'):
        return {'real':serialize(x.real),'imag':serialize(x.imag)}
    return x


def run_case(name,R=64,N=32,degree=192):
    if degree<1:
        raise ValueError('Taylor degree must be at least one')
    started=time.time()
    data=dictionary(CASES[name])
    q,qp,seedinfo=seed(data)
    center=(F(1,64),F(1,64))
    targets=[(F(1,2),F(1,2)),(2,2)]
    if name=='complex_A':
        ray=1
        targets.extend([(4,2),(4,0),(R,0)])
    else:
        ray=-1
        targets.extend([(-2,2),(-4,0),(-R,0)])
    state,transportinfo=transport(data,(q,qp),center,targets,degree)
    value,tailinfo=coefficient(data,state,R,ray,N)
    result=dict(name=name,parameters=CASES[name],path=[center]+targets,
                coefficient=value,seed=seedinfo,transport=transportinfo,
                outward_decimal_enclosure_40_places={
                    'real':outward_fixed_interval(value.real),
                    'imag':outward_fixed_interval(value.imag)},
                coefficient_certificate=tailinfo,elapsed_seconds=time.time()-started)
    print(name,'real',decimal_interval(value.real,30),'imag',decimal_interval(value.imag,30),flush=True)
    print('tail',decimal_interval(tailinfo['tail_bound'],10),'elapsed',time.time()-started,flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=list(CASES))
    parser.add_argument('--R',type=int,default=64)
    parser.add_argument('--N',type=int,default=32)
    parser.add_argument('--degree',type=int,default=192)
    parser.add_argument('--dps',type=int,default=160)
    parser.add_argument('--output',default='cert/results_generic.json')
    args=parser.parse_args()
    iv.dps=args.dps
    result={'arithmetic':backend_info(),'working_decimal_digits':args.dps,
      'certificate':'Exact rational endpoint enclosures; decimal strings are display only.',
      'cases':[]}
    for name in ([args.case] if args.case else list(CASES)):
        result['cases'].append(run_case(name,args.R,args.N,args.degree))
        Path(args.output).write_text(json.dumps(serialize(result),indent=2)+'\n')
    return result


if __name__=='__main__':main()
