"""Validated real-axis matching for the resolved apparent-origin CHE example.

Run from project root with python -m cert.resonant. Each accessory bracket is
verified by opposite strict signs of the real incoming-coefficient functional.
The Sturm identity in the manuscript supplies uniqueness and simplicity.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
from mpmath import iv

from .elementary import endpoints, rexp, rlog, backend_info

A0 = F(69, 100)
ALPHA = F(7, 10)
GAMMA = F(2, 5)
CASES = [(F(0), '0.7474537865533673590405589901'),
         (F(1, 10**7), '0.7474538155896238054048569455'),
         (F(1, 10**5), '0.7474566901786279866189859707'),
         (F(1, 10**3), '0.7477441452401204816551441447')]


def interval(value):
    value = F(value)
    return iv.mpf(value.numerator) / value.denominator


def rational_ends(value):
    return [str(x) for x in endpoints(value)]


def symmetric(bound):
    return iv.mpf([-bound.b, bound.b])


def initial(eps, t, degree=270):
    alpha = interval(ALPHA)
    gamma = interval(GAMMA)
    A = interval(A0) + eps*alpha/2
    B = t*eps
    rho = -A/alpha
    u, v = eps + 2*rho, 1+gamma
    C = rho*(rho-1+eps)
    D = A+B+rho*v
    z = interval(F(-1, 8))
    h = [iv.mpf(1), t]
    for n in range(1, degree):
        h.append(((n*(n+eps+gamma-alpha)+B)*h[-1]+
                  (alpha*(n-1)+A)*h[-2])/((n+1)*(n+eps)))
    value, deriv = h[-1], iv.mpf(0)
    for coefficient in reversed(h[:-1]):
        deriv, value = deriv*z+value, value*z+coefficient
    ratio = interval(F(3, 8))
    value += symmetric(ratio**(degree+1)/(1-ratio))
    deriv += symmetric(3*ratio**degree*((degree+1)-degree*ratio)/(1-ratio)**2)
    gauge = rexp(alpha*z+(u-rho)*rlog(-z)+v*rlog(1-z))
    q0 = gauge*value
    d0 = gauge*(deriv+(alpha+(u-rho)/z+v/(z-1))*value)
    return (alpha,u,v,C,D), (q0,d0), (value,deriv)


def continue_negative_ray(data, state, radius, degree=210):
    alpha,u,v,C,D = data
    polys = ([0,0,1,-2,1],[0,-u,-alpha+2*u+v,2*alpha-u-v,-alpha],
             [C+u,-2*(C+u)-D,C+u+D+v])
    x = F(1, 8)
    steps = 0
    while x < radius:
        d = x/2
        h = min(d/4, radius-x)
        z = interval(-x)
        dz = interval(-h)
        disk = interval(d)
        r0, r1 = interval(x-d), interval(x+1-d)
        P = abs(alpha)+abs(u)/r0+abs(v)/r1
        Q = abs(C+u)/r0**2+abs(D)/(r0*r1)+abs(v)/r1**2
        norm = max(iv.mpf(1),(disk*P+disk**2*Q).b)
        q,dq = state
        tau = interval(h/d)
        tail = max(abs(q).b,(disk*abs(dq)).b)*rexp(norm)*tau**(degree+1)/(1-tau)
        shifted = [[sum((poly[k]*comb(k,j)*z**(k-j)
                        for k in range(j,len(poly))),iv.mpf(0))
                    for j in range(len(poly))] for poly in polys]
        a,b,c = shifted
        aa = [q,dq]
        for n in range(degree):
            total = sum((a[j]*(n-j+2)*(n-j+1)*aa[n-j+2]
                         for j in range(1,min(4,n)+1)), iv.mpf(0))
            total += sum((b[j]*(n-j+1)*aa[n-j+1] for j in range(min(4,n)+1)),iv.mpf(0))
            total += sum((c[j]*aa[n-j] for j in range(min(2,n)+1)),iv.mpf(0))
            aa.append(-total/(a[0]*(n+2)*(n+1)))
        q, dq = aa[degree], (degree+1)*aa[degree+1]
        for n in range(degree-1,-1,-1):
            q = q*dz+aa[n]
            dq = dq*dz+(n+1)*aa[n+1]
        state = (q+symmetric(tail),dq+symmetric(tail/disk))
        x += h
        steps += 1
    return state,steps


def coefficient(eps, t, radius=60, degree=42):
    if abs(eps)>F(1,1000) or abs(t)>1:
        raise ValueError("Seed majorant requires |epsilon|<=1/1000 and |t|<=1")
    data,state,_ = initial(interval(eps),interval(t))
    (q,qd),steps=continue_negative_ray(data,state,radius)
    alpha,u,v,C,D = data
    z=iv.mpf(-radius)
    p=alpha+u/z+v/(z-1)
    aa=[iv.mpf(1)]
    for n in range(degree):
        prev=aa[n-1] if n else 0
        aa.append(((n*(n+1-u-v+alpha)+C+D)*aa[n]-
                  ((n-1)*(n-u)+C)*prev)/(alpha*(n+1)))
    P,Pd=aa[-1],iv.mpf(0)
    for a in reversed(aa[:-1]):
        Pd,P=Pd/z+P,P/z+a
    Pd=-Pd/z**2
    finite=((p*P+Pd)*q-P*qd)/alpha
    n=degree
    E=(n*(n+1-u-v+alpha)+C+D)*aa[n]-((n-1)*(n-u)+C)*aa[n-1]
    G=-(n*(n+1-u)+C)*aa[n]
    sigma=(alpha-abs(u)/radius-abs(v)/(radius-1)).a
    if sigma<=0: raise ArithmeticError("Inadequate damping margin")
    defect=interval(F(radius,radius-1))*(abs(E)*interval(F(1,(n+1)*radius**(n+1)))+
                                          abs(G)*interval(F(1,(n+2)*radius**(n+2))))
    W=abs(C+u)/radius+abs(D)*rlog(interval(F(radius,radius-1)))+abs(v)/(radius-1)
    remainder=(abs(q)+abs(qd)/sigma)*rexp(W/sigma)*defect/alpha
    out=finite+symmetric(remainder)
    return out,dict(steps=steps, R=radius, N=degree, finite_part=rational_ends(finite),
                    tail_bound=rational_ends(remainder)[1],complete_enclosure=rational_ends(out))


def main():
    iv.dps=170
    data={'arithmetic':backend_info(),'ray':'negative real, constant gauge phase removed',
          'seed_degree':270,'Taylor_degree':210,'brackets':[]}
    # Uniform origin solution at a genuine parameter box containing resonance.
    eps=interval(F(1,1000))*iv.mpf([-1,1])
    t=iv.mpf([interval(F(74,100)).a,interval(F(75,100)).b])
    _,_,(h,hp)=initial(eps,t)
    data['uniform_seed']={'epsilon':['-1/1000','1/1000'],'t':['74/100','75/100'],
                        'z':'-1/8','H':rational_ends(h),'H_prime':rational_ends(hp)}
    for eps,mid in CASES:
        delta=F(1,10**28)
        lo,hi=F(mid)-delta,F(mid)+delta
        left,ld=coefficient(eps,lo)
        right,rd=coefficient(eps,hi)
        if not ((left.a>0 and right.b<0) or (left.b<0 and right.a>0)):
            raise ArithmeticError(f'Endpoint signs unresolved at epsilon={eps}: {left} vs {right}')
        data['brackets'].append(dict(epsilon=str(eps),midpoint_decimal=mid,radius=str(delta),t_lower=str(lo),
                                     t_upper=str(hi),left=ld,right=rd))
        print('epsilon', eps,'center', mid,'radius',delta, 'left',left, 'right',right,flush=True)
    Path(__file__).with_name('results_resonant.json').write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    main()
