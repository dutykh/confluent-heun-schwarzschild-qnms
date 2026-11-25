"""Validated scalar-Schwarzschild horizon initialization and finite transport.

State: h=exp(i*Omega*x)u_H, v=h_x, x=z+log(z-1), Omega=2M*omega.
All arithmetic enclosures use mpmath.iv; elementary bounds use cert.elementary's
MPFI bridge. No numerical Taylor remainder or roundoff estimate is inferred from
agreement. Fundamental matrices use an O(N) polynomial-ODE Taylor recurrence,
explicit Cauchy remainders, and complex-disk radii to control wrapping.

transport_horizon(omega, ell, anchor, R, dps=100, ...) encloses h,h_x at
z=anchor+iR. omega may be an iv.mpc rectangle or exact decimal box dict.
The default path is 1.5 -> anchor -> anchor+iR, with anchor>1, R>=0.
The returned h,hx are rectangular iv.mpc enclosures of the complex-disk results.
"""
from __future__ import annotations
from contextlib import contextmanager
import time
from mpmath import mp,iv
from cert.elementary import cexp,rexp,rsqrt,cabs


def lower(x):
    return mp.make_mpf(iv.mpf(x)._mpi_[0])


def upper(x):
    return mp.make_mpf(iv.mpf(x)._mpi_[1])


def point(x):
    return iv.mpf(x)


def up_abs(z):
    return upper(cabs(as_complex(z)))


def down_abs(z):
    return lower(cabs(as_complex(z)))


@contextmanager
def precision(dps):
    oldmp,oldiv=mp.prec,iv.prec
    mp.dps=int(dps);iv.dps=int(dps)
    try: yield
    finally:mp.prec=oldmp;iv.prec=oldiv


def as_complex(z):
    if isinstance(z,dict):
        real=z.get('re',z.get('real'));imag=z.get('im',z.get('imag'))
        return iv.mpc(iv.mpf(real),iv.mpf(imag))
    if hasattr(z,'_mpci_'):
        return z
    if hasattr(z,'real') and hasattr(z,'imag'):
        return iv.mpc(z.real,z.imag)
    return iv.mpc(z,0)


def midpoint(z):
    z=as_complex(z)
    # mp midpoint arithmetic only chooses an exactly representable center.
    # Any center-rounding error is included by ball_from_rectangle below.
    re=(lower(z.real)+upper(z.real))/2
    im=(lower(z.imag)+upper(z.imag))/2
    return mp.mpc(re,im)


def ball_from_rectangle(z,extra=0):
    z=as_complex(z);c=midpoint(z)
    r=upper(cabs(z-as_complex(c))+iv.mpf(extra))
    return c,r


def ball_rectangle(center,radius):
    rr=iv.mpf([-radius,radius])
    return iv.mpc(iv.mpf(center.real)+rr,iv.mpf(center.imag)+rr)


def _iv_max(*xs):
    return iv.mpf(max(upper(x) for x in xs))


def horizon_enclosure(omega,ell,dps=100,terms=220,rho='0.9'):
    """Certified h,h_x at z=3/2, uniformly for omega in its input rectangle.

    The infinite Frobenius tail is bounded by a proved geometric coefficient
    majorant. A resonance box containing any recurrence denominator zero is
    rejected. Returns rectangular enclosures plus disk centers/radii.
    """
    with precision(dps):
        om=as_complex(omega);L=int(ell)*(int(ell)+1)
        if ell<0 or int(ell)!=ell:raise ValueError('ell must be a nonnegative integer')
        N=int(terms)
        if N<3:raise ValueError('Need at least three horizon terms')
        rho=iv.mpf(rho);w=-iv.mpf(1)/2;z0=iv.mpf(3)/2
        if not(lower(rho)>upper(-w) and upper(rho)<1):raise ValueError('Need .5<rho<1')
        lam=2*iv.j*om
        coeff=[iv.mpc(1)];prev=iv.mpc(0);cur=iv.mpc(1)
        H=cur;Hz=iv.mpc(0)
        for n in range(1,N+1):
            den=n*(n-lam)
            if down_abs(den)<=0:raise ValueError(f'Horizon denominator may vanish at n={n}')
            nxt=((n*(n-1)+4*om**2-L)*cur-(lam*(n-1)+4*om**2)*prev)/den
            coeff.append(nxt);H+=nxt*w**n;Hz-=n*nxt*w**(n-1)
            prev,cur=cur,nxt
        # For all n>=N+1, |A_n|<=abar, |B_n|<=bbar.
        O=iv.mpf(up_abs(om));lamabs=2*O;C=4*O**2+L;n0=iv.mpf(N+1)
        den=1-lamabs/n0
        if lower(den)<=0:raise ValueError('Too few horizon terms for the tail proof')
        abar=(1+C/n0**2)/den
        bbar=(lamabs/n0+C/n0**2)/den
        contraction=abar*rho+bbar*rho**2
        if upper(contraction)>=1:raise ValueError('Geometric coefficient majorant not contractive; increase terms')
        CH=_iv_max(cabs(coeff[N-1])*rho**(N-1),cabs(coeff[N])*rho**N)
        q=(-w)/rho
        tailH=CH*q**(N+1)/(1-q)
        tailHz=CH/(-w)*q**(N+1)*((N+1)/(1-q)+q/(1-q)**2)
        G=z0*cexp(2*iv.j*om*(z0-1));f=(z0-1)/z0
        h=G*H;v=f*G*(Hz+(2*iv.j*om+1/z0)*H)
        rh=upper(cabs(G)*tailH)
        rv=upper(cabs(f*G)*(tailHz+cabs(2*iv.j*om+1/z0)*tailH))
        hc,hr=ball_from_rectangle(h,rh);vc,vr=ball_from_rectangle(v,rv)
        return {'h':ball_rectangle(hc,hr),'hx':ball_rectangle(vc,vr),
                'centers':(hc,vc),'radii':(hr,vr),'z':iv.mpc(z0),
                'horizon_terms':N,'horizon_tail_H':upper(tailH),
                'horizon_tail_Hz':upper(tailHz),'majorant_contraction':upper(contraction),
                'horizon_coefficient_majorant':upper(CH),'dps':dps}


def _fundamental_polynomial(z,omega,ell,step,order):
    """Interval fundamental polynomial for (h,h_x), using O(order) recurrence.

    z^2(z-1)h''+(z-2iOmega z^3)h'-(Lz+1)h=0.
    """
    z=as_complex(z);om=as_complex(omega);hh=as_complex(step);L=ell*(ell+1)
    aa=[z**3-z**2,3*z**2-2*z,3*z-1,iv.mpc(1)]
    bb=[z-2*iv.j*om*z**3,1-6*iv.j*om*z**2,-6*iv.j*om*z,-2*iv.j*om]
    cc=[-L*z-1,iv.mpc(-L)]
    if down_abs(aa[0])<=0:raise ValueError('Taylor center box meets a singularity')
    columns=[]
    for col in range(2):
        y=[iv.mpc(1 if col==0 else 0),iv.mpc(0) if col==0 else z/(z-1)]
        for n in range(order-1):
            rhs=(aa[1]*n*(n+1)+bb[0]*(n+1))*y[n+1]
            rhs+=(aa[2]*n*(n-1)+bb[1]*n+cc[0])*y[n]
            if n>=1:rhs+=(aa[3]*(n-1)*(n-2)+bb[2]*(n-1)+cc[1])*y[n-1]
            if n>=2:rhs+=bb[3]*(n-2)*y[n-2]
            y.append(-rhs/(aa[0]*(n+2)*(n+1)))
        val=y[-1];der=order*y[-1]
        for k in range(order-1,-1,-1):val=val*hh+y[k]
        for k in range(order-1,0,-1):der=der*hh+k*y[k]
        v=(z+hh-1)/(z+hh)*der
        columns.append((val,v))
    return ((columns[0][0],columns[1][0]),(columns[0][1],columns[1][1]))


def _cauchy_remainder(z,omega,ell,step,order,radius_fraction='0.8'):
    z=as_complex(z);om=as_complex(omega);step=as_complex(step)
    d0,d1=down_abs(z),down_abs(z-1)
    r=iv.mpf(lower(iv.mpf(radius_fraction)*min(d0,d1)))
    if lower(r)<=0:raise ValueError('No positive analytic disk')
    bmax=1+1/(iv.mpf(d1)-r)
    vpot=ell*(ell+1)/(iv.mpf(d0)-r)**2+1/(iv.mpf(d0)-r)**3
    M=_iv_max(bmax,vpot+2*cabs(om)*bmax)
    S=rexp(M*r);q=cabs(step)/r
    if upper(q)>=1:raise ValueError('Taylor step leaves the certified analytic disk')
    remh=S*q**(order+1)/(1-q)
    remhz=S/r*q**order*((order+1)/(1-q)+q/(1-q)**2)
    fnext=cabs((z+step-1)/(z+step))
    return (upper(remh),upper(fnext*remhz)),{'analytic_radius':lower(r),'q':upper(q),'Mr':upper(M*r)}


def validated_step(z,omega,ell,step,centers,radii,order=96):
    T=_fundamental_polynomial(z,omega,ell,step,order)
    remainder,meta=_cauchy_remainder(z,omega,ell,step,order)
    totalnorm=_iv_max(cabs(as_complex(centers[0]))+radii[0],cabs(as_complex(centers[1]))+radii[1])
    newc=[];newr=[]
    for i in range(2):
        box=T[i][0]*as_complex(centers[0])+T[i][1]*as_complex(centers[1])
        extra=cabs(T[i][0])*radii[0]+cabs(T[i][1])*radii[1]+iv.mpf(remainder[i])*totalnorm
        c,r=ball_from_rectangle(box,upper(extra));newc.append(c);newr.append(r)
    return tuple(newc),tuple(newr),{'remainder':remainder,**meta}


def _segment(zstart,zend,omega,ell,centers,radii,order,step_fraction,records):
    zstart,zend=as_complex(zstart),as_complex(zend)
    delta=zend-zstart;length=up_abs(delta)
    if length==0:return centers,radii
    omabs=up_abs(omega);t=mp.mpf(0);z=zstart
    while t<1:
        zm=midpoint(z)
        desired=min(mp.mpf(step_fraction)*min(abs(zm),abs(zm-1)),1/omabs if omabs else mp.inf)
        dt=min(mp.mpf(1)-t,desired/length)
        if dt<=0:raise ArithmeticError('Nonpositive transport step')
        tn=min(mp.mpf(1),t+dt)
        if tn==t:raise ArithmeticError('Transport parameter stagnated')
        zn=zstart+iv.mpf(tn)*delta;step=zn-z
        centers,radii,meta=validated_step(z,omega,ell,step,centers,radii,order)
        records.append({'z_mid':midpoint(zn),'radii':radii,**meta})
        t,z=tn,zn
    return centers,radii


def transport_horizon(omega,ell,anchor='1.5',R=80,dps=100,order=96,
                      step_fraction='0.18',horizon_terms=220,return_steps=False):
    """Rigorous h,h_x boxes at anchor+iR; accepts an omega rectangle.

    Default settings target absolute radii far below1e-43 in tested QNM examples.
    Actual radii are returned; no tolerance is silently assumed achieved.
    """
    begin=time.perf_counter()
    with precision(dps):
        om=as_complex(omega);a=iv.mpf(anchor);rr=iv.mpf(R)
        if lower(a)<=1 or lower(rr)<0:raise ValueError('Need anchor>1 and R>=0')
        ini=horizon_enclosure(om,ell,dps,horizon_terms)
        centers,radii=ini['centers'],ini['radii'];records=[]
        z0=iv.mpc(iv.mpf(3)/2)
        centers,radii=_segment(z0,iv.mpc(a),om,ell,centers,radii,order,step_fraction,records)
        centers,radii=_segment(iv.mpc(a),iv.mpc(a,rr),om,ell,centers,radii,order,step_fraction,records)
        result={'h':ball_rectangle(centers[0],radii[0]),'hx':ball_rectangle(centers[1],radii[1]),
                'centers':centers,'radii':radii,'target':iv.mpc(a,rr),
                'omega':om,'ell':ell,'dps':dps,'order':order,'step_fraction':step_fraction,
                'steps':len(records),'horizon_terms':horizon_terms,
                'initial_radii':ini['radii'],'horizon_tail_H':ini['horizon_tail_H'],
                'horizon_tail_Hz':ini['horizon_tail_Hz'],
                'max_local_remainder':max((max(x['remainder']) for x in records),default=mp.mpf(0)),
                'elapsed_seconds':time.perf_counter()-begin}
        if return_steps:result['step_records']=records
        return result


def exact_box(z):
    """Lossless dyadic endpoints, for transport between Python certificate layers.

    Decimal display conversions are deliberately not used here. Each endpoint
    is an exact rational integer / power-of-two denominator string.
    """
    def rational(raw):
        sign,man,exp,bc=raw
        num=(-1 if sign else 1)*man
        if exp>=0:return str(num<<exp)
        return f'{num}/{1<<(-exp)}'
    z=as_complex(z)
    return {'re':[rational(x) for x in z.real._mpi_],
            'im':[rational(x) for x in z.imag._mpi_]}
