#!/usr/bin/env python3
"""Exact-rational three-point/Cauchy Rouché root isolation.

The checker itself uses only Python's standard library. mpmath.iv rectangles
can be supplied via complex_box_from_iv(); their binary endpoints are converted
EXACTLY, never through formatted decimals.

Analytic hypotheses (must be supplied/proved outside the point evaluator):
  F is holomorphic on a neighborhood of the closed disk |z-c| <= outer_radius,
  and |F| <= outer_sup on its boundary. Three point rectangles enclose the exact
  values at c, c+h, c-h, where h is a positive real rational.

A successful strict rational inequality proves exactly one zero, counted with
multiplicity, in |z-c| < root_radius. In particular that zero is simple.

CLI: python3 cert/isolate.py input.json [--output certificate.json]
See make_certificate() and schwarzschild_outer_bound() for programmatic use.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import argparse
import json


def rat(x) -> F:
    if isinstance(x, F):
        return x
    if isinstance(x, bool) or isinstance(x, float):
        raise TypeError('Use exact integers, decimal strings, or rational strings; floats are rejected')
    return F(x)


def ceil_fraction(x: F) -> int:
    return -((-x.numerator)//x.denominator)


def ceil_dyadic(x: F, bits: int = 128) -> F:
    """Round upward exactly, including for negative x."""
    unit = 1 << bits
    return F(ceil_fraction(x*unit), unit)


def floor_dyadic(x: F, bits: int = 128) -> F:
    unit = 1 << bits
    return F((x*unit).numerator//(x*unit).denominator, unit)


def sqrt_bounds(x: F, bits: int = 128) -> tuple[F,F]:
    """Exact dyadic enclosure, via integer square root; no floating arithmetic."""
    x = rat(x)
    if x < 0:
        raise ValueError('square root of negative rational')
    scale = 1 << bits
    n = isqrt((x.numerator*scale*scale)//x.denominator)
    lo = F(n,scale)
    hi = lo if lo*lo == x else F(n+1,scale)
    assert lo*lo <= x <= hi*hi
    return lo,hi


def exp_upper(x: F, bits: int = 128, terms: int = 48) -> F:
    """Rigorous rational upper bound for exp(x), x>=0.

    Reduce to x/2^k <= 1/2; sum Taylor terms 0..terms; bound remaining
    positive tail geometrically; square back with outward dyadic rounding.
    """
    x=rat(x)
    if x < 0:
        raise ValueError('exp_upper only needs nonnegative inputs')
    if x == 0:
        return F(1)
    y=x
    squarings=0
    while y > F(1,2):
        y /= 2
        squarings += 1
    term=F(1)
    total=F(1)
    for n in range(1,terms+1):
        term *= y/n
        total += term
    nextterm=term*y/(terms+1)
    tail=nextterm/(1-y/F(terms+2))
    upper=ceil_dyadic(total+tail,bits)
    for _ in range(squarings):
        upper=ceil_dyadic(upper*upper,bits)
    return upper


def _log_1_to_2_upper(x: F, terms: int = 48) -> F:
    if not 1 <= x <= 2:
        raise ValueError('internal log interval error')
    t=(x-1)/(x+1)
    power=t
    total=F(0)
    for n in range(terms):
        total += power/F(2*n+1)
        power *= t*t
    # First omitted denominator 2*terms+1; all later denominators are larger.
    tail=power/(F(2*terms+1)*(1-t*t))
    return 2*(total+tail)


def log_upper_ge_one(x: F, bits: int = 128, terms: int = 48) -> F:
    """Rigorous log(x) upper bound for rational x>=1, using atanh series."""
    x=rat(x)
    if x < 1:
        raise ValueError('log argument must be >=1')
    k=0
    while x > 2:
        x /= 2
        k += 1
    val=k*_log_1_to_2_upper(F(2),terms)+_log_1_to_2_upper(x,terms)
    return ceil_dyadic(val,bits)


@dataclass(frozen=True)
class RealBox:
    lo:F
    hi:F

    def __post_init__(self):
        object.__setattr__(self,'lo',rat(self.lo))
        object.__setattr__(self,'hi',rat(self.hi))
        if self.lo > self.hi:
            raise ValueError('reversed interval')

    @property
    def mid(self): return (self.lo+self.hi)/2

    @property
    def radius(self): return (self.hi-self.lo)/2

    @property
    def maxabs(self): return max(abs(self.lo),abs(self.hi))

    def __sub__(self,other): return RealBox(self.lo-other.hi,self.hi-other.lo)

    def divide_positive(self,d):
        d=rat(d)
        if d<=0: raise ValueError('positive denominator required')
        return RealBox(self.lo/d,self.hi/d)


@dataclass(frozen=True)
class ComplexBox:
    re:RealBox
    im:RealBox

    @classmethod
    def from_json(cls,obj):
        return cls(RealBox(*obj['re']),RealBox(*obj['im']))

    @property
    def midpoint(self): return self.re.mid,self.im.mid

    def abs_upper(self,bits=128):
        return sqrt_bounds(self.re.maxabs**2+self.im.maxabs**2,bits)[1]

    def radius_upper(self,bits=128):
        return sqrt_bounds(self.re.radius**2+self.im.radius**2,bits)[1]

    def __sub__(self,other): return ComplexBox(self.re-other.re,self.im-other.im)

    def divide_positive(self,d):
        return ComplexBox(self.re.divide_positive(d),self.im.divide_positive(d))

    def as_json(self):
        return {'re':[str(self.re.lo),str(self.re.hi)],
                'im':[str(self.im.lo),str(self.im.hi)]}


def _mpf_tuple_fraction(tup) -> F:
    sign,man,exponent,bitcount=tup
    if bitcount < 0:
        raise ValueError('nonfinite mpmath endpoint')
    ans=F(int(man))
    if exponent >= 0:
        ans *= 1 << int(exponent)
    else:
        ans /= 1 << int(-exponent)
    return -ans if sign else ans


def complex_box_from_iv(z) -> ComplexBox:
    """Convert finite mpmath.iv.mpc box to exact Fraction endpoints."""
    rlo,rhi=z.real._mpi_
    ilo,ihi=z.imag._mpi_
    return ComplexBox(RealBox(_mpf_tuple_fraction(rlo),_mpf_tuple_fraction(rhi)),
                      RealBox(_mpf_tuple_fraction(ilo),_mpf_tuple_fraction(ihi)))


def _fraction_decimal_hint(x: F, digits: int = 18) -> str:
    """Display-only decimal; it is never used to make a certification decision."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec=digits
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def _measure(x:F):
    return {'exact':str(x),'decimal_hint':_fraction_decimal_hint(x)}


def make_certificate(center,root_radius,outer_radius,outer_sup,h,
                     f0:ComplexBox,fplus:ComplexBox,fminus:ComplexBox,
                     *,bits=192,metadata=None):
    """Check a Rouché inequality using exact Fraction arithmetic.

    center = (exact real, exact imaginary); h is along the real direction.
    No numerical function evaluation is performed here.
    """
    cr,ci=map(rat,center)
    r,R,M,step=map(rat,(root_radius,outer_radius,outer_sup,h))
    if not 0 < r < R: raise ValueError('Require 0 < root_radius < outer_radius')
    if not 0 < step < R: raise ValueError('Require 0 < h < outer_radius')
    if M < 0: raise ValueError('outer_sup must be nonnegative')
    slope_box=(fplus-fminus).divide_positive(2*step)
    dr,di=slope_box.midpoint
    slope_abs_lower=sqrt_bounds(dr*dr+di*di,bits)[0]
    slope_sample_error=slope_box.radius_upper(bits)
    cauchy_derivative_error=M*step*step/(R*(R*R-step*step))
    derivative_error=slope_sample_error+cauchy_derivative_error
    center_abs_upper=f0.abs_upper(bits)
    higher_order_error=M*r*r/(R*(R-r))
    lhs=center_abs_upper+derivative_error*r+higher_order_error
    rhs=slope_abs_lower*r
    gap=rhs-lhs
    success=gap>0
    return {
      'method':'Exact rational three-point Cauchy/Rouche certificate',
      'proved':success,
      'conclusion':('Exactly one zero, counted with multiplicity, in the open root disk; hence it is simple.'
                    if success else 'The sufficient inequality did not pass; no zero count is asserted.'),
      'center':{'re':str(cr),'im':str(ci)},
      'root_radius':str(r),'outer_radius':str(R),'outer_sup':str(M),'sample_offset':str(step),
      'assumptions':[
        'F is holomorphic on an open neighborhood of the closed outer disk.',
        'The supplied M bounds |F| on the entire outer circle.',
        'The three supplied rectangles enclose the exact values at c,c+h,c-h.',
        'Any normalization pole or branch cut is excluded from that outer disk.'
      ],
      'sample_enclosures':{'center':f0.as_json(),'plus':fplus.as_json(),'minus':fminus.as_json()},
      'linear_slope':{'re':str(dr),'im':str(di)},
      'inequality':{
        'center_abs_upper':_measure(center_abs_upper),
        'slope_abs_lower':_measure(slope_abs_lower),
        'sample_slope_error':_measure(slope_sample_error),
        'cauchy_slope_error':_measure(cauchy_derivative_error),
        'higher_order_error':_measure(higher_order_error),
        'left_side_upper':_measure(lhs),'right_side_lower':_measure(rhs),
        'strict_gap':_measure(gap)},
      'metadata':metadata or {}
    }


def _distance_from_real_interval(x,lo,hi):
    return lo-x if x<lo else x-hi if x>hi else F(0)


def schwarzschild_outer_bound(center,outer_radius,ell,anchor='1.5',
                              *,horizon_terms=64,bits=128,pi_upper=F(22,7)):
    """Conservative rational M for |A_in| on a complete frequency disk.

    Assumptions: massless scalar Schwarzschild, Omega=2M*omega, unit horizon
    amplitude, positive-real/lower-half disk, anchor>=1.5, strict vertical-ray
    phase inequality throughout its enclosing frequency box. No infinity
    polynomial or numerical transport is used in this bound.

    Proof: uniform Frobenius majorants at z0=3/2; real-segment Gronwall bound;
    then the bounded forward h Volterra limit A_in. The default pi bound is the
    classical strict inequality pi<22/7. All other transcendental upper bounds
    are constructed by rational series with positive explicit tails.
    """
    pr,im=map(rat,center)
    radius=rat(outer_radius)
    a=rat(anchor)
    pi_upper=rat(pi_upper)
    if pi_upper<F(22,7):
        raise ValueError('This self-contained proof requires pi_upper >= 22/7')
    if not isinstance(ell,int) or ell<0: raise ValueError('ell must be nonnegative integer')
    if not radius>0: raise ValueError('outer radius must be positive')
    if a<F(3,2): raise ValueError('anchor must be at least 3/2')
    pmin=pr-radius
    qcenter=-im
    qmin=qcenter-radius
    qmax=qcenter+radius
    if pmin<=0 or qmin<0:
        raise ValueError('This bound requires the disk inside the lower-right quadrant')
    if qmax*qmax>=4*pmin*pmin*a*(a-1):
        raise ValueError('Uniform strict phase condition failed; increase anchor or shrink disk')
    norm_lo,norm_hi=sqrt_bounds(pr*pr+im*im,bits)
    omega_min=norm_lo-radius
    omega_max=norm_hi+radius
    if omega_min<=0: raise ValueError('disk includes zero')
    k=ceil_dyadic(1/omega_min,bits)
    lam_abs=ceil_dyadic(2*omega_max,bits)
    lam_re_lo,lam_re_hi=2*qmin,2*qmax
    lam_im_min=2*pmin
    L=F(ell*(ell+1))
    N=int(horizon_terms)
    if N<1 or F(N+1)<=lam_abs:
        raise ValueError('Increase horizon_terms so N+1 > sup|lambda|')
    bprev=F(0)
    bcur=F(1)
    coef=[bcur]
    for n in range(1,N+1):
        dx=_distance_from_real_interval(F(n),lam_re_lo,lam_re_hi)
        denroot=sqrt_bounds(dx*dx+lam_im_min*lam_im_min,bits)[0]
        if denroot<=0: raise ValueError('coefficient denominator lower bound is zero')
        num=(F(n*(n-1))+lam_abs*lam_abs+L)*bcur+lam_abs*(n-1+lam_abs)*bprev
        bn=ceil_dyadic(num/(n*denroot),bits)
        coef.append(bn)
        bprev,bcur=bcur,bn
    J=F(N+1)
    aa=(1+(lam_abs*lam_abs+L)/(J*J))/(1-lam_abs/J)
    bb=((lam_abs/J)*(1+lam_abs/J))/(1-lam_abs/J)
    # A>=1, so rho=A+B satisfies rho²−A*rho−B=B(A+B−1)>=0.
    rho=ceil_dyadic(aa+bb,bits)
    C=max(coef[N]/rho**N,coef[N-1]/rho**(N-1))
    C=ceil_dyadic(C,bits)
    t=rho/2
    if t>=1: raise ValueError('geometric horizon tail does not converge; increase horizon_terms')
    Htail=C*t**(N+1)/(1-t)
    Hztail=C*rho*t**N*((N+1)-N*t)/(1-t)**2
    H=ceil_dyadic(sum((coef[n]/2**n for n in range(N+1)),F(0))+Htail,bits)
    Hz=ceil_dyadic(sum((n*coef[n]/2**(n-1) for n in range(1,N+1)),F(0))+Hztail,bits)
    z0=F(3,2)
    G=z0*exp_upper(lam_re_hi*(z0-1),bits)
    h0=ceil_dyadic(G*H,bits)
    hx0=ceil_dyadic(F(1,3)*G*(Hz+(lam_abs+1/z0)*H),bits)
    if a==z0:
        hbound,hxbound=h0,hx0
        real_exponent=F(0)
    else:
        dx=(a-z0)+log_upper_ge_one((a-1)/(z0-1),bits)
        intv=L*(1/z0-1/a)+(1/z0**2-1/a**2)/2
        real_exponent=max(F(1),2*omega_max)*dx+intv
        both=ceil_dyadic(max(h0,hx0)*exp_upper(real_exponent,bits),bits)
        hbound=hxbound=both
    Vray=rat(pi_upper)*L/(2*a)+1/(a*a)
    Mraw=ceil_dyadic((hbound+k*hxbound)*exp_upper(k*Vray,bits),bits)
    M=F(ceil_fraction(Mraw))
    return {
       'outer_sup':str(M),
       'method':'Uniform horizon majorant + real-segment Gronwall + forward ray Volterra',
       'center':{'re':str(pr),'im':str(im)},'outer_radius':str(radius),
       'ell':ell,'anchor':str(a),'horizon_terms':N,'rational_rounding_bits':bits,
       'pi_upper':str(rat(pi_upper)),
       'details':{
          'omega_min_lower':_measure(omega_min),'omega_max_upper':_measure(omega_max),
          'phase_margin_squared':_measure(4*pmin*pmin*a*(a-1)-qmax*qmax),
          'kernel_upper':_measure(k),'lambda_abs_upper':_measure(lam_abs),
          'horizon_H_upper':_measure(H),'horizon_Hz_upper':_measure(Hz),
          'horizon_tail_upper':_measure(Htail),'horizon_derivative_tail_upper':_measure(Hztail),
          'h_at_anchor_upper':_measure(hbound),'hx_at_anchor_upper':_measure(hxbound),
          'real_transport_exponent_upper':_measure(real_exponent),
          'ray_potential_integral_upper':_measure(Vray),'unrounded_M_upper':_measure(Mraw)
       },
       'analytic_domain_justification':[
          'The closed disk lies strictly in Re(Omega)>0 and away from Omega=0.',
          'All unit-horizon normalization poles Omega=−i*j/2 lie outside it.',
          'The outgoing continuation and forward coefficient are holomorphic in this half-plane.',
          'Strict phase decay holds on the entire frequency box enclosing the disk.'
       ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    data=json.loads(args.input.read_text())
    if data.get('operation')=='outer_bound':
        result=schwarzschild_outer_bound((data['center']['re'],data['center']['im']),
             data['outer_radius'],int(data['ell']),data.get('anchor','1.5'),
             horizon_terms=data.get('horizon_terms',64),bits=data.get('bits',128))
    else:
        rect=data['sample_enclosures']
        result=make_certificate((data['center']['re'],data['center']['im']),
          data['root_radius'],data['outer_radius'],data['outer_sup'],data['sample_offset'],
          ComplexBox.from_json(rect['center']),ComplexBox.from_json(rect['plus']),
          ComplexBox.from_json(rect['minus']),bits=data.get('bits',192),metadata=data.get('metadata'))
    output=json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output,end='')

if __name__=='__main__':
    main()
