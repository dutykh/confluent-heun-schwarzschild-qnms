"""Uniform scalar horizon enclosures when the unit Frobenius solution has a pole.

The physical parameter is lambda=2i Omega. Multiplication by (N-lambda)
is performed symbolically at the obstructed recurrence index N. No division
by an interval containing zero occurs. Run: python -m cert.resonant_horizon.
"""
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
from mpmath import iv
from .elementary import backend_info, cabs, cexp, endpoints


def coefficient_numerator(n, lam, ell_factor, current, previous):
    return (n*(n-1)-lam**2-ell_factor)*current-lam*(n-1-lam)*previous


def divide_at_ordinary_index(n, lam, numerator):
    divisor = n*(n-lam)
    if cabs(divisor).a <= 0:
        raise ValueError(f'Regularized chart meets another resonance at n={n}')
    return numerator/divisor


def primitive_state(lam, index, ell, degree=220):
    """Return h, h_x for h=(index-lambda) exp(i Omega x)u_H at z=3/2."""
    if not (1<=index<degree and ell>=0 and isinstance(ell,int)):
        raise ValueError('Require integer 1<=index<degree and nonnegative integer ell')
    L=ell*(ell+1)
    unscaled=[iv.mpc(1)]
    for n in range(1,index):
        old=unscaled[-2] if n>=2 else 0
        unscaled.append(divide_at_ordinary_index(n,lam,
            coefficient_numerator(n,lam,L,unscaled[-1],old)))
    coeff=[(index-lam)*a for a in unscaled]
    old=unscaled[-2] if index>=2 else 0
    coeff.append(coefficient_numerator(index,lam,L,unscaled[-1],old)/index)
    for n in range(index+1,degree+1):
        coeff.append(divide_at_ordinary_index(n,lam,
             coefficient_numerator(n,lam,L,coeff[-1],coeff[-2])))
    w=-iv.mpf(1)/2
    H=coeff[-1]
    Hw=iv.mpc(0)
    for a in reversed(coeff[:-1]):
        Hw,H=Hw*w+H,H*w+a
    rho=iv.mpf(9)/10
    u=cabs(lam).b
    j=degree+1
    lower=1-u/j
    if lower.a<=0: raise ValueError('Increase horizon degree')
    A=(1+(u*u+L)/j**2)/lower
    B=(u/j+(u*u+L)/j**2)/lower
    contraction=A*rho+B*rho**2
    if contraction.b>=1: raise ValueError('Increase horizon degree')
    C=max((cabs(coeff[-2])*rho**(degree-1)).b,(cabs(coeff[-1])*rho**degree).b)
    q=-w/rho
    tail=C*q**(degree+1)/(1-q)
    dtail=C/(-w)*q**(degree+1)*((degree+1)/(1-q)+q/(1-q)**2)
    def square(e):
        side=iv.mpf([-e.b,e.b])
        return iv.mpc(side,side)
    H+=square(tail)
    Hw+=square(dtail)
    z=iv.mpf(3)/2
    g=z*cexp(lam*(z-1))
    h=g*H
    hx=((z-1)/z)*g*(-Hw+(lam+1/z)*H)
    return h,hx,tail,dtail,coeff[index],contraction


def determinant(index,ell):
    """Integer continuant for the obstruction matrix at lambda=N."""
    L=ell*(ell+1)
    previous,current=0,1
    for k in range(index):
        diagonal=index**2+L-k*(k+1)
        offdiag_product=k*(index-k)*index*(index-k)
        previous,current=current,diagonal*current-offdiag_product*previous
    return current


def serialize(z):
    if hasattr(z,'_mpi_'): return [str(a) for a in endpoints(z)]
    return {'real':serialize(z.real),'imag':serialize(z.imag)}


def main():
    iv.dps=110
    radius=iv.mpf(1)/1024
    offset=iv.mpf([-radius,radius])
    cases=[]
    for ell in range(3):
        for N in range(1,4):
            lam=iv.mpc(iv.mpf(N)+offset,offset)
            h,hx,e,de,p,c=primitive_state(lam,N,ell)
            # This component of the regularized state is separated from zero
            # throughout the full parameter neighborhood.
            if cabs(h).a <= 0: raise ArithmeticError('Frame nonvanishing unresolved')
            obstruction=F(-determinant(N,ell),N*factorial(N-1)**2)
            if not (endpoints(p.real)[0]<=obstruction<=endpoints(p.real)[1]):
                raise ArithmeticError('Rational obstruction not contained')
            cases.append({'ell':ell,'N':N,'h':serialize(h),'hx':serialize(hx),
                          'h_abs_lower':serialize(cabs(h))[0],'H_tail':serialize(e)[1],
                          'Hw_tail':serialize(de)[1], 'contraction':serialize(c)[1],
                          'at_central_resonance_hN_exact':str(obstruction)})
    result={'arithmetic':backend_info(),'working_decimal_digits':110,'horizon_degree':220,
            'chart_half_side_in_lambda':'1/1024','z':'3/2','cases':cases}
    Path(__file__).with_name('results_resonant_horizon.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Uniform lower moduli:',[float(F(a['h_abs_lower'])) for a in cases])


if __name__=='__main__': main()
