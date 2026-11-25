"""Verify the comparison residual identity and uniform rational inequalities."""
from fractions import Fraction
from sympy import symbols, simplify

x, alpha, eps, gamma, b, A = symbols('x alpha eps gamma b A', positive=True)
phi = (1+x)**(-b)
logPprime = -alpha+eps/x+(1+gamma)/(1+x)
ratio = -x*(1+x)*(phi.diff(x,2)+logPprime*phi.diff(x))/phi+A*x
target = b*eps + (A-alpha*b)*x + b*(gamma-b)*x/(1+x)
assert simplify(ratio-target) == 0

alpha=Fraction(7,10)
gamma=Fraction(2,5)
Amin=Fraction(69,100)
epsmax=Fraction(1,1000)
Amax=Amin+alpha*epsmax/2
assert 0 < gamma < Amin/alpha <= Amax/alpha < 1
assert Amin-alpha*gamma>0
assert Amin/8-epsmax>0
assert 3*Fraction(1,8)/(1-3*Fraction(1,8))==Fraction(3,5)<1
print('Comparison residual identity and interval inequalities passed.')
