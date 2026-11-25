"""Enclose the scalar incoming amplitude on the continued positive-real-part sheet.

All operations in this module are interval rational arithmetic.  Complex
moduli in the tail majorants are bounded by their Cartesian one-norms.
"""
from fractions import Fraction
from mpmath import iv


def rational(value):
    value = Fraction(value)
    return iv.mpf(value.numerator) / value.denominator


def norm_upper(z):
    return (abs(z.real) + abs(z.imag)).b


def error_box(bound):
    extent = bound.b
    side = iv.mpf([-extent, extent])
    return iv.mpc(side, side)


def polynomial(omega, ell, degree, z):
    """P, dP/dz and the two coefficients of the exact differential residual."""
    angular = ell * (ell + 1)
    nu = iv.mpc(0, 2) * omega
    terms = [iv.mpc(1)]
    previous = iv.mpc(0)
    for k in range(degree):
        new = ((k * (k + 1) - angular) * terms[-1] - k*k*previous) / (nu * (k+1))
        previous = terms[-1]
        terms.append(new)
    tail_c = (degree*(degree+1) - angular)*terms[-1] - degree*degree*previous
    tail_d = -(degree+1)**2*terms[-1]
    s = 1/z
    value = terms[-1]
    ds = iv.mpc(0)
    for a_k in reversed(terms[:-1]):
        ds = value + s*ds
        value = a_k + s*value
    return value, -ds*s*s, tail_c, tail_d


def from_endpoint(omega, ell, anchor, radius, degree, h, hx, corrected=True,
                  quadrature_options=None):
    """h,hx are interval enclosures of the exact horizon state at a+iR."""
    a, R = rational(anchor), rational(radius)
    z = iv.mpc(a, R)
    P, Pprime, c, d = polynomial(omega, ell, degree, z)
    nu = iv.mpc(0, 2)*omega
    if not (a > 1 and omega.real.a > 0 and omega.imag.b < 0 and
            omega.real.a*R > -omega.imag.a):
        raise ValueError('Strict tail damping test failed')
    omin = max(omega.real.a, -omega.imag.b)
    k = 1/omin
    # Along ray |z| >= s. Bound exp(x)-1 by x/(1-x).
    potential = ell*(ell+1)/R + 1/(2*R*R)
    x = k*potential
    if not x < 1:
        raise ValueError('Increase terminal radius so the tail majorant is less than one')
    residual = norm_upper(c)/(R**(degree+1)*(degree+1)) + norm_upper(d)/(R**(degree+2)*(degree+2))
    I0, Ie, quad = iv.mpc(0), iv.mpc(0), {}
    if corrected:
        from .quadrature import defect_integrals
        I0, Ie, quad = defect_integrals(omega, c, d, degree, anchor, radius,
                                        return_diagnostics=True, **(quadrature_options or {}))
    # Coefficients of h,hx are collected before interval evaluation.
    state_part = ((nu*P + ((z-1)/z)*Pprime + I0)*h + (-P+(Ie-I0)/nu)*hx)/nu
    free_b = hx/nu
    free_a = h-free_b
    factor = x/(1-x) if corrected else 1/(1-x)
    tail = ((norm_upper(free_a)+norm_upper(free_b))*residual*factor/(2*omin)).b
    return state_part + error_box(tail), dict(tail=tail, finite_part=state_part,
                                                quad=quad, I0=I0, Ie=Ie)


def amplitude(omega, ell, anchor, R, N, corrected=True, transport_options=None,
              quadrature_options=None):
    from .transport import transport_horizon
    info = transport_horizon(omega,ell,anchor,R,**(transport_options or {}))
    h, hx = info['h'], info['hx']
    result, info_a = from_endpoint(omega,ell,anchor,R,N,h,hx,corrected,quadrature_options)
    return result, dict(transport=info, h=h, hx=hx, **info_a)
