"""Outward interval enclosures for the Schwarzschild free-tail integrals.

For z_R=a+iR, m=N+2, this module encloses

    I0 = integral_[z_R,z_R+i infinity] (c z^-m+d z^(-m-1)) dz,
    Ie = integral_[z_R,z_R+i infinity]
           exp(2i omega (z-z_R+Log((z-1)/(z_R-1))))
           (c z^-m+d z^(-m-1)) dz.

The ray has a>1,R>0 and 0<p=inf Re(omega), with p R>q,
q=max(0,-inf Im(omega)). These sufficient restrictions are checked.
The logarithm above is principal: its argument is in the right half-plane.

I0 is evaluated by its exact rational antiderivative. Ie uses symmetric
analytic Taylor quadrature and an explicit infinite-tail bound. ALL error
bounds include their own outward arithmetic. Transcendentals are delegated to
the MPFI bridge in elementary.py, not to mpmath.iv's experimental functions.
The working precision is the caller's current iv.dps.

Proof of the local enclosure. Write g_j(t)=E(t)(z_R+it)^(-j). Its analytic
Taylor coefficients about t=t0 satisfy a rational recurrence using three
previous coefficients,
since P(w) g_j'=Q_j(w) g_j with P=z(z-1), Q_j=-2 omega z^2-i j(z-1).
On |w|=r, with d0<=|z(t0)| and d1<=|z(t0)-1|, r<min(d0,d1),

 |E(t0+w)| <= exp(-mu*t0)
                *exp(2 ||omega||_1 [r-log(1-r/d1)]),
 mu=2p-2q/R>0.

This follows by writing the local logarithm as log(1+iw/(z(t0)-1))
and bounding its absolutely convergent series. Hence a Cauchy majorant M
for c*g_m+d*g_(m+1) is obtained by multiplication by
|c|/(d0-r)^m+|d|/(d0-r)^(m+1). Integrating the Taylor series over [-h,h]
annihilates odd terms. If k is the first even degree omitted, its integrated
remainder is <=2 M h (h/r)^k/[(k+1)(1-(h/r)^2)].

Finally |E(t)|<=exp(-mu*t) on the full real ray, because the derivative of
its log modulus is at most -2p+2q/R. For t>=T, |z_R+it|>=R+T, giving the
tail bound exp(-mu*T)(|c|/(R+T)^m+|d|/(R+T)^(m+1))/mu.
Disk error bounds are enclosed by Cartesian rectangles. They are added to
the interval Taylor sum, so rounding, quadrature truncation and tail errors
are all included. No assertion about an upstream ODE transport is made here.
"""
from __future__ import annotations

from fractions import Fraction
from mpmath import iv

try:
    from .elementary import cexp, clog, rexp, rlog
except ImportError:
    from elementary import cexp, clog, rexp, rlog


def _fraction(value):
    """Interpret an int/string/Fraction as an exact rational input."""
    if isinstance(value, float):
        raise TypeError('Use a decimal string or Fraction for exact geometry')
    return Fraction(value)


def _ivq(value):
    value = _fraction(value)
    return iv.mpf(value.numerator) / value.denominator


def _endpoint_fraction(raw):
    sign, mantissa, exponent, bitcount = raw
    if bitcount < 0 or not isinstance(exponent, int):
        raise ValueError('A finite interval endpoint is required')
    value = Fraction((-1 if sign else 1) * mantissa)
    return value * (2**exponent) if exponent >= 0 else value / (2**(-exponent))


def _lower_fraction(x):
    return _endpoint_fraction(x._mpi_[0])


def _upper_fraction(x):
    return _endpoint_fraction(x._mpi_[1])


def _norm_upper(z):
    """Cartesian one-norm, requiring no transcendental modulus."""
    return (abs(z.real) + abs(z.imag)).b


def _power(z, n):
    """Integer powers by rational interval operations only."""
    if n < 0:
        return 1 / _power(z, -n)
    out = iv.mpf(1) if hasattr(z, '_mpi_') else iv.mpc(1)
    while n:
        if n & 1:
            out *= z
        z *= z
        n >>= 1
    return out


def _error_box(bound):
    side = iv.mpf([-1, 1]) * bound.b
    return iv.mpc(side, side)


def _dyadic_below(value):
    """Largest power of two <= the positive exact Fraction value."""
    if value <= 0:
        raise ValueError('Positive step bound required')
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    candidate = Fraction(2**exponent) if exponent >= 0 else Fraction(1, 2**(-exponent))
    if candidate > value:
        candidate /= 2
    return candidate


def _series(omega, z, exponent, initial, order):
    """Taylor coefficients of E(t)*z(t)^(-exponent), exact interval recurrence."""
    p0 = z*(z-1)
    p1 = iv.j*(2*z-1)
    q0 = -2*omega*z*z-iv.j*exponent*(z-1)
    q1 = -4*iv.j*omega*z+exponent
    q2 = 2*omega
    coefficients = [initial]
    for k in range(order):
        numerator = (q0-k*p1)*coefficients[k]
        if k:
            numerator += (q1+k-1)*coefficients[k-1]
        if k > 1:
            numerator += q2*coefficients[k-2]
        coefficients.append(numerator / ((k+1)*p0))
    return coefficients


def _panel(omega, c, d, m, anchor, radius, left, half_width,
           disk_radius, mu, omega_norm, c_norm, d_norm, order, budget):
    center = left+half_width
    z = iv.mpc(_ivq(anchor), _ivq(radius+center))
    zR = iv.mpc(_ivq(anchor), _ivq(radius))
    h, r, t = map(_ivq, (half_width, disk_radius, center))
    # Lower bounds on distances to the singularities, using exact geometry.
    d0 = _ivq(max(anchor, radius+center))
    d1 = _ivq(max(anchor-1, radius+center))
    if not (r < d0 and r < d1 and h < r):
        raise ArithmeticError('Invalid analytic quadrature disk')
    ratio = h/r
    majorant = rexp(-mu*t + 2*omega_norm*(r-rlog(1-r/d1)))
    majorant *= c_norm/_power(d0-r,m) + d_norm/_power(d0-r,m+1)
    majorant = majorant.b
    while True:
        omitted = order+1 if (order+1) % 2 == 0 else order+2
        remainder = (2*majorant*h*_power(ratio,omitted)
                     / ((omitted+1)*(1-ratio*ratio))).b
        if remainder <= budget:
            break
        order += 8
        if order > 1024:
            raise ArithmeticError('Taylor order exceeds 1024; increase precision or relax tolerance')
    phase = cexp(2*iv.j*omega*(iv.j*t+clog((z-1)/(zR-1))))
    gm0 = phase * _power(z,-m)
    gm = _series(omega,z,m,gm0,order)
    gm1 = _series(omega,z,m+1,gm0/z,order)
    total = iv.mpc(0)
    hp = h
    hh = h*h
    for k in range(0,order+1,2):
        total += (c*gm[k]+d*gm1[k]) * (2*hp/(k+1))
        hp *= hh
    return iv.j*total, remainder, order


def defect_integrals(omega, c, d, N, anchor, radius, *, order=32,
                     rtol='1e-30', atol='0', max_panels=100000,
                     return_diagnostics=False):
    """Return (I0,Iexp), or (I0,Iexp,diagnostics), as complex IV rectangles.

    omega,c,d may be point or rectangular iv.mpc values. Geometry must be
    exact ints, decimal strings or Fractions. N is an integer >=0.
    rtol is relative to the explicit *unweighted defect majorant*
      S=||c||_1/((N+1)R^(N+1))+||d||_1/((N+2)R^(N+2)),
    not to the possibly cancelling integral. The requested absolute analytic
    quadrature error is atol+rtol*S; half is allocated to each of local Taylor
    truncation and the infinite tail. Interval dependency and rounding widths
    are included in the returned enclosure but are not controlled by rtol.
    No global precision is changed, and no point-evaluation oracle is used.
    """
    omega, c, d = tuple(x if hasattr(x,'_mpci_') else iv.mpc(x) for x in (omega,c,d))
    anchor, radius = _fraction(anchor), _fraction(radius)
    if isinstance(N,bool) or int(N) != N or N < 0:
        raise ValueError('N must be an integer >=0')
    N = int(N)
    if int(order) != order or order < 2:
        raise ValueError('order must be an integer >=2')
    order = int(order)
    if anchor <= 1 or radius <= 0:
        raise ValueError('Require exact anchor>1 and radius>0')
    p = _lower_fraction(omega.real)
    q = max(Fraction(0),-_lower_fraction(omega.imag))
    if p <= 0 or p*radius <= q:
        raise ValueError('Require inf Re(omega)>0 and inf Re(omega)*R>max(0,-inf Im(omega))')
    omega_norm = _norm_upper(omega)
    c_norm, d_norm = _norm_upper(c), _norm_upper(d)
    m = N+2
    zR = iv.mpc(_ivq(anchor),_ivq(radius))
    I0 = c*_power(zR,-N-1)/(N+1)+d*_power(zR,-N-2)/(N+2)
    scale = c_norm/((N+1)*_power(_ivq(radius),N+1))
    scale += d_norm/((N+2)*_power(_ivq(radius),N+2))
    target = _ivq(atol)+_ivq(rtol)*scale
    if _fraction(atol) < 0 or _fraction(rtol) < 0:
        raise ValueError('Tolerances must be nonnegative')
    if c_norm == 0 and d_norm == 0:
        diagnostics = dict(panels=0,max_order=0,local_remainder=iv.mpf(0),
                           tail_remainder=iv.mpf(0),analytic_error=iv.mpf(0),
                           truncation_endpoint=Fraction(0),reference_scale=scale,
                           target=target,precision=iv.dps)
        return (I0,iv.mpc(0),diagnostics) if return_diagnostics else (I0,iv.mpc(0))
    if not target.a > 0:
        raise ValueError('A positive quadrature error target is required')
    budget = target.a/2
    mu = _ivq(2*p-2*q/radius)
    T = Fraction(1)
    for _ in range(10000):
        tR = _ivq(radius+T)
        tail = (rexp(-mu*_ivq(T))*(c_norm/_power(tR,m)+d_norm/_power(tR,m+1))/mu).b
        if tail <= budget:
            break
        T *= 2
    else:
        raise ArithmeticError('Could not bound infinite quadrature tail')
    total, local = iv.mpc(0), iv.mpf(0)
    left, panels, max_order = Fraction(0), 0, 0
    # This cap is an efficiency choice only; the full Cauchy majorant is
    # evaluated regardless of it.  Larger disks reduce the panel count.
    frequency_limit = Fraction(2)/_upper_fraction(omega_norm)
    while left < T:
        distance0 = max(anchor,radius+left)
        distance1 = max(anchor-1,radius+left)
        disk_radius = _dyadic_below(min(distance0/2,distance1/2,frequency_limit))
        half_width = min(disk_radius/8,(T-left)/2)
        panel_budget = budget*_ivq(2*half_width/T)
        value, err, used_order = _panel(omega,c,d,m,anchor,radius,left,half_width,
                                       disk_radius,mu,omega_norm,c_norm,d_norm,
                                       order,panel_budget)
        total += value
        local += err
        left += 2*half_width
        panels += 1
        max_order = max(max_order,used_order)
        if panels >= max_panels and left < T:
            raise ArithmeticError('Quadrature max_panels exceeded')
    analytic_error = (local+tail).b
    result = total+_error_box(analytic_error)
    diagnostics = dict(panels=panels,max_order=max_order,local_remainder=local.b,
                       tail_remainder=tail,analytic_error=analytic_error,
                       truncation_endpoint=T,reference_scale=scale,target=target,
                       phase_decay=mu,precision=iv.dps,
                       finite_taylor_sum=total,
                       branch='principal Log((a-1+i(R+t))/(a-1+iR)), t>=0')
    return (I0,result,diagnostics) if return_diagnostics else (I0,result)
