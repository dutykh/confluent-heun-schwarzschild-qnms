"""Scalar Schwarzschild QNM incoming amplitude from a vertical outgoing contour.

Conventions: Omega=2 M omega, exp(-i omega t), z=r/(2M),
x=z+log(z-1). State (m,v)=(exp(-i Omega x)u_out, d m/dx).
Only the outgoing endpoint is initialized by inverse powers; the truncated
polynomial is NOT assumed convergent. Its two-term exact differential residual
gives the reported rigorous analytic terminal-error bounds. The SciPy ODE and
root errors are empirical, not interval-certified.

Run: python analysis/qnm_contour.py --campaign
Dependencies: numpy scipy mpmath (the latter only for optional high precision).
"""
import argparse
import cmath
import json
import math
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


def asymptotic_coefficients(omega, ell, order):
    L = ell * (ell + 1)
    a = [1.0 + 0j]
    for k in range(order):
        previous = a[k - 1] if k else 0j
        a.append(((k * (k + 1) - L) * a[k] - k * k * previous)
                 / (2j * omega * (k + 1)))
    return a


def terminal_state(omega, ell, anchor, radius, order):
    a = asymptotic_coefficients(omega, ell, order)
    z = anchor + 1j * radius
    p = sum(c * z ** (-n) for n, c in enumerate(a))
    pz = sum(-n * c * z ** (-n - 1) for n, c in enumerate(a) if n)
    return np.array([p, (z - 1) / z * pz], dtype=complex)


def analytic_tail_bounds(omega, ell, anchor, radius, order):
    """Absolute m,m_x errors for exact propagation; roundoff is not included."""
    p = omega.real
    q = max(0.0, -omega.imag)
    monotone = p > 0 and q * q <= 4 * p * p * anchor * (anchor - 1)
    if not monotone:
        raise ValueError('Anchor does not satisfy monotone phase condition')
    L = ell * (ell + 1)
    a = asymptotic_coefficients(omega, ell, order)
    N = order
    c = (N * (N + 1) - L) * a[N] - N * N * (a[N-1] if N else 0j)
    d = -(N + 1) ** 2 * a[N]
    F = abs(c) * radius ** (-N - 1) / (N + 1)
    F += abs(d) * radius ** (-N - 2) / (N + 2)
    V0 = math.pi * L / (2 * anchor) + 1 / anchor ** 2
    # atan(anchor/radius) avoids cancellation for a long tail.
    VR = L / anchor * math.atan(anchor / radius)
    VR += 1 / (math.sqrt(anchor * anchor + radius * radius)
               * (math.sqrt(anchor * anchor + radius * radius) + radius))
    E0 = F / abs(omega) * math.exp(VR / abs(omega))
    E1 = abs(omega) * E0
    bseg = (V0 - VR) / abs(omega)
    A = E0 + E1 / abs(omega)
    return dict(terminal_value=E0, terminal_derivative=E1,
                matching_value=A * math.exp(bseg),
                matching_derivative=E1 + abs(omega) * A * math.expm1(bseg),
                residual_integral=F, phase_monotone=True,
                total_potential_integral_bound=V0)


def horizon_state(omega, ell, anchor=1.5, terms=140, rtol=2e-12):
    """Unit horizon state (u_H, du_H/dx), from its convergent local series."""
    z0 = 1.5
    w = 1 - z0
    L = ell * (ell + 1)
    lam = 2j * omega
    h0 = 1.0 + 0j
    hprev = 0j
    H = h0
    Hz = 0j
    for n in range(1, terms + 1):
        den = n * (n - lam)
        if abs(den) < 1e-12:
            raise ValueError('Horizon resonance: use regularized recurrence')
        hn = ((n * (n - 1) + 4 * omega ** 2 - L) * h0
              - (lam * (n - 1) + 4 * omega ** 2) * hprev) / den
        H += hn * w ** n
        Hz -= n * hn * w ** (n - 1)
        hprev, h0 = h0, hn
    f = (z0 - 1) / z0
    g = cmath.exp(-2j * omega + 1j * omega * z0) * z0 * (z0 - 1) ** (-1j * omega)
    qg = 1j * omega + 1 / z0 - 1j * omega / (z0 - 1)
    state = np.array([g * H, g * f * (Hz + qg * H)], dtype=complex)
    if anchor != z0:
        def rhs(z, y):
            f = (z - 1) / z
            U = f * (L / z ** 2 + 1 / z ** 3)
            return np.array([y[1] / f, (U - omega ** 2) * y[0] / f])
        sol = solve_ivp(rhs, (z0, anchor), state, method='DOP853',
                        rtol=rtol, atol=rtol * 0.01)
        if not sol.success:
            raise RuntimeError(sol.message)
        state = sol.y[:, -1]
    return state


def incoming_amplitude(omega, ell, anchor=1.5, radius=80.0, order=40,
                       rtol=2e-12, horizon_terms=140, diagnostics=False):
    """The physical A_in, not just a scaled equation with the same zeros."""
    omega = complex(omega)
    if omega.real <= 0 or not (math.isfinite(omega.real) and math.isfinite(omega.imag)):
        raise ValueError('The implemented continuation requires Re(Omega)>0')
    if ell < 0 or int(ell) != ell or anchor <= 1 or radius <= 0 or order < 0:
        raise ValueError('Require integer ell>=0, anchor>1, radius>0, order>=0')
    L = ell * (ell + 1)
    terminal = terminal_state(omega, ell, anchor, radius, order)
    def rhs(s, y):
        z = anchor + 1j * s
        b = z / (z - 1)
        v = L / z ** 2 + 1 / z ** 3
        return 1j * np.array([b * y[1], v * y[0] - 2j * omega * b * y[1]])
    sol = solve_ivp(rhs, (radius, 0.0), terminal, method='DOP853',
                    rtol=rtol, atol=rtol * 0.01)
    if not sol.success:
        raise RuntimeError(sol.message)
    m, mx = sol.y[:, -1]
    uH, uHx = horizon_state(omega, ell, anchor, horizon_terms, rtol)
    x = anchor + math.log(anchor - 1)
    phase = cmath.exp(1j * omega * x)
    D = phase * (uH * mx + (1j * omega * uH - uHx) * m)
    A = D / (2j * omega)
    if diagnostics:
        bounds = analytic_tail_bounds(omega, ell, anchor, radius, order)
        bounds['determinant_tail_only'] = abs(phase) * (
            abs(uH) * bounds['matching_derivative']
            + abs(1j * omega * uH - uHx) * bounds['matching_value'])
        return A, dict(bounds=bounds, ode_evaluations=sol.nfev,
                       horizon_value=[uH.real,uH.imag],
                       outgoing_scaled_value=[m.real,m.imag],
                       outgoing_scaled_derivative=[mx.real,mx.imag])
    return A


CASES = [
    (0,0,.22-.21j,1.5),(0,1,.17-.70j,3.0),(0,2,.15-1.20j,5.0),
    (1,0,.59-.20j,1.5),(1,1,.53-.61j,1.5),(1,2,.46-1.08j,2.0),
    (2,0,.97-.19j,1.5),(2,1,.93-.59j,1.5),(2,2,.86-1.02j,1.5),
]


def solve_mode(ell, overtone, guess, anchor, radius, order, rtol):
    t0=time.perf_counter()
    def residual(vec):
        val=incoming_amplitude(complex(*vec),ell,anchor,radius,order,rtol)
        return [val.real,val.imag]
    answer=root(residual,[guess.real,guess.imag],tol=1e-10)
    w=complex(*answer.x)
    value,diag=incoming_amplitude(w,ell,anchor,radius,order,rtol,diagnostics=True)
    return dict(ell=ell,n=overtone,omega=[w.real,w.imag],anchor=anchor,
                radius=radius,order=order,rtol=rtol,amplitude_abs=abs(value),
                root_success=bool(answer.success),root_message=str(answer.message),
                elapsed_seconds=time.perf_counter()-t0,**diag)


def campaign():
    results=[]
    path=Path(__file__).with_name('qnm_contour_results.json')
    for ell,n,guess,anchor in CASES:
        record=dict(ell=ell,n=n,runs=[])
        for R,N,tol in [(60.,28,1e-10),(80.,40,2e-12),(100.,48,3e-14)]:
            result=solve_mode(ell,n,guess,anchor,R,N,tol)
            record['runs'].append(result)
            print(ell,n,R,N,result['omega'],result['amplitude_abs'],flush=True)
            guess=complex(*result['omega'])
        # Deliberately change the matching point/ray as a separate diagnostic.
        result=solve_mode(ell,n,guess,anchor+0.4,100.,48,3e-14)
        record['runs'].append(result)
        results.append(record)
        path.write_text(json.dumps(dict(convention='Omega=2Momega',
            certified=False,tail_bounds='analytic bounds; arithmetic and propagation empirical',
            modes=results),indent=2))
    return results


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--campaign',action='store_true')
    args=ap.parse_args()
    if args.campaign:
        campaign()
    else:
        print(incoming_amplitude(.22-.21j,0,diagnostics=True))
