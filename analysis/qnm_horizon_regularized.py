#!/usr/bin/env python3
"""Gamma-regularized scalar Schwarzschild horizon data.

Conventions: x=z+log(z-1), Omega=2*M*omega, time dependence exp(-i*omega*t),
L=ell*(ell+1), u_H/exp(-i*Omega*x)->1 as z->1+.

Returns u_hat=u_H/Gamma(1-2*i*Omega) and d u_hat/dx at an exterior matching
point (default z=1.5). The local series uses t=1-z and requires |t|<1.
It remains evaluable at Omega=-i*N/2, N=1,2,...; no singular Heun evaluation
or division by n-2*i*Omega is performed. mpmath is the only dependency.

The stopping and roundoff checks are numerical diagnostics, NOT certified
interval error bounds. Increase dps and compare results for production work.
Complex z uses mpmath's principal log(z-1); continuation to another spatial
sheet must be performed separately and explicitly.

Examples:
    from qnm_horizon_regularized import horizon_regularized
    u, ux = horizon_regularized('0.7-0.2j', ell=2, dps=70)
    u, ux = horizon_regularized(-1j/2, ell=0, dps=70)
    python qnm_horizon_regularized.py --test --output results.json
    python qnm_horizon_regularized.py --omega-re 0.7 --omega-im -0.2 --ell 2
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any
import mpmath as mp


@dataclass(frozen=True)
class HorizonData:
    u: Any
    ux: Any
    omega: Any
    z: Any
    ell: int
    terms: int
    working_dps: int
    last_term_scaled: Any
    roundoff_indicator: Any
    ode_residual_scaled: Any
    exact_resonance_n: int | None


def _complex(value):
    if isinstance(value, str):
        # mpmath's two-component parser does not accept Python complex strings.
        s = value.replace(' ', '').replace('I', 'j').replace('i', 'j')
        if 'j' in s:
            s = s.removesuffix('j')
            split = [i for i in range(1, len(s))
                     if s[i] in '+-' and s[i-1] not in 'eE']
            if split:
                i = split[-1]
                return mp.mpc(mp.mpf(s[:i]), mp.mpf(s[i:]))
            return mp.mpc(0, mp.mpf(s or '1'))
    return mp.mpc(value)


def _valid_ell(ell):
    if isinstance(ell, bool) or int(ell) != ell or ell < 0:
        raise ValueError('ell must be a nonnegative integer')
    return int(ell)


def resonance_determinant(n: int, ell: int) -> int:
    """Exact integer det A_n(ell*(ell+1)); positive for n>=1, ell>=0."""
    ell = _valid_ell(ell)
    if isinstance(n, bool) or int(n) != n or n < 1:
        raise ValueError('n must be a positive integer')
    n = int(n)
    ll = ell*(ell+1)
    dm2, dm1 = 0, 1
    for k in range(n):
        determinant = (n*n+ll-k*(k+1))*dm1 - n*k*(n-k)**2*dm2
        dm2, dm1 = dm1, determinant
    return dm1


def _evaluate_once(omega, ell, z, dps, working_dps, max_terms):
    with mp.workdps(working_dps):
        omega, z = _complex(omega), _complex(z)
        t, lam, ll = 1-z, 2j*omega, ell*(ell+1)
        if not mp.isfinite(omega) or not mp.isfinite(z):
            raise ValueError('omega and z must be finite')
        if abs(t) >= 1 or z == 1 or z == 0:
            raise ValueError('local matching point must obey 0<|1-z|<1')
        exact_n = (int(lam.real) if lam.imag == 0 and lam.real >= 1
                   and lam.real == int(lam.real) else None)
        minimum_n = max(32, int(abs(lam))+16)
        if max_terms <= minimum_n+12:
            raise ValueError('max_terms too small for this frequency')
        # All three sums use t derivatives. ux is formed afterwards by chain rule.
        sums = [mp.rgamma(1-lam), mp.mpc(0), mp.mpc(0)]
        sum_abs = [abs(sums[0]), mp.mpf(0), mp.mpf(0)]
        pm2, pm1 = mp.mpc(0), mp.mpc(1)
        factorial, t_power = mp.mpf(1), mp.mpc(1)
        small_count, last_scaled = 0, mp.inf
        tol = mp.power(10, -dps-8)
        for n in range(1, max_terms):
            p = (n*(n-1)-lam*lam-ll)*pm1 - lam*(n-1)*(n-1-lam)**2*pm2
            factorial *= n
            t_power *= t
            term = p*t_power*mp.rgamma(n+1-lam)/factorial
            terms = [term, n*term/t, n*(n-1)*term/(t*t)]
            for k in range(3):
                sums[k] += terms[k]
                sum_abs[k] += abs(terms[k])
            last_scaled = max(abs(terms[k])/max(1, abs(sums[k])) for k in range(3))
            pm2, pm1 = pm1, p
            small_count = small_count+1 if n >= minimum_n and last_scaled < tol else 0
            if small_count >= 12:
                break
        else:
            raise ArithmeticError(f'series did not meet stopping criterion in {max_terms} terms')
        roundoff = max(mp.eps*sum_abs[k]/max(1, abs(sums[k])) for k in range(3))
        H, Ht, Htt = sums
        prefactor = mp.exp(-lam+lam*z/2)*z*mp.exp(-lam*mp.log(z-1)/2)
        logder = lam/2+1/z-lam/(2*(z-1))
        logder_prime = -1/(z*z)+lam/(2*(z-1)**2)
        f, fp = (z-1)/z, 1/(z*z)
        u = prefactor*H
        uz = prefactor*(logder*H-Ht)
        uzz = prefactor*((logder_prime+logder**2)*H-2*logder*Ht+Htt)
        ux = f*uz
        uxx = f*f*uzz+f*fp*uz
        potential = f*(ll/z**2+1/z**3)
        rhs = (omega**2-potential)*u
        ode_residual = abs(uxx+rhs)/max(1, abs(uxx), abs(rhs))
        return HorizonData(u=u, ux=ux, omega=omega, z=z, ell=ell, terms=n+1,
                           working_dps=working_dps, last_term_scaled=last_scaled,
                           roundoff_indicator=roundoff,
                           ode_residual_scaled=ode_residual,
                           exact_resonance_n=exact_n)


def evaluate_horizon_regularized(omega, ell=0, z='1.5', *, dps=70,
                                  max_terms=4000, guard_digits=30) -> HorizonData:
    """Evaluate regularized u and ux, with convergence diagnostics.

    Inputs given as decimal strings preserve their supplied precision. Python
    float/complex inputs already carry binary-double rounding. Three attempts
    with increased internal precision are allowed when cancellation is detected.
    The nonzero scalar resonances are detected only by exact equality; nearby
    frequencies are never snapped to the resonance grid.
    """
    ell = _valid_ell(ell)
    if dps < 20 or guard_digits < 10:
        raise ValueError('use dps>=20 and guard_digits>=10')
    work = int(dps+guard_digits)
    for attempt in range(3):
        data = _evaluate_once(omega, ell, z, int(dps), work, int(max_terms))
        with mp.workdps(work):
            target = mp.power(10, -dps-5)
            if data.roundoff_indicator < target and data.ode_residual_scaled < target:
                return data
            # More precision is warranted only when diagnostics actually fail.
            work += max(guard_digits, int(mp.ceil(mp.log10(max(
                data.roundoff_indicator, data.ode_residual_scaled)/target)))+15)
    raise ArithmeticError('precision/ODE residual checks failed after precision escalation')


def horizon_regularized(omega, ell=0, z='1.5', **kwargs):
    """Return just (u_H/Gamma(1-2iOmega), d_x[u_H/Gamma(1-2iOmega)])."""
    data = evaluate_horizon_regularized(omega, ell, z, **kwargs)
    return data.u, data.ux


def _ordinary_horizon(omega, ell, z, digits=100):
    """Independent nonsingular h_n recurrence, used only by the test suite."""
    with mp.workdps(digits):
        omega, z = _complex(omega), _complex(z)
        lam, t, ll = 2j*omega, 1-z, ell*(ell+1)
        hm2, hm1 = mp.mpc(0), mp.mpc(1)
        H, Ht = hm1, mp.mpc(0)
        tpow = mp.mpc(1)
        for n in range(1, 1600):
            h = ((n*(n-1)-lam**2-ll)*hm1-lam*(n-1-lam)*hm2)/(n*(n-lam))
            tpow *= t
            term = h*tpow
            H += term
            Ht += n*term/t
            hm2, hm1 = hm1, h
            if n > 200 and max(abs(term), abs(n*term/t)) < mp.power(10, -digits-5):
                break
        prefactor = mp.exp(-lam+lam*z/2)*z*mp.exp(-lam*mp.log(z-1)/2)
        a = lam/2+1/z-lam/(2*(z-1))
        return prefactor*H, (z-1)/z*prefactor*(a*H-Ht)


def run_tests():
    """Regression tests, not a validation of any QNM frequency or infinity sheet."""
    rows=[]
    with mp.workdps(120):
        for n in [1,2,3,4,8,12,16,24,32]:
            for ell in [0,1,2,5]:
                data=evaluate_horizon_regularized(-mp.j*n/2,ell,dps=80)
                du,dv=_ordinary_horizon(mp.j*n/2,ell,mp.mpf('1.5'),digits=120)
                factor=mp.exp(-n)*resonance_determinant(n,ell)/mp.factorial(n)
                err=max(abs(data.u-factor*du)/max(1,abs(factor*du)),
                        abs(data.ux-factor*dv)/max(1,abs(factor*dv)))
                assert err < mp.mpf('1e-75'), (n,ell,err)
                assert abs(data.u)+abs(data.ux)>0
                assert data.exact_resonance_n == n
                rows.append({'N':n,'ell':ell,'determinant':str(resonance_determinant(n,ell)),
                             'identity_error':mp.nstr(err,8),'ode_residual':mp.nstr(data.ode_residual_scaled,8),
                             'terms':data.terms,'working_dps':data.working_dps})
        # Independent ordinary recurrence, generic complex frequencies.
        generic=[]
        for wr,wi in [('0.7','-0.2'),('2.0','-1.2'),('-0.4','-0.3'),('1.2','0.8')]:
            w=mp.mpc(wr,wi)
            for ell in [0,2]:
                data=evaluate_horizon_regularized(w,ell,dps=80)
                a,b=_ordinary_horizon(w,ell,mp.mpf('1.5'),digits=120)
                scale=mp.rgamma(1-2j*w)
                err=max(abs(data.u-scale*a)/max(1,abs(scale*a)),
                        abs(data.ux-scale*b)/max(1,abs(scale*b)))
                assert err < mp.mpf('1e-75')
                generic.append({'omega':[wr,wi],'ell':ell,'error':mp.nstr(err,8)})
        # Near resonance: no snapping, and no multiply-infinity-by-zero operation.
        near=[]
        for n in [1,4,12]:
            w=-mp.j*n/2+mp.mpf('1e-30')
            data=evaluate_horizon_regularized(w,2,dps=80)
            a,b=_ordinary_horizon(w,2,mp.mpf('1.5'),digits=120)
            scale=mp.rgamma(1-2j*w)
            err=max(abs(data.u-scale*a)/max(1,abs(scale*a)),
                    abs(data.ux-scale*b)/max(1,abs(scale*b)))
            assert err < mp.mpf('1e-75')
            assert data.exact_resonance_n is None
            near.append({'N':n,'offset':'1e-30','error':mp.nstr(err,8)})
        # Zero frequency is regular for the local horizon evaluator, although the
        # infinity coefficient formula 1/(2*i*Omega) needs separate analysis.
        static=[]
        for ell in range(6):
            data=evaluate_horizon_regularized(0,ell,dps=80)
            z=mp.mpf('1.5')
            fun=lambda r:r*mp.legendre(ell,2*r-1)
            eu,ev=fun(z),(z-1)/z*mp.diff(fun,z)
            err=max(abs(data.u-eu)/max(1,abs(eu)),abs(data.ux-ev)/max(1,abs(ev)))
            assert err < mp.mpf('1e-75')
            static.append({'ell':ell,'error':mp.nstr(err,8)})
    return {'status':'PASS','normalization':'u_H/Gamma(1-2iOmega)',
            'matching_z':'1.5','requested_dps':80,'certified_error_bound':False,
            'exceptional_identity_tests':rows,'generic_tests':generic,
            'near_resonance_tests':near,'static_legendre_tests':static}


def _json_data(data, digits):
    c=lambda value:{'real':mp.nstr(mp.re(value),digits),'imag':mp.nstr(mp.im(value),digits)}
    return {'u':c(data.u),'ux':c(data.ux),'omega':c(data.omega),'z':c(data.z),
            'ell':data.ell,'terms':data.terms,'working_dps':data.working_dps,
            'exact_resonance_n':data.exact_resonance_n,
            'last_term_scaled':mp.nstr(data.last_term_scaled,8),
            'roundoff_indicator':mp.nstr(data.roundoff_indicator,8),
            'ode_residual_scaled':mp.nstr(data.ode_residual_scaled,8),
            'certified_error_bound':False}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test',action='store_true')
    parser.add_argument('--output')
    parser.add_argument('--omega-re',default='0.7')
    parser.add_argument('--omega-im',default='-0.2')
    parser.add_argument('--ell',type=int,default=0)
    parser.add_argument('--z',default='1.5')
    parser.add_argument('--dps',type=int,default=70)
    args=parser.parse_args()
    if args.test:
        result=run_tests()
    else:
        with mp.workdps(args.dps+30):
            w=mp.mpc(args.omega_re,args.omega_im)
            result=_json_data(evaluate_horizon_regularized(w,args.ell,args.z,dps=args.dps),args.dps)
    payload=json.dumps(result,indent=2)+'\n'
    if args.output:
        Path(args.output).write_text(payload)
        print(json.dumps({'status':result.get('status','OK'),'output':str(Path(args.output).resolve())}))
    else:
        print(payload,end='')
