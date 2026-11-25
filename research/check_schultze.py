#!/usr/bin/env python3
"""Schultze (2006), (6.37), (6.41), (6.42), applied to generic case B.

Run: python research/check_schultze.py --output research/schultze_comparison.json
Check the saved comparison: python research/check_schultze.py --check-saved
Requires Python 3 and mpmath (version included in JSON). No ODE solver used.
The frozen intervals come from the independent residual/transport proof archive
cert/results_generic.json, complex_B, outward_decimal_enclosure_40_places.
"""
import argparse
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path

import mpmath as mp

PARAMETERS = [[".8", ".6"], ["-.2", ".15"], [".4", "-.25"],
              [".35", ".2"], ["-.1", ".4"]]
BOX = [
    ["0.5537762471318039728382188972524773924480",
     "0.5537762471318039728382471409990674875430"],
    ["0.8443876927568560589000227480816169483819",
     "0.8443876927568560589000509918282070434769"]]


def serialize(z):
    return [mp.nstr(z.real, 43), mp.nstr(z.imag, 43)]


def error_bounds(z):
    """Rigorous modulus lower/upper bounds for the *printed decimal pair*."""
    lo2, hi2 = (Fraction(0) for _ in range(2))
    for text, interval in zip(z, BOX):
        q = Fraction(text)
        a, b = map(Fraction, interval)
        lo2 += max(Fraction(0), a-q, q-b)**2
        hi2 += max(abs(a-q), abs(b-q))**2
    scale = 10**44
    low = isqrt((lo2.numerator*scale**2)//lo2.denominator)
    high = 1+isqrt((hi2.numerator*scale**2)//hi2.denominator)
    return [f"{low//scale}.{low%scale:044d}", f"{high//scale}.{high%scale:044d}"]


def extrapolate(points):
    """Neville value at h=0, for h=1/sqrt(k)."""
    xs = [1/mp.sqrt(k) for k, _ in points]
    vals = [value for _, value in points]
    for j in range(1, len(points)):
        vals = [(xs[i]*vals[i+1]-xs[i+j]*vals[i])/(xs[i]-xs[i+j])
                for i in range(len(vals)-1)]
    return vals[0]


def compute(dps, max_power):
    with mp.workdps(dps):
        g, beta, gamma, delta, eta = [mp.mpc(x, y) for x, y in PARAMETERS]
        A = delta + g*(beta+gamma+2)/2
        B = eta + beta/2 + (gamma-g)*(beta+1)/2
        a0, a1 = -gamma, -beta
        b0 = -(A+B) + (gamma+1)*(beta+1+g)/2
        b1 = B + (beta+1)*(g-gamma-1)/2
        r = (b0+b1)/g
        mu = (a0+a1)/2+r-1
        nu = (a0+a1)/2-r-1
        prefactor = 2*mp.sqrt(mp.pi)*mp.exp(g/2-mp.j*mp.pi*nu)*g**(-r-mp.mpf(".25"))
        assert abs(r+delta/g) < 100*mp.eps

        # Polynomial coefficients with a k-independent part removed only once.
        al = -2*(a1+mu)+a0+g
        ac = -b1-(1-a1)*(mu+(1-a0-g)/2)
        bl = 4*mu-2*a0+a1-g
        bc = (1-a0)*(1-a1+g)/2-mu**2+mu*(1+a0+g-a1)-b0
        cl = a0-2*mu
        cc = mu*(mu-a0)
        e0, em1, em2 = mp.mpc(1), mp.mpc(0), mp.mpc(0)
        samples = []
        rows = []
        for k in range(2**max_power):
            ka = 3*k*k+k*al+ac
            kb = -3*(k-1)**2+(k-1)*bl+bc
            kc = (k-2)**2+(k-2)*cl+cc
            n = k+1
            e = (ka*e0+kb*em1+kc*em2)/(n*(n-a1))
            if n == 1:
                # v(t)=(1-t)^mu H0(-t/(1-t)).
                assert abs(e+mu+B/(1+beta)) < 100*mp.eps
            if n >= 128 and n & (n-1) == 0:
                scale = mp.exp(-2*mp.sqrt(g)*mp.sqrt(n))*n**(r+mp.mpf(".75"))
                raw = prefactor*scale*e
                samples.append((n, raw))
                raw_text = serialize(raw)
                acc_text = serialize(extrapolate(samples))
                rows.append({"k": n, "raw": raw_text,
                             "raw_error_interval": error_bounds(raw_text),
                             "neville_all_nodes": acc_text,
                             "neville_error_interval": error_bounds(acc_text)})
            em2, em1, e0 = em1, e0, e
        return {"working_dps": dps, "reflected_parameters": [
                    serialize(z) for z in [a0, b0, a1, b1, g]],
                "r_mu_nu": [serialize(z) for z in [r, mu, nu]],
                "rows": rows}


def check_saved(saved_path, reference_path):
    """Read-only verification of archived decimal-distance bounds and a short rerun."""
    data = json.loads(saved_path.read_text(encoding="utf-8"))
    if data["parameters_alpha_beta_gamma_delta_eta"] != PARAMETERS:
        raise ValueError("Saved parameter set differs")
    source = json.loads(reference_path.read_text(encoding="utf-8"))
    target = next(row for row in source["cases"] if row["name"] == "complex_B")
    if target["parameters"] != PARAMETERS:
        raise ValueError("Reference parameter set differs")
    dims = ["real", "imag"]
    archived = [target["outward_decimal_enclosure_40_places"][j] for j in dims]
    if data["reference_real_imag_outward"] != archived or BOX != archived:
        raise ValueError("Reference rectangle differs")
    for idx, dim in enumerate(dims):
        exact = target["coefficient"][dim]
        if not (Fraction(BOX[idx][0]) <= Fraction(exact["lower"])
                <= Fraction(exact["upper"]) <= Fraction(BOX[idx][1])):
            raise ValueError("Printed rectangle fails to contain archived rational bounds")
    runs = data["runs"]
    if [r["working_dps"] for r in runs] != [58, 88]:
        raise ValueError("The release comparison requires both working precisions")
    low_rows, high_rows = [r["rows"] for r in runs]
    if low_rows != high_rows:
        raise ValueError("Saved raw or extrapolated 43-digit values differ between runs")
    if [r["k"] for r in high_rows] != [2**j for j in range(7, 18)]:
        raise ValueError("Stored dyadic sample is incomplete")
    for row in high_rows:
        for field, bound in [("raw", "raw_error_interval"),
                             ("neville_all_nodes", "neville_error_interval")]:
            if error_bounds(row[field]) != row[bound]:
                raise ValueError(f"Decimal comparator distance differs at k={row['k']}")
    if high_rows[:6] != compute(58, 12)["rows"]:
        raise ValueError("Short independent rerun failed")
    print("PASS: canonical rational rectangle, eleven dyadic samples, both precision "
          "runs, exact decimal distance bounds and recurrence replay through k=4096.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write full recomputation here (default: stdout)")
    parser.add_argument("--check-saved", action="store_true", help="read-only fast verification")
    parser.add_argument("--saved", type=Path, default=Path(__file__).with_name("schultze_comparison.json"))
    parser.add_argument("--reference", type=Path, default=Path(__file__).resolve().parent.parent/"cert"/"results_generic.json")
    parser.add_argument("--max-power", type=int, default=17)
    parser.add_argument("--precision", nargs="+", type=int, default=[58, 88])
    args = parser.parse_args()
    if args.check_saved:
        check_saved(args.saved, args.reference)
        return
    if not 7 <= args.max_power:
        parser.error("--max-power must be at least 7")
    data = {"source": "Stephan Schultze (Duisburg-Essen, 2006), pp. 99-100",
            "mpmath_version": mp.__version__,
            "parameters_alpha_beta_gamma_delta_eta": PARAMETERS,
            "reference_real_imag_outward": BOX,
            "extrapolation": "Neville at h=0 through all dyadic k=128,...,kmax in h=k^(-1/2). No intrinsic truncation bound.",
            "reference_error_interpretation":
                "Each pair bounds the modulus of the difference of the displayed decimal and any point in the archived certified box. The bounds use exact rational arithmetic and integer square roots.",
            "runs": [compute(p, args.max_power) for p in args.precision]}
    out = json.dumps(data, indent=2)+"\n"
    if args.output:
        args.output.write_text(out, encoding="utf-8")
    else:
        print(out, end="")


if __name__ == "__main__":
    main()
