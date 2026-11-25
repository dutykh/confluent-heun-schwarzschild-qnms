"""Rigorous uniform frequency enclosures via nine-node complex interpolation.

Run from the project root. The default command recomputes 72 additional
validated point enclosures. The --check-saved option reruns the rational
interpolation and all inequalities from the stored additional samples, the
three original samples per mode and an independently regenerated outer bound.
"""
from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from math import prod
from pathlib import Path

from .amplitude import amplitude, rational
from .campaign import DEST, MODES, PRECISION, SHIFT, exact_json
from .check import as_box
from .isolate import ComplexBox, RealBox, make_certificate, schwarzschild_outer_bound, sqrt_bounds

S = Fraction(1, 10_000)
D = Fraction(1, 100)
ROOT_RADIUS = Fraction(1, 10**28)
NODES = tuple(Fraction(k, 4) for k in range(-4, 5))
HERE = Path(__file__).resolve().parent
SAMPLES = HERE / 'neighborhood_samples.json'
OUTPUT = HERE / 'results_neighborhoods.json'
BITS = 256


def evaluate(task):
    """A rationally specified sample of the exact amplitude on the physical sheet."""
    from mpmath import iv
    mode, node = task
    ell, real, imag, anchor, height, degree = MODES[mode]
    real, imag = Fraction(real) + S*node, Fraction(imag)
    iv.dps = PRECISION
    omega = iv.mpc(rational(real), rational(imag))
    enclosure, info = amplitude(omega, ell, anchor, height, degree, corrected=False,
                                transport_options={'dps': PRECISION})
    return {'mode': mode, 'ell': ell, 'anchor': anchor, 'height': height, 'degree': degree,
            'node': str(node), 'omega_exact': [str(real), str(imag)],
            'enclosure': exact_json(enclosure), 'endpoint_h': exact_json(info['h']),
            'endpoint_hx': exact_json(info['hx']), 'tail_upper': exact_json(info['tail']),
            'steps': info['transport']['steps'], 'working_dps': PRECISION}


def mul_linear(coeffs, node):
    """Multiply an ascending coefficient vector by (t-node)."""
    result = [Fraction(0)]*(len(coeffs)+1)
    for j, value in enumerate(coeffs):
        result[j] -= node*value
        result[j+1] += value
    return result


def lagrange_weights():
    """Exact power-basis coefficients of each real-node Lagrange basis polynomial."""
    basis = []
    for node in NODES:
        polynomial = [Fraction(1)]
        denominator = Fraction(1)
        for other in NODES:
            if other != node:
                polynomial = mul_linear(polynomial, other)
                denominator *= node-other
        basis.append(tuple(value/denominator for value in polynomial))
    return basis


def interval_sum(weights, intervals):
    lower = upper = Fraction(0)
    for weight, interval in zip(weights, intervals, strict=True):
        endpoints = (weight*interval.lo, weight*interval.hi)
        lower += min(endpoints)
        upper += max(endpoints)
    return RealBox(lower, upper)


def coefficient_boxes(samples):
    boxes = [samples[t] for t in NODES]
    weights = lagrange_weights()
    return [ComplexBox(interval_sum(column, [p.re for p in boxes]),
                       interval_sum(column, [p.im for p in boxes]))
            for column in zip(*weights, strict=True)]


def interpolation_errors(M, inner_radius):
    """Analytic bounds for |G-P| and |G'-P'| in the scaled t variable."""
    outer = D/S
    inner = Fraction(inner_radius)
    if not (0 < inner < outer and max(NODES) < outer):
        raise ValueError('The neighborhoods and nodes must lie inside the outer circle')
    # Q(t) = t * product over positive nodes of (t*t-node*node).
    outer_Q_lower = outer * prod(outer**2 - t**2 for t in NODES if t > 0)
    Q = [Fraction(1)]
    for t in NODES:
        Q = mul_linear(Q, t)
    Q_upper = sum(abs(q)*inner**j for j, q in enumerate(Q))
    Qprime_upper = sum(j*abs(q)*inner**(j-1) for j, q in enumerate(Q) if j)
    H_upper = outer*M/(outer_Q_lower*(outer-inner))
    Hprime_upper = H_upper/(outer-inner)
    return Q_upper*H_upper, Qprime_upper*H_upper+Q_upper*Hprime_upper


def existing_inputs(mode):
    ell, real, imag, anchor, height, degree = MODES[mode]
    center = Fraction(real), Fraction(imag)
    boxes = {}
    for pos, shift in [('center',0),('plus',SHIFT),('minus',-SHIFT)]:
        data = json.loads((DEST/f'{mode}_{pos}.json').read_text())
        if (tuple(map(Fraction, data['omega_exact'])) != (center[0]+shift,center[1])
                or (data['mode'],data['position'],data['degree'],data['radius'],data['anchor'])
                != (mode,pos,degree,height,anchor)):
            raise ValueError('Original sample does not match expected frequency and settings')
        boxes[pos] = as_box(data['enclosure'])
    return center, boxes, Fraction(schwarzschild_outer_bound(center, D, ell, anchor)['outer_sup'])


def decimal_directed(value, places, *, upper):
    """A public outward bound with exactly 'places' digits after the point."""
    value = Fraction(value)
    scaled = value*10**places
    index = -((-scaled.numerator)//scaled.denominator) if upper else scaled.numerator//scaled.denominator
    sign = '-' if index < 0 else ''
    integer, tail = divmod(abs(index), 10**places)
    return f'{sign}{integer}.{tail:0{places}d}'


def certify(mode, rows):
    center, old_boxes, M = existing_inputs(mode)
    old = make_certificate(center, ROOT_RADIUS, D, M, SHIFT, old_boxes['center'],
                           old_boxes['plus'], old_boxes['minus'], bits=224)
    if not old['proved']:
        raise ArithmeticError('Original Rouché proof failed')
    values = {Fraction(0):old_boxes['center']}
    ell, _, _, anchor, height, degree = MODES[mode]
    for row in rows:
        node = Fraction(row['node'])
        if (node in values or node not in NODES or row['mode'] != mode
                or tuple(map(Fraction,row['omega_exact'])) != (center[0]+S*node,center[1])
                or (row['ell'],row['anchor'],row['height'],row['degree'])
                != (ell,anchor,height,degree)):
            raise ValueError('Unexpected or duplicated interpolation sample')
        values[node] = as_box(row['enclosure'])
    coeffs = coefficient_boxes(values)
    scale_one = coeffs[1].divide_positive(S)
    beta_re, beta_im = scale_one.midpoint
    blo, bhi = sqrt_bounds(beta_re**2+beta_im**2, BITS)
    center_abs = values[0].abs_upper(BITS)
    Ef, Edt = interpolation_errors(M,1)
    spread = scale_one.radius_upper(BITS) + sum(j*co.abs_upper(BITS)/S
                  for j, co in enumerate(coeffs) if j>=2) + Edt/S
    derivative_lower = blo-spread
    amplitude_upper = center_abs+(bhi+spread)*S
    boundary_inner = derivative_lower*ROOT_RADIUS-center_abs
    boundary_outer = derivative_lower*S-center_abs
    original_gap = Fraction(old['inequality']['strict_gap']['exact'])
    contraction = spread/blo if blo>0 else Fraction(1)
    success = contraction<1 and boundary_inner>0 and boundary_outer>0
    if not success:
        raise ArithmeticError(f'Uniform test failed at mode {mode}')
    return {
        'mode':mode,'center':list(map(str,center)),'outer_radius':str(D),'outer_sup':str(M),
        'fixed_radius':str(S),'root_radius':str(ROOT_RADIUS),'scaled_nodes':list(map(str,NODES)),
        'power_coefficient_boxes':[co.as_json() for co in coeffs],
        'scale_in_polynomial':'t=(Omega-c)/fixed_radius',
        'interpolation_F_error_upper':str(Ef),'interpolation_dOmega_error_upper':str(Edt/S),
        'beta':list(map(str,(beta_re,beta_im))),'beta_modulus_lower':str(blo),
        'beta_modulus_upper':str(bhi),'derivative_disk_radius':str(spread),
        'derivative_modulus_lower':str(derivative_lower),'F_center_modulus_upper':str(center_abs),
        'F_neighborhood_modulus_upper':str(amplitude_upper),'contraction_upper':str(contraction),
        'inner_boundary_lower_via_uniform':str(boundary_inner),
        'inner_boundary_lower_via_original':str(original_gap),
        'inner_boundary_lower':str(max(original_gap,boundary_inner)),
        'outer_neighborhood_boundary_lower':str(boundary_outer),'proved':success}


def run(check_saved=False, workers=1):
    if workers<1:
        raise ValueError('A positive worker count is required')
    if check_saved:
        rows = json.loads(SAMPLES.read_text())
    else:
        tasks = [(mode,t) for mode in sorted(MODES) for t in NODES if t]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(evaluate,tasks))
        SAMPLES.write_text(json.dumps(rows,indent=2)+'\n')
    if len(rows) != len(MODES)*(len(NODES)-1):
        raise ValueError('Wrong number of interpolation samples')
    results = [certify(mode,[r for r in rows if r['mode']==mode]) for mode in sorted(MODES)]
    if check_saved:
        if json.loads(OUTPUT.read_text()) != results:
            raise ArithmeticError('Saved neighborhood bounds differ from exact replay')
    else:
        OUTPUT.write_text(json.dumps(results, indent=2)+'\n')
    for record in results:
        r = {k:Fraction(value) for k,value in record.items() if k in
             ('derivative_disk_radius','derivative_modulus_lower','F_neighborhood_modulus_upper',
              'inner_boundary_lower','contraction_upper','outer_neighborhood_boundary_lower')}
        # All display bounds are produced by integer division with known direction.
        re,im = map(Fraction,record['beta'])
        mid_re,mid_im = Fraction(round(re*10000),10000),Fraction(round(im*10000),10000)
        display_spread = r['derivative_disk_radius'] + sqrt_bounds((re-mid_re)**2+(im-mid_im)**2,BITS)[1]
        fmt = lambda x,k,up:decimal_directed(x,k,upper=up)
        print(record['mode'],fmt(mid_re,4,True),fmt(mid_im,4,True),fmt(display_spread,5,True),
              fmt(r['derivative_modulus_lower'],3,False),fmt(r['F_neighborhood_modulus_upper']*10000,3,True),
              fmt(r['inner_boundary_lower']*10**28,3,False),
              fmt(r['outer_neighborhood_boundary_lower']*10000,3,False),
              fmt(r['contraction_upper'],4,True),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-saved',action='store_true')
    parser.add_argument('--workers',type=int,default=1)
    args = parser.parse_args()
    run(args.check_saved,args.workers)
