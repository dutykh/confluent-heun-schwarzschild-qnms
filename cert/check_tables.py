"""Read the release TeX and check the numerical entries in its eight tables."""
import argparse
from decimal import Decimal
from fractions import Fraction as F
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TEX = ROOT / 'reports/Confluent_Heun_Generic_Progress.tex'


def load(stem):
    return json.loads((ROOT/(stem+'.json')).read_text(), parse_float=Decimal)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def table_panels(tex, label):
    tables = re.findall(r'\\begin\{table\*?\}(.*?)\\end\{table\*?\}', tex, re.S)
    matches = [t for t in tables if r'\label{tab:'+label+'}' in t]
    require(len(matches) == 1, f'Expected exactly one table labelled {label}')
    panels = re.findall(r'\\begin\{tabular\}\{[^}]+\}(.*?)\\end\{tabular\}',
                        matches[0], re.S)
    result = []
    for panel in panels:
        clean = re.sub(r'\\(?:hline|toprule|midrule|bottomrule)\b', '', panel)
        rows = [r for r in clean.split(r'\\') if '&' in r]
        result.append([[c.strip().strip('$').strip() for c in row.split('&')]
                       for row in rows[1:]])
    return result


def scalar(cell):
    """Return the number and the decimal unit in its last printed place."""
    text = cell.replace(r'\,', '').replace(r'\times', '').strip()
    m = re.fullmatch(r'(.*)10\^\{(-?\d+)\}', text)
    exponent = 0
    if m:
        text, exponent = m.group(1) or '1', int(m.group(2))
    value = Decimal(text)
    return F(value)*F(10)**exponent, F(10)**(value.as_tuple().exponent+exponent)


def rounded(cell, source, direction='nearest', *, places=None, significant=None):
    value, unit = scalar(cell)
    if places is not None:
        require(unit == F(10)**-places, f'{cell}: expected {places} decimal places')
    if significant is not None:
        require(10**(significant-1) <= abs(value)/unit < 10**significant,
                f'{cell}: expected {significant} significant digits')
    difference = value - F(source)
    if direction == 'up':
        valid = 0 <= difference < unit
    elif direction == 'down':
        valid = -unit < difference <= 0
    else:
        valid = value == round(F(source)/unit)*unit
    require(valid, f'{cell} is not {direction} rounding of {source}')


def pair(text):
    return tuple(map(int, text.strip('()').split(',')))


def settings(rec):
    return rec['radius'], rec['asympt_order']


def mode_rows(stem):
    return {(r['ell'], r['n']):r for r in load('analysis/'+stem)['modes']}


def accessory(panels):
    rows, = panels
    cases = {F(c['epsilon']):c for c in load('cert/results_resonant')['brackets']}
    require(len(rows) == len(cases) == 4 and {scalar(r[0])[0] for r in rows} == cases.keys(),
            'Accessory rows differ')
    for eps, midpoint in rows:
        record = cases[scalar(eps)[0]]
        require(F(midpoint) == F(record['midpoint_decimal']), 'Wrong accessory midpoint')
        radius = F('1e-28')
        require(F(record['radius']) == radius, 'Wrong accessory radius')
        require((F(record['t_lower']), F(record['t_upper'])) ==
                (F(midpoint)-radius, F(midpoint)+radius), 'Wrong accessory bracket')
        for side in ('left', 'right'):
            endpoint = record[side]
            low, high = map(F, endpoint['complete_enclosure'])
            finite_low, finite_high = map(F, endpoint['finite_part'])
            tail = F(endpoint['tail_bound'])
            require(tail >= 0 and low <= finite_low-tail <= finite_high+tail <= high,
                    'Accessory interval omits transport or tail')
        require(F(record['left']['complete_enclosure'][0]) > 0 >
                F(record['right']['complete_enclosure'][1]), 'Unresolved endpoint signs')


def generic(panels):
    parameters, bounds = panels
    cases = {r['name'][-1]:r for r in load('cert/results_generic')['cases']}
    observed = {}
    name = None
    for label, symbols, values in parameters:
        if label:
            name = label
            require(name not in observed, 'Duplicate generic parameter group')
            observed[name] = []
        expected_symbols = (r'(\alpha,\beta,\gamma)' if not observed[name]
                            else r'(\delta,\eta)')
        require(symbols == expected_symbols, 'Wrong parameter labels for '+name)
        components = values.strip('()').replace('\\', '').replace(' ', '').split(',')
        for term in components:
            m = re.fullmatch(r'([+-]?[\d.]+)([+-][\d.]+)i', term)
            require(m, f'Invalid complex number {term}')
            observed[name].append(tuple(map(F,m.groups())))
    require(observed.keys() == cases.keys(), 'Missing generic parameter group')
    for name, case in cases.items():
        require(observed[name] == [tuple(map(F,z)) for z in case['parameters']],
                'Generic parameters differ: '+name)
    seen = set()
    for name, symbol, lower, upper in bounds:
        axis = {r'\Re\Bcal':'real', r'\Im\Bcal':'imag'}[symbol]
        entry = cases[name]['coefficient'][axis]
        lo, hi = F(entry['lower']), F(entry['upper'])
        certificate = cases[name]['coefficient_certificate']
        finite = certificate['finite_part'][axis]
        tail = F(certificate['tail_bound']['upper'])
        require(tail >= 0 and lo <= F(finite['lower'])-tail <= F(finite['upper'])+tail <= hi,
                f'Generic coefficient omits transport or tail: {name} {axis}')
        places = {'A':30, 'B':25}[name]
        rounded(lower,lo,'down',places=places)
        rounded(upper,hi,'up',places=places)
        seen.add((name,axis))
    require(len(bounds) == len(seen) == 4, 'Missing generic enclosure row')


def historical(tex, table):
    rows, = table_panels(tex,table)
    if table == 'adaptive':
        data = mode_rows('qnm_tail_optimization_results')
        previous, improved = None, None
        keys = ('previous_fixed','legacy_adaptive','adaptive_fixed_domain')
        for mode, *columns in rows:
            record = data[pair(mode)]
            times = columns[-1].split(r'\to')
            require(len(times) == len(keys), 'Missing an adaptive timing')
            for col, timing, key in zip(columns[:-1],times,keys,strict=True):
                require(pair(col) == settings(record['evaluations'][key]['settings']),
                        'Adaptive settings differ')
                rounded(timing,record['root_searches'][key]['wall_seconds'],places=2)
            if pair(mode) == (2,0):
                previous,improved = map(F,times[1:])
        for expression, expected in [(r'1-8.74/9.87=',100*(1-improved/previous)),
                                      (r'9.87/8.74-1=',100*(previous/improved-1))]:
            found = re.search(re.escape(expression)+r'([\d.]+)\\%',tex)
            require(found, 'Missing historical timing comparison '+expression)
            rounded(found.group(1),expected,places=2)
    elif table == 'forward':
        data = mode_rows('qnm_forward_coefficient_results')
        for mode, radius, degree, timing, discrepancy, tail in rows:
            record = data[pair(mode)]
            run = record['root_searches'][record['selected_root_run_index']]
            require((int(radius),int(degree)) == settings(run['settings']),
                    'Forward settings differ')
            rounded(timing,run['wall_seconds'],places=2)
            rounded(discrepancy,run['difference_from_CF'],significant=3)
            rounded(tail,run['amplitude_tail_bound_conditional'],'up',significant=3)
    else:
        records = mode_rows('qnm_free_correction_results')[(0,0)]['radius_diagnostics']
        data = {settings(r['settings']):r for r in records}
        for grid, uncorrected, corrected, tail, gain in rows:
            run = data[pair(grid)]
            rounded(uncorrected,run['uncorrected_discrepancy_from_reference'],significant=3)
            rounded(corrected,run['corrected_discrepancy_from_reference'],significant=3)
            rounded(tail,run['remaining_tail_bound_conditional'],'up',significant=3)
            rounded(gain,run['bound_gain_vs_best_uncorrected'],places=0)
    require(len(rows) == 3 and {pair(row[0]) for row in rows} == data.keys(),
            'Historical rows incomplete: '+table)


def certified(tex, table):
    rows, = table_panels(tex,table)
    neighborhoods = {r['mode']:r for r in load('cert/results_neighborhoods')}
    require(len(rows) == 9 and {r[0]+r[1] for r in rows} == neighborhoods.keys(),
            'Certified modes incomplete: '+table)
    for ell, overtone, *columns in rows:
        mode = ell+overtone
        record = neighborhoods[mode]
        if table == 'centers':
            point = load('cert/results/'+mode+'_center')
            disk = load('cert/results/'+mode+'_disk')
            require(point['mode'] == mode and point['position'] == 'center' and
                    disk['proved'] is True and record['proved'] is True,
                    'Missing successful certificate '+mode)
            proof = disk['outer_proof']
            require(disk['center'] == proof['center'] and proof['ell'] == int(ell) and
                    F(point['anchor']) == F(proof['anchor']) and
                    F(proof['outer_sup']) == F(disk['outer_sup']) == F(record['outer_sup']),
                    'Outer-disk inputs differ '+mode)
            require(list(map(F,columns)) == list(map(F,point['omega_exact'])) ==
                    [F(disk['center'][axis]) for axis in ('re','im')] == list(map(F,record['center'])),
                    'Incorrect center '+mode)
            require(F(disk['root_radius']) == F(record['root_radius']) == F('1e-28') and
                    F(proof['outer_radius']) == F(disk['outer_radius']) ==
                    F(record['outer_radius']) == F('1e-2') and F(record['fixed_radius']) == F('1e-4'),
                    'Incorrect disk radii '+mode)
            for short, long in [('re','real'),('im','imag')]:
                require(list(map(F,point['enclosure'][long])) ==
                        list(map(F,disk['sample_enclosures']['center'][short])) ==
                        list(map(F,record['power_coefficient_boxes'][0][short])),
                        'Inconsistent center sample '+mode)
        elif table == 'rouche':
            anchor, grid, outer, upper, lower = columns
            point = load('cert/results/'+mode+'_center')
            disk = load('cert/results/'+mode+'_disk')
            require(F(anchor) == F(point['anchor']) and
                    pair(grid) == (point['radius'],point['degree']) and
                    F(outer) == F(disk['outer_sup']), 'Incorrect setting '+mode)
            iq = disk['inequality']
            u = F(iq['left_side_upper']['exact'])*10**28
            v = F(iq['right_side_lower']['exact'])*10**28
            rounded(upper,u,'up',places=3)
            rounded(lower,v,'down',places=3)
            require(F(upper) < F(lower), 'Nonpositive disk margin '+mode)
        else:
            bre, bim, radius, dmin, fmax, inner, outer = columns
            x, y = map(F,record['beta'])
            rounded(bre,x,places=4)
            rounded(bim,y,places=4)
            display, unit = scalar(radius)
            excess = display - F(record['derivative_disk_radius'])
            shift_squared = (x-F(bre))**2 + (y-F(bim))**2
            # Characterize ceil(spread + |beta-b_star|) using rational comparisons.
            require(unit == F('1e-5') and excess >= 0 and excess**2 >= shift_squared and
                    (excess < unit or (excess-unit)**2 < shift_squared),
                    'Incorrect rounded derivative-disk radius '+mode)
            rounded(dmin,record['derivative_modulus_lower'],'down',places=3)
            rounded(fmax,10**4*F(record['F_neighborhood_modulus_upper']),'up',places=3)
            rounded(inner,10**28*F(record['inner_boundary_lower']),'down',places=3)
            rounded(outer,10**4*F(record['outer_neighborhood_boundary_lower']),'down',places=3)
            require(min(map(F,(dmin,inner,outer))) > 0, 'Nonpositive neighborhood bound '+mode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tex', type=Path, default=DEFAULT_TEX)
    args = parser.parse_args()
    tex = re.sub(r'(?<!\\)%[^\n]*', '', args.tex.read_text())
    tables = re.findall(r'\\begin\{table\*?\}(.*?)\\end\{table\*?\}',tex,re.S)
    labels = [re.findall(r'\\label\{tab:([^}]+)\}',t) for t in tables]
    expected = ('accessory','generic','adaptive','forward','correction','centers','rouche','wide')
    require(labels == [[name] for name in expected], 'Unexpected table list, order or labels')
    accessory(table_panels(tex,'accessory'))
    generic(table_panels(tex,'generic'))
    for name in ('adaptive','forward','correction'):
        historical(tex,name)
    for name in ('centers','rouche','wide'):
        certified(tex,name)
    match = re.search(r'At the three radii respectively we obtain(.*?)\\\]', tex, re.S)
    require(match, 'Missing fixed-endpoint ratio bounds')
    pairs = re.findall(r'\(([\d.]+),\s*([\d.]+)\)', match.group(1))
    gains = load('cert/results/00_fixed_endpoint_errors')
    require(len(pairs) == len(gains) == 3, 'Missing fixed-endpoint ratio')
    for (lo,hi), record in zip(pairs,gains,strict=True):
        actual_lo,actual_hi = map(F,record['improvement'])
        rounded(lo,actual_lo,'down',places=2)
        rounded(hi,actual_hi,'up',places=2)
        require(F(lo) < actual_lo <= actual_hi < F(hi), 'Incorrect actual-error ratio')
    print('PASS: all eight manuscript tables, including both generic panels, '
          'error directions, rounded slope centers and both historical timing ratios.')


if __name__ == '__main__':
    main()
