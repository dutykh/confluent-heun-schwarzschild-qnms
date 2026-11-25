"""Recompute exact rational disk tests from the three saved coefficient boxes."""
import argparse
import json
from fractions import Fraction as F
from .campaign import DEST, MODES, SHIFT
from .isolate import ComplexBox, RealBox, make_certificate, schwarzschild_outer_bound


def as_box(box):
    return ComplexBox(RealBox(*box['real']), RealBox(*box['imag']))


def check(mode, root_radius='1e-28', *, check_saved=False):
    ell, real, imag, anchor, R, N = MODES[mode]
    center = (real, imag)
    boxes = {}
    for position, shift in [('center',0),('plus',SHIFT),('minus',-SHIFT)]:
        point = json.loads((DEST/f'{mode}_{position}.json').read_text())
        if (list(map(F,point['omega_exact'])) != [F(real)+shift,F(imag)]
            or (point['mode'],point['position'],point['degree'],point['radius'],point['anchor'])
            != (mode,position,N,R,anchor)):
            raise ValueError('Sample does not match specified mode, frequency or settings')
        boxes[position] = as_box(point['enclosure'])
    outer = schwarzschild_outer_bound(center,'0.01',ell,anchor)
    result = make_certificate(center,root_radius,'0.01',outer['outer_sup'],SHIFT,
                              boxes['center'], boxes['plus'], boxes['minus'], bits=224)
    result['outer_proof'] = outer
    path = DEST/f'{mode}_disk.json'
    if check_saved:
        if json.loads(path.read_text()) != result:
            raise ArithmeticError(f'Saved disk certificate differs at {mode}')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(mode,'certified',result['proved'], 'M',outer['outer_sup'], 'margin',
          result['inequality']['strict_gap']['decimal_hint'], flush=True)
    if not result['proved']: raise ArithmeticError('Strict inequality failed')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=MODES)
    parser.add_argument('--radius',default='1e-28')
    args=parser.parse_args()
    check(args.mode,args.radius)
