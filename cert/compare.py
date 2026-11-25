"""Enclose true approximation errors at one common rational frequency."""
import json
from fractions import Fraction as F
from .campaign import DEST
from .check import as_box
from .isolate import sqrt_bounds

def modulus(box):
    def nearest(part):
        return max(F(0),part.lo,-part.hi)
    return (sqrt_bounds(nearest(box.re)**2+nearest(box.im)**2,256)[0],box.abs_upper(256))

def run(*, check_saved=False):
    reference=json.loads((DEST/'00_center.json').read_text())
    exact=as_box(reference['enclosure'])
    records=[]
    for r in (20,30,40):
        source=json.loads((DEST/f'00_compare{r}.json').read_text())
        if source['omega_exact']!=reference['omega_exact']: raise ValueError('Mismatched frequencies')
        info=source['diagnostics']
        u=modulus(as_box(info['plain_finite_part'])-exact)
        v=modulus(as_box(info['finite_part'])-exact)
        if v[0]<=0: raise ArithmeticError('Cannot bound improvement quotient')
        gain=(u[0]/v[1],u[1]/v[0])
        records.append(dict(radius=r,degree=source['degree'],
                            uncorrected_error=list(map(str,u)),corrected_error=list(map(str,v)),
                            improvement=list(map(str,gain))))
        print(r,[(float(a),float(b)) for a,b in (u,v,gain)])
    path = DEST/'00_fixed_endpoint_errors.json'
    if check_saved:
        if json.loads(path.read_text()) != records:
            raise ArithmeticError('Saved fixed endpoint error comparison differs')
    else:
        path.write_text(json.dumps(records, indent=2)+'\n')

if __name__=='__main__': run()
