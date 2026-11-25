"""Independent normalization and contour tests for the lateral enclosures."""
from fractions import Fraction as Q
from pathlib import Path
import json
from mpmath import iv
from .amplitude import rational
from .elementary import cabs, cexp
from .lateral import incoming,regularized_seed
from .transport import horizon_enclosure,ball_rectangle,precision


def intersects(a,b):
    return all(max(x.a,y.a)<=min(x.b,y.b)
               for x,y in ((a.real,b.real),(a.imag,b.imag)))


def unpack(z):
    def component(ab):
        return iv.mpf([rational(ab[0]).a,rational(ab[1]).b])
    return iv.mpc(component(z['real']),component(z['imag']))


def main():
    recorded=json.loads(Path(__file__).with_name('results_lateral.json').read_text())
    with precision(85):
        # At the first resonance, the residue line equals the unit outgoing
        # horizon column. This seed is independent of the regularized recurrence.
        q=rational('0.5')
        seed=horizon_enclosure(iv.j*q,0,dps=85)
        regular=regularized_seed(-iv.j*q,0)
        H,v=seed['h'],seed['hx']
        # e^(x-1) at z=3/2, normalized to leading coefficient +1 at z=1.
        factor=cexp(iv.mpc(q))*q
        targetH,targetv=factor*H,factor*(H+v)
        for i,target in enumerate((targetH,targetv)):
            direct=ball_rectangle(regular['centers'][i],regular['radii'][i])
            assert intersects(direct,target)
        for value,dec in zip((Q(1,4),Q(1,2)),('0.25','0.5'),strict=True):
            point,box=[r for r in recorded['tested'] if r['frequency_center_imag']=='-'+dec]
            assert box['box_halfwidth']=='1e-12' and point['box_halfwidth']=='0'
            min_mod=Q(322,100) if value==Q(1,4) else Q(104,100)
            assert Q(box['min_lateral_modulus'][0])>min_mod
            for row in point['cases']:
                old=unpack(row['enclosure'])
                other=incoming(iv.mpc(0,-rational(value)),0,row['shore'],R=72,
                               degree=31,dps=85,order=96,height=3,regularized=value==Q(1,2))
                assert intersects(old,other['enclosure'])
                assert cabs(other['enclosure']).a>min_mod.numerator/rational(min_mod.denominator)
        print('Independent outgoing seed, both modulus bounds, and four contour perturbations passed.')

if __name__=='__main__':main()
