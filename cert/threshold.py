"""Certify a zero-free punctured neighborhood of the monopole threshold.

At z=-R the outgoing solution is written u_+=exp(i Omega x)m.
Both spatial continuations use the negative real direction on their far tail.
The Volterra bounds for Im(Omega)<=0 apply to the full closed frequency disk
intersected with that half-plane. The forward horizon state is propagated on
an enclosing frequency rectangle, independently on each finite spatial path.
Run from the project root: python -m cert.threshold
"""
import json
from pathlib import Path

from mpmath import iv

from .amplitude import error_box
from .elementary import backend_info, rexp, rlog
from .transport import (_segment, ball_rectangle, exact_box, horizon_enclosure,
                        precision)


def enclose_threshold(radius_denominator=1024, terminal_radius=2, dps=80):
    if not isinstance(radius_denominator, int) or not isinstance(terminal_radius, int):
        raise ValueError('The denominator and terminal radius must be integers')
    if not radius_denominator > terminal_radius > 0:
        raise ValueError('Require zero < R delta < 1')
    with precision(dps):
        delta = iv.mpf(1)/radius_denominator
        R = iv.mpf(terminal_radius)
        omega = iv.mpc(iv.mpf([-delta,delta]), iv.mpf([-delta,0]))
        moment = 1/(2*R)
        mass = 1/(2*R**2)
        defect = delta*(rlog(1/(R*delta))+2)
        dm = rexp(2*moment)*defect
        dv = mass*dm + 2*delta*moment*rexp(moment)
        log_ratio = rlog(1+1/R)
        m = R*log_ratio + error_box(dm)
        mx = 1/R-(1+1/R)*log_ratio + error_box(dv)
        initial = horizon_enclosure(omega,0,dps=dps)
        cases = []
        for name, sign in (('upper_spatial',1),('lower_spatial',-1)):
            source = iv.mpc('1.5',0)
            vertices = (source,iv.mpc('1.5',sign),iv.mpc(-R,sign),iv.mpc(-R,0))
            centers,radii = initial['centers'],initial['radii']
            records = []
            for start,end in zip(vertices,vertices[1:]):
                centers,radii = _segment(start,end,omega,0,centers,radii,96,'0.16',records)
            h,v = [ball_rectangle(c,r) for c,r in zip(centers,radii)]
            W = (2*iv.j*omega*h-v)*m+h*mx
            if not W.real.b < 0:
                raise ArithmeticError('The Wronskian enclosure failed to exclude zero')
            cases.append(dict(sheet=name,steps=len(records),h=exact_box(h),hx=exact_box(v),
                              wronskian=exact_box(W)))
        return dict(ell=0,frequency_radius=f'1/{radius_denominator}',ray_start=-terminal_radius,
                    domain='0 < abs(Omega) <= radius and Im(Omega) <= 0; respective physical bank',
                    working_decimal_digits=dps,arithmetic=backend_info(),
                    m=exact_box(m),mx=exact_box(mx),cases=cases)


def main():
    result = enclose_threshold()
    output = Path(__file__).with_name('results_threshold.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    from fractions import Fraction
    for row in result['cases']:
        print(row['sheet'],row['steps'],'W enclosure',
              {key:[float(Fraction(x)) for x in value]
               for key,value in row['wronskian'].items()})
    print(output)


if __name__=='__main__':
    main()
