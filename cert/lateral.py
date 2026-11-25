"""Two lateral incoming amplitudes on the negative imaginary frequency axis.

Run from the project root with python -m cert.lateral. The infinity tail runs
to the left above or below both regular singularities. The finite path from
z=3/2 and the logarithm continued from the exterior specify the lateral sheet.
"""
from argparse import ArgumentParser
from fractions import Fraction
import json
from pathlib import Path
from mpmath import iv
from .amplitude import error_box, norm_upper, polynomial, rational
from .campaign import exact_json
from .elementary import backend_info, cabs, cexp
from .transport import (_segment, ball_rectangle, ball_from_rectangle,
                        horizon_enclosure, precision, upper)


def regularized_seed(omega,ell,terms=220):
    """Enclose (1-2iΩ)h across the first horizon resonance at fixed z=3/2."""
    nu=2*iv.j*omega
    L=ell*(ell+1)
    w=-rational('0.5')
    s=[1-nu,-nu**2-L]
    for j in range(2,terms+1):
        denominator=j*(j-nu)
        if cabs(denominator).a <= 0:
            raise ValueError('Another horizon resonance meets the frequency box.')
        s.append(((j*(j-1)-nu**2-L)*s[-1]-nu*(j-1-nu)*s[-2])/denominator)
    H=s[-1]
    Hz=-terms*s[-1]
    for j in range(terms-1,-1,-1):
        H=H*w+s[j]
    for j in range(terms-1,0,-1):
        Hz=Hz*w-j*s[j]
    rho=rational('0.9')
    J=terms+1
    u=cabs(nu).b
    if not u<J:
        raise ValueError('Horizon series too short.')
    A=(1+(u*u+L)/J**2)/(1-u/J)
    B=(u/J+(u*u+L)/J**2)/(1-u/J)
    if not A*rho+B*rho**2 < 1:
        raise ValueError('Horizon majorant contraction failed.')
    C=max((cabs(s[-1])*rho**terms).b,(cabs(s[-2])*rho**(terms-1)).b)
    t=-w/rho
    E=C*t**(terms+1)/(1-t)
    Ez=C/(-w)*t**(terms+1)*((terms+1)/(1-t)+t/(1-t)**2)
    z=rational('1.5')
    f=(z-1)/z
    gauge=z*cexp(nu*(z-1))
    h=gauge*H
    v=f*gauge*(Hz+(nu+1/z)*H)
    c0,r0=ball_from_rectangle(h,upper(cabs(gauge)*E))
    c1,r1=ball_from_rectangle(v,upper(cabs(f*gauge)*(Ez+cabs(nu+1/z)*E)))
    return {'centers':(c0,c1),'radii':(r0,r1),'horizon_terms':terms}


def incoming(omega, ell, shore, R=80, degree=34, dps=85, order=96, height=2, regularized=False):
    if shore not in (-1, 1) or R < 2 or height <= 0:
        raise ValueError('Expected shore +/-1, R>=2 and positive height.')
    with precision(dps):
        nu = 2*iv.j*omega
        omega_lower = (-omega.imag).a
        phase = nu.real.a - cabs(nu).b / (R+1)
        if not (omega_lower > 0 and phase > 0):
            raise ValueError('The horizontal tail damping is not established.')
        seed = (regularized_seed(omega,ell) if regularized else
                horizon_enclosure(omega, ell, dps=dps))
        centers, radii = seed['centers'], seed['radii']
        z0 = iv.mpc(rational(Fraction(3,2)))
        corner = z0 + iv.j*(shore*height)
        target = iv.mpc(-R,shore*height)
        records = []
        for start, end in ((z0,corner),(corner,target)):
            centers,radii=_segment(start,end,omega,ell,centers,radii,order,'0.18',records)
        h,hx = [ball_rectangle(c,r) for c,r in zip(centers,radii,strict=True)]
        P,Pz,c,d = polynomial(omega,ell,degree,target)
        finite = ((nu*h-hx)*P + h*(target-1)/target*Pz)/nu
        k = 1/omega_lower
        V = rational(Fraction(ell*(ell+1),R) + Fraction(1,2*R**2))
        if not k*V < 1:
            raise ValueError('The rational exponential majorant needs kV < 1.')
        G = (norm_upper(c)/(R**(degree+1)*(degree+1))
             + norm_upper(d)/(R**(degree+2)*(degree+2)))
        free_b = hx/nu
        free_a = h-free_b
        tail = ((norm_upper(free_a)+norm_upper(free_b))*G/(2*omega_lower*(1-k*V))).b
        box = finite+error_box(tail)
        return {'first_resonance_regularized':regularized,'shore':shore,'path':[[1.5,0],[1.5,shore*height],[-R,shore*height]],
                'endpoint':target,'ell':ell,'degree':degree,'dps':dps,'order':order,
                'horizon_terms':seed['horizon_terms'],'horizon_radii':seed['radii'],
                'terminal_state_radii':radii,'steps':len(records),
                'minimum_damping':phase,'finite':finite,'tail':tail,'enclosure':box}


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='cert/results_lateral.json')
    args=parser.parse_args()
    with precision(85):
        output=[]
        for frequency,regularized in (('0.25',False),('0.5',True)):
            for halfwidth in ('0','1e-12'):
                radius=rational(halfwidth)
                offset=iv.mpf([-radius.b,radius.b])
                omega=iv.mpc(offset,-rational(frequency)+offset)
                cases=[incoming(omega,0,shore,regularized=regularized) for shore in (1,-1)]
                plus,minus=[x['enclosure'] for x in cases]
                jump=plus-minus
                if not(cabs(plus).a>0 and cabs(minus).a>0):
                    raise ArithmeticError('Lateral nonzero test failed.')
                if not (regularized or jump.imag.a>0 or jump.imag.b<0):
                    raise ArithmeticError('No separation between lateral limits established.')
                print(frequency,halfwidth,'right',plus,flush=True)
                print(frequency,halfwidth,'left ',minus,flush=True)
                output.append({'frequency_center_imag':'-'+frequency,
                          'box_halfwidth':halfwidth,'omega':omega,'cases':cases,
                          'jump':jump,'min_lateral_modulus':min(cabs(plus).a,cabs(minus).a)})
        result={'backend':backend_info(),'tested':output}
        Path(args.output).write_text(json.dumps(exact_json(result),indent=2)+'\n')

if __name__=='__main__':main()
