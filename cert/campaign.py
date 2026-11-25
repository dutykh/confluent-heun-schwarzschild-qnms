"""Recompute certified point amplitudes for the Schwarzschild test disks."""
import argparse
import json
import os
from fractions import Fraction
from pathlib import Path
from mpmath import iv
from .amplitude import amplitude, rational

# Centers are proposed decimal rationals. Their provenance is immaterial to the proof.
MODES = {
    '00': (0, '0.220909878160839371750923012336', '-0.209791434173761917563478133711', '1.5', 80, 46),
    '20': (2, '0.967287744421425973455558870304', '-0.193517551956575725317307466702', '1.5', 24, 44),
    '02': (0, '0.151483871070351525940443288239', '-1.202157180071607028187312037419', '5', 30, 70),
    '01': (0, '0.172233836672798344939925390787', '-0.696104893612920991050153640965', '3', 72, 80),
    '10': (1, '0.585872266534565417244491491108', '-0.195319977827156444311113978345', '1.5', 65, 72),
    '11': (1, '0.528897301209665079267177324478', '-0.612514783118094246451411572602', '1.5', 55, 75),
    '12': (1, '0.459078669862603343695496136573', '-1.080266850038214426940243559470', '2', 50, 78),
    '21': (2, '0.927701158039531126436544113024', '-0.591207873975925052419212578699', '1.5', 48, 70),
    '22': (2, '0.861088108753151536213124505866', '-1.017116804308548974940009868965', '1.5', 50, 75),
}
PRECISION = 85
SHIFT = Fraction(1, 10**12)
DEST = Path(__file__).resolve().parent / 'results'

def endpoint(raw):
    sign, mantissa, exponent, bits = raw
    if bits < 0:
        raise ArithmeticError('nonfinite interval')
    return Fraction((-1)**sign * mantissa) * Fraction(2)**exponent

def exact_json(obj):
    if hasattr(obj, '_mpi_'):
        return [str(endpoint(e)) for e in obj._mpi_]
    if hasattr(obj, '_mpci_'):
        return {'real':exact_json(obj.real), 'imag':exact_json(obj.imag)}
    if hasattr(obj, '_mpf_'):
        return str(endpoint(obj._mpf_))
    if hasattr(obj, '_mpc_'):
        return list(map(lambda raw: str(endpoint(raw)),obj._mpc_))
    if isinstance(obj, Fraction):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): exact_json(v) for k,v in obj.items()}
    if isinstance(obj, (tuple,list)):
        return [exact_json(v) for v in obj]
    if isinstance(obj,(int,str,bool,float)) or obj is None:
        return obj
    return str(obj)

def run(mode, position, comparison_radius=None):
    iv.dps=PRECISION
    ell, re, im, anchor, radius, degree = MODES[mode]
    if comparison_radius is not None:
        if mode!='00' or position!='center': raise ValueError('Fundamental center required')
        radius=comparison_radius
        degree={20:11,30:17,40:22}[radius]
    delta={'center':0,'plus':SHIFT,'minus':-SHIFT}[position]
    omega=iv.mpc(rational(Fraction(re)+delta),rational(im))
    corrected=position=='center'
    z, info=amplitude(omega,ell,anchor,radius,degree,corrected,
                      transport_options={'dps':PRECISION},
                      quadrature_options={'rtol':'1e-15'})
    if comparison_radius is not None:
        from .amplitude import from_endpoint
        _, plain = from_endpoint(omega,ell,anchor,radius,degree,info['h'], info['hx'],False)
        info['plain_finite_part']=plain['finite_part']
        info['plain_proof_tail']=plain['tail']
    from .elementary import backend_info
    result=dict(backend=backend_info(),mode=mode,position=position,omega_exact=[str(Fraction(re)+delta),im],
                precision=PRECISION,anchor=anchor,radius=radius,degree=degree,
                corrected=corrected,enclosure=z,diagnostics=info)
    DEST.mkdir(exist_ok=True)
    tag=f'compare{comparison_radius}' if comparison_radius else position
    path=DEST / f'{mode}_{tag}.json'
    tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(exact_json(result),indent=2)+'\n')
    os.replace(tmp,path)
    print(mode,tag,'A =',z,'tail =',info['tail'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=MODES); ap.add_argument('position',choices=['center','plus','minus'])
    ap.add_argument('--comparison-radius',type=int,choices=[20,30,40])
    args=ap.parse_args()
    run(args.mode,args.position,args.comparison_radius)
