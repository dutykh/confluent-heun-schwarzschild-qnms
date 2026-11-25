"""Real-frequency flux check of the full amplitude normalization."""
import cmath
import json
import math
from pathlib import Path
import qnm_contour as c

def run():
    out=[]
    for ell in (0,1,2):
        for w in (.5,1.):
            A,d=c.incoming_amplitude(complex(w),ell,radius=100,order=48,
                                      rtol=3e-14,diagnostics=True)
            m=complex(*d['outgoing_scaled_value'])
            mx=complex(*d['outgoing_scaled_derivative'])
            uH,uHx=c.horizon_state(w,ell,rtol=3e-14)
            x=1.5+math.log(.5)
            phase=cmath.exp(1j*w*x)
            u=phase*m
            ux=phase*(mx+1j*w*m)
            Ao=(u.conjugate()*uHx-ux.conjugate()*uH)/(2j*w)
            err=abs(abs(A)**2-abs(Ao)**2-1)
            out.append(dict(ell=ell,omega=w,A_in=[A.real,A.imag],
                            A_out=[Ao.real,Ao.imag],flux_identity_error=err))
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':run()
