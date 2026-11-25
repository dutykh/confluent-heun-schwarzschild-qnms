"""Focused validation checks for the rigorous finite transport implementation."""
import json,time
from pathlib import Path
from mpmath import mp,iv
from cert.transport import transport_horizon,horizon_enclosure,as_complex,up_abs,lower,upper


def overlap(a,b):
    a,b=as_complex(a),as_complex(b)
    return (lower(a.real)<=upper(b.real) and lower(b.real)<=upper(a.real)
            and lower(a.imag)<=upper(b.imag) and lower(b.imag)<=upper(a.imag))


def contains(a,b):
    a,b=as_complex(a),as_complex(b)
    return lower(a.real)<=lower(b.real)<=upper(b.real)<=upper(a.real) and lower(a.imag)<=lower(b.imag)<=upper(b.imag)<=upper(a.imag)


def main():
    mp.dps=100;iv.dps=100
    records=[]
    # Exact static scalar solutions h=z P_l(2z-1), l=0,1,2.
    for ell in (0,1,2):
        ans=transport_horizon(iv.mpc(0,0),ell,'1.5',3,dps=100,order=72,step_fraction='.18')
        z=iv.mpc(iv.mpf('1.5'),iv.mpf(3));f=(z-1)/z
        if ell==0:h=z;hx=f
        elif ell==1:h=z*(2*z-1);hx=f*(4*z-1)
        else:h=z*(6*z*z-6*z+1);hx=f*(18*z*z-12*z+1)
        assert contains(ans['h'],h) and contains(ans['hx'],hx),(ell,ans,h,hx)
        records.append({'test':'exact static polynomial','ell':ell,'steps':ans['steps'],
                        'seconds':ans['elapsed_seconds'],'radii':[mp.nstr(x,20) for x in ans['radii']]})
    cases=[(0,'0.22090987816083937175092301233587754','-0.20979143417376191756347813371075195','1.5',80),
           (2,'0.96728774442142597345555887030351716','-0.19351755195657572531730746670174375','1.5',20),
           (0,'0.15148387107035152594044328824134263','-1.2021571800716070281873120374150582','5',30)]
    for ell,re,im,a,R in cases:
        om=iv.mpc(iv.mpf(re),iv.mpf(im))
        lo=transport_horizon(om,ell,a,R,dps=100,order=96,step_fraction='.18')
        hi=transport_horizon(om,ell,a,R,dps=120,order=112,step_fraction='.16',horizon_terms=260)
        assert overlap(lo['h'],hi['h']) and overlap(lo['hx'],hi['hx'])
        assert max(lo['radii'])<mp.mpf('1e-43'),lo
        # The refined rectangle may extend beyond a low-order disk rectangle by
        # rounding at their edges; intersection is the generic rigorous test.
        records.append({'test':'QNM parameter point and refinement','ell':ell,'omega':[re,im],'anchor':a,'R':R,
                        'steps':lo['steps'],'seconds':lo['elapsed_seconds'],
                        'radii':[mp.nstr(x,25) for x in lo['radii']],
                        'refined_steps':hi['steps'],'refined_seconds':hi['elapsed_seconds'],
                        'refined_radii':[mp.nstr(x,25) for x in hi['radii']],
                        'refined_box_contained':contains(lo['h'],hi['h']) and contains(lo['hx'],hi['hx'])})
    # Nonpoint parameter rectangle: proof must remain valid, not treat midpoint as exact.
    om=iv.mpc(iv.mpf(['.21','.23']),iv.mpf(['-.22','-.20']))
    hz=horizon_enclosure(om,0,dps=100)
    assert hz['majorant_contraction']<1
    for re,im in [('.21','-.22'),('.23','-.20'),('.22','-.21')]:
        p=horizon_enclosure(iv.mpc(iv.mpf(re),iv.mpf(im)),0,dps=100)
        assert contains(hz['h'],p['h']) and contains(hz['hx'],p['hx'])
    records.append({'test':'horizon parameter rectangle','contraction':mp.nstr(hz['majorant_contraction'],20),
                    'radii':[mp.nstr(x,20) for x in hz['radii']]})
    Path('cert/transport_test_results.json').write_text(json.dumps({'passed':True,'records':records},indent=2)+'\n')
    print(json.dumps({'passed':True,'records':records},indent=2))

if __name__=='__main__':main()
