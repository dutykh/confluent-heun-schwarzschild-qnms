"""Identical endpoint bound comparison and forward/backward numerical controls."""
import json,time
from pathlib import Path
import mpmath as mp
import qnm_forward_coefficient as fwd
import qnm_contour_mp as contour
import qnm_tail_optimization as bounds
mp.mp.dps=40
path=Path(__file__).with_name('qnm_forward_coefficient_results.json');data=json.loads(path.read_text())
def zz(v):return mp.mpc(v['real'],v['imag'])
for rec in data['modes']:
    ell,n,c=rec['ell'],rec['n'],rec['anchor'];wr=zz(rec['CF_frequency']);seed=rec['seed']
    center=mp.mpc(*seed);dom=bounds.Domain(center,'.01')
    for diag in rec['radius_diagnostics']:
        cfg=diag['settings'];R,N=cfg['radius'],cfg['asympt_order']
        h,hx=zz(diag['h_R']),zz(diag['hx_R'])
        diag['identical_endpoint_bounds_at_CF']=fwd.coefficient_tail_comparison(wr,ell,c,R,N,h,hx)
        A,dd=fwd.forward_coefficient(center,ell,diagnostics=True,**cfg)
        diag['offroot_seed_diagnostic']={'A_N':A,'bounds':fwd.coefficient_tail_comparison(center,ell,c,R,N,dd['h_R'],dd['hx_R'])}
        print('bound ratios',ell,n,'R',R,'nearroot state/direct',mp.nstr(diag['identical_endpoint_bounds_at_CF']['state_over_direct_ratio'],6),'state/weighted',mp.nstr(diag['identical_endpoint_bounds_at_CF']['state_over_weighted_ratio'],6),'offroot state/direct',mp.nstr(diag['offroot_seed_diagnostic']['bounds']['state_over_direct_ratio'],6),flush=True)
    for run in rec['root_searches']:
        w=zz(run['root']);cfg=run['settings']
        run['coefficient_tail_bound_comparison']=fwd.coefficient_tail_comparison(w,ell,c,cfg['radius'],cfg['asympt_order'],zz(run['h_R']),zz(run['hx_R']))
    chosen=rec['root_searches'][rec['selected_root_run_index']];cfg=chosen['settings']
    bc={**cfg,'scaled_horizon':True}
    st=time.perf_counter();Ab,db=contour.incoming_amplitude(wr,ell,diagnostics=True,**bc);etime=time.perf_counter()-st
    Af,df=fwd.forward_coefficient(wr,ell,diagnostics=True,**cfg)
    st=time.perf_counter();wb,calls=contour.solve_root(ell,seed,bc,40);rtime=time.perf_counter()-st
    rec['identical_endpoint_backward_control']={'settings':bc,'A_at_CF':Ab,
       'difference_from_forward_A_at_CF':abs(Ab-Af),'evaluation_wall_seconds':etime,
       'outgoing_backward_steps':db['steps'],'root':wb,'root_calls':calls,'root_wall_seconds':rtime,
       'difference_from_CF':abs(wb-wr),'difference_from_forward_root':abs(wb-zz(chosen['root']))}
    print('same endpoint backward root',ell,n,'seconds',rtime,'CFdiff',mp.nstr(abs(wb-wr),8),'fwd rootdiff',mp.nstr(abs(wb-zz(chosen['root'])),8),flush=True)
    bounds.dump(path,data)
