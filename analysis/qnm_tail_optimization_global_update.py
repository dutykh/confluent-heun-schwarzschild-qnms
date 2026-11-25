"""Append single-global-residual estimates without changing timed computations."""
import json
from pathlib import Path
import mpmath as mp
import qnm_tail_optimization as opt
mp.mp.dps=40
path=Path(__file__).with_name('qnm_tail_optimization_results.json');data=json.loads(path.read_text())
for rec in data['modes']:
    ell,n,c=rec['ell'],rec['n'],rec['anchor']
    center=mp.mpc(rec['frequency_disk_center']['real'],rec['frequency_disk_center']['imag'])
    dom=opt.Domain(center,rec['frequency_disk_radius'])
    chosen=rec['selection_improved'];R,N=chosen['radius'],chosen['order']
    rec['uniform_global_residual_bound_at_selected_settings']=opt.global_residual_bounds(dom,ell,c,R,N)
    # Record the complete searched order profile at the chosen radius (unrefined
    # selection bound), rather than conflating smallest feasible with best N.
    nmax=max(80,int(mp.ceil(4*dom.omin*mp.sqrt(mp.mpf(c)**2+R*R)))+20)
    aa=opt.coefficient_data(dom,ell,nmax)
    rec['order_profile_at_selected_radius']=[{'N':nn,'uniform_state_bound':opt.tail_bound(dom,ell,c,R,nn,aa,level='improved',refine=False)['state_norm']} for nn in range(6,nmax+1)]
    for tag,run in rec['root_searches'].items():
        w=mp.mpc(run['root']['real'],run['root']['imag']);cfg=rec['evaluations'][tag]['settings']
        run['global_residual_diagnostics']=opt.direct_amplitude_bounds(w,ell,c,cfg['radius'],cfg['asympt_order'])
        print('global',ell,n,tag,'ordinary_amp',run['amplitude_tail_bound_at_root_diagnostic'],'global_direct',mp.nstr(run['global_residual_diagnostics']['direct_kernel_identity_amplitude'],8),'global_state',mp.nstr(run['global_residual_diagnostics']['global_state_amplitude'],8),flush=True)
    # Check whether extra refinement could make the preceding sampled radius
    # feasible at its best unrefined order; this is diagnostic only.
    prior=chosen['search'][-2] if len(chosen['search'])>1 else None
    if prior:
        pr,pn=prior['radius'],prior['best_order']
        rec['preceding_radius_refinement_check']={'R':pr,'N':pn,
            'two_stage_refined_state':opt.tail_bound(dom,ell,c,pr,pn,level='improved',refine=True)['state_norm'],
            'global_state':opt.global_residual_bounds(dom,ell,c,pr,pn)['state_norm']}
opt.dump(path,data)
