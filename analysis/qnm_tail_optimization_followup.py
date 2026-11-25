"""Recompute guarded weighted bounds, and add legacy-adaptive controlled timings."""
import json,time
from pathlib import Path
import mpmath as mp
import qnm_tail_optimization as opt
import qnm_contour_mp as contour
mp.mp.dps=40
path=Path(__file__).with_name('qnm_tail_optimization_results.json');data=json.loads(path.read_text())
cf=json.loads(Path(__file__).with_name('qnm_leaver_results.json').read_text())
for rec in data['modes']:
    ell,n,c=rec['ell'],rec['n'],rec['anchor']
    center=mp.mpc(rec['frequency_disk_center']['real'],rec['frequency_disk_center']['imag']);dom=opt.Domain(center,rec['frequency_disk_radius'])
    ref=next(v for v in cf['modes'] if v['ell']==ell and v['overtone']==n)
    wr=mp.mpc(ref['omega']['real'],ref['omega']['imag']);point=opt.Domain(wr)
    # Selection itself used no expint; only refined reported bounds need refreshing.
    for key,level in [('selection_legacy','legacy'),('selection_improved','improved')]:
        ss=rec[key];ss['uniform_bound']=opt.tail_bound(dom,ell,c,ss['radius'],ss['order'],level=level,refine=True)
    s=rec['selection_improved'];R,N=s['radius'],s['order']
    rec['point_bound_comparison']={lvl:opt.tail_bound(point,ell,c,R,N,level=lvl,refine=True) for lvl in ('legacy','geometry','improved')}
    for tag,ev in rec['evaluations'].items():
        cfg=ev['settings'];bd=opt.tail_bound(point,ell,c,cfg['radius'],cfg['asympt_order'],level='improved',refine=True)
        ev['improved_tail_bound']=bd;ev['amplitude_tail_bound_diagnostic']=opt.amplitude_bound(bd,wr,ell,c)
    # Quantify selector setup time separately; it can be amortized over root calls.
    rec['selector_timing_seconds']={}
    for level in ('legacy','improved'):
        st=time.perf_counter();selected=opt.optimize(dom,ell,c,mp.mpf(data['target_uniform_outgoing_state_tail_norm']),level)
        rec['selector_timing_seconds'][level]=time.perf_counter()-st
        assert (selected['radius'],selected['order'])==(rec['selection_'+level]['radius'],rec['selection_'+level]['order'])
    legacy=rec['selection_legacy'];cfg={**rec['evaluations']['adaptive_fixed_domain']['settings'],'radius':legacy['radius'],'asympt_order':legacy['order']}
    timings=[]
    for _ in range(2):
        st=time.perf_counter();A,diag=contour.incoming_amplitude(wr,ell,diagnostics=True,**cfg);timings.append(time.perf_counter()-st)
    bd=opt.tail_bound(point,ell,c,cfg['radius'],cfg['asympt_order'],level='improved',refine=True)
    rec['evaluations']['legacy_adaptive']={'settings':cfg,'A_at_CF':A,'abs_A_at_CF':abs(A),'wall_seconds':timings,'outgoing_steps':diag['steps'],'horizon_steps':diag['horizon_steps'],'improved_tail_bound':bd,'amplitude_tail_bound_diagnostic':opt.amplitude_bound(bd,wr,ell,c)}
    print('legacy adaptive eval',ell,n,'R,N',legacy['radius'],legacy['order'],'seconds',timings,'steps',diag['steps'],'|A|',mp.nstr(abs(A),8),flush=True)
    # Independent same-seed, same-tolerance root search. Enforce the selected disk.
    calls=0;distmax=mp.mpf(0);trace=[]
    def residual(w):
        global calls,distmax
        calls+=1;distmax=max(distmax,abs(w-dom.center))
        if abs(w-dom.center)>dom.radius: raise ValueError('Secant iterate left selected disk')
        value=contour.incoming_amplitude(w,ell,**cfg)
        trace.append({'omega':w,'abs_A':abs(value)})
        return value
    st=time.perf_counter()
    with mp.workdps(40):
        root=mp.findroot(residual,(center,center+mp.mpc('.001','.001')),solver='secant',tol=mp.mpf('1e-29'),maxsteps=20,verify=True)
    elapsed=time.perf_counter()-st
    A,diag=contour.incoming_amplitude(root,ell,diagnostics=True,**cfg)
    rec['root_searches']['legacy_adaptive']={'root':root,'calls':calls,'wall_seconds':elapsed,'abs_A':abs(A),'distance_from_disk_center':abs(root-dom.center),'maximum_iterate_distance_from_disk_center':distmax,'difference_from_cf':abs(root-wr),'frequency_domain_enforced':True,'trace':trace}
    print('legacy adaptive root',ell,n,'seconds',elapsed,'CFdiff',mp.nstr(abs(root-wr),8),flush=True)
    for tag,rr in rec['root_searches'].items():
        rt=rr['root'];omega=mp.mpc(rt['real'],rt['imag']) if isinstance(rt,dict) else rt
        ecfg=rec['evaluations'][tag]['settings']
        sb=opt.tail_bound(opt.Domain(omega),ell,c,ecfg['radius'],ecfg['asympt_order'],level='improved',refine=True)
        rr['amplitude_tail_bound_at_root_diagnostic']=opt.amplitude_bound(sb,omega,ell,c)
    rec['additional_bound_evaluation_speedup']=sum(timings)/sum(rec['evaluations']['adaptive_fixed_domain']['wall_seconds'])
    rec['additional_bound_root_search_speedup']=elapsed/rec['root_searches']['adaptive_fixed_domain']['wall_seconds']
    data['weighted_integral_evaluation']='Guarded extra precision for exp(a)E_p(a), positivity checked; not outward rounded.'
    opt.dump(path,data)
# Independent positive quadrature spot checks for the most cancellation-prone expint evaluation.
checks=[]
for c,R,p,mu in [('1.5',140,67,'1.92'),('1.5',35,64,'1.90'),('5',30,41,'.20')]:
    c=mp.mpf(c);R=mp.mpf(R);mu=mp.mpf(mu);r=mp.sqrt(c*c+R*R)
    with mp.workdps(70):
        integral=r**(-p)*mp.quad(lambda u:mp.exp(-mu*u)*(r*r/(c*c+(R+u)**2))**(mp.mpf(p)/2),[0,1,10,mp.inf])
        bound=opt.weighted_power_upper(c,R,p,mu)
    assert bound>=integral
    checks.append({'c':c,'R':R,'p':p,'mu':mu,'upper_to_quadrature_ratio':bound/integral})
data['weighted_integral_spot_checks']=checks;opt.dump(path,data)
print('Weighted-integral quadrature checks passed',flush=True)
