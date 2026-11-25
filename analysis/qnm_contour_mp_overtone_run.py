"""Additional high precision l=0,n=2 diagnostic; modifies only MP results JSON."""
import json,time
from pathlib import Path
import mpmath as mp
from qnm_contour_mp import incoming_amplitude,solve_root,complex_record,serialize
mp.mp.dps=40
path=Path(__file__).with_name('qnm_contour_mp_results.json')
data=json.loads(path.read_text())
cf=json.loads(Path(__file__).with_name('qnm_leaver_results.json').read_text())
ref=next(x for x in cf['modes'] if x['ell']==0 and x['overtone']==2)
wr=mp.mpc(ref['omega']['real'],ref['omega']['imag'])
mode={'ell':0,'n':2,'reference_cf':complex_record(wr),'reference_diagnostics':[]}; start=time.time()
configs=[dict(anchor='5',radius=100,asympt_order=50,taylor_order=32,step_fraction='.25',scaled_horizon=True),
         dict(anchor='5',radius=140,asympt_order=65,taylor_order=40,step_fraction='.25',scaled_horizon=True),
         dict(anchor='5',radius=140,asympt_order=65,taylor_order=48,step_fraction='.20',scaled_horizon=True)]
for cfg in configs:
    A,diag=incoming_amplitude(wr,0,diagnostics=True,**cfg)
    mode['reference_diagnostics'].append({'settings':cfg,'A':A,'abs_A':abs(A),**diag})
    print('overtone diagnostic',cfg,'|A|',mp.nstr(abs(A),12),'tail',mp.nstr(diag['tail_bounds']['amplitude_tail_only'],12),flush=True)
    data['modes']=[x for x in data['modes'] if not(x['ell']==0 and x['n']==2)]+[mode]
    path.write_text(json.dumps(serialize(data),indent=2)+'\n')
settings=configs[-1]
root,calls=solve_root(0,('.15','-1.20'),settings,40)
ar,diag=incoming_amplitude(root,0,diagnostics=True,**settings)
mode.update({'root':root,'root_seed':['.15','-1.20'],'root_calls':calls,'root_settings':settings,
             'root_amplitude':ar,'root_amplitude_abs':abs(ar),'difference_from_cf':abs(root-wr),
             'root_diagnostics':diag})
with mp.workdps(55):
    ah,dh=incoming_amplitude(root,0,diagnostics=True,**settings)
mode['precision_check']={'dps':55,'A':ah,'amplitude_change':abs(ah-ar)}
# Change both Taylor accuracy and anchor; this is a fixed-root diagnostic.
with mp.workdps(55):
    changed={**settings,'anchor':'5.4','taylor_order':56,'step_fraction':'.18'}
    ac,dc=incoming_amplitude(root,0,diagnostics=True,**changed)
mode['changed_anchor_order_check']={'settings':changed,'dps':55,'A':ac,'abs_A':abs(ac),**dc}
mode['elapsed_seconds']=round(time.time()-start,3)
path.write_text(json.dumps(serialize(data),indent=2)+'\n')
print('overtone root',mp.nstr(root,38),'difference CF',mp.nstr(abs(root-wr),12),'changed anchor |A|',mp.nstr(abs(ac),12),flush=True)
