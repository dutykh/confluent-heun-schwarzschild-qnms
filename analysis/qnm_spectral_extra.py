import json,time
from pathlib import Path
from qnm_spectral_refine import refine
src=json.load(open('analysis/qnm_spectral_highprecision.json'))
rows=[r for r in src['results'] if r['ell']==0 and r['N']==320]
seeds={r['provisional_overtone']:complex(*map(float,r['omega_strings'])) for r in rows}
out={'method':'additional high-resolution ell=0 continuation','results':[]}
for N in [400,480]:
 for n in [1,2]:
  start=time.time()
  r=refine(0,N,seeds[n],50)
  r['provisional_overtone']=n
  r['elapsed_seconds']=time.time()-start
  seeds[n]=complex(*map(float,r['omega_strings']))
  out['results'].append(r)
  Path('analysis/qnm_spectral_highprecision_extra.json').write_text(json.dumps(out,indent=2)+'\n')
  print(N,n,r['omega_strings'],'seconds',r['elapsed_seconds'],flush=True)
