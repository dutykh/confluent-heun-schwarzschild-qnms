"""Checks supporting generic_che_resonance_notes.md; no novelty claim."""
import json
from pathlib import Path
import mpmath as mp
mp.mp.dps=100

def reg(alpha,beta,gamma,delta,eta,z,order=400):
    A=delta+alpha*(beta+gamma+2)/2
    B=eta+beta/2+(gamma-alpha)*(beta+1)/2
    pm2,pm1=mp.mpc(0),mp.mpc(1)
    H,Hd=mp.rgamma(beta+1),mp.mpc(0)
    fact,tp=mp.mpf(1),mp.mpc(1)
    for n in range(1,order):
        an=(n-1)*(n+beta+gamma-alpha)+B
        bn=alpha*(n-2)+A
        p=an*pm1+(n-1)*(n-1+beta)*bn*pm2
        fact*=n;tp*=z
        term=p*tp*mp.rgamma(n+beta+1)/fact
        H+=term;Hd+=n*term/z
        pm2,pm1=pm1,p
    return H,Hd

z=mp.mpf('0.3');a=mp.mpf(1);g=mp.mpf('.5');d=-mp.mpf('.75')
V=mp.quad(lambda t:mp.exp(-t)*(1-t)**(-mp.mpf('1.5')),[0,z])
Vd=mp.exp(-z)*(1-z)**(-mp.mpf('1.5'))
rows=[]
for c in [mp.mpf(0),mp.mpf(1),mp.mpc(-2,1)]:
    for epstext in ['1e-10','1e-20','1e-30']:
        eps=mp.mpf(epstext);b=-1+eps;k=c*eps;e=mp.mpf('.5')+k-eps/4
        H,Hd=reg(a,b,g,d,e,z)
        U,Ud=H*mp.gamma(eps),Hd*mp.gamma(eps)
        error=max(abs(U-(1+c*V)),abs(Ud-c*Vd))
        assert error<100*eps
        T,Td=reg(a,-b,g,d,e,z)
        T*=mp.gamma(1-b);Td*=mp.gamma(1-b)
        S=z**(-b)*T;Sd=z**(-b)*(Td-b*T/z)
        F,Fd=(H-k*S)/eps,(Hd-k*Sd)/eps
        W=F*Sd-Fd*S
        exact=-b*mp.rgamma(1+b)/eps*mp.exp(-a*z)*z**(-b-1)*(1-z)**(-g-1)
        assert abs(W-exact)<mp.mpf('1e-65')
        rows.append({'slope':str(c),'epsilon':epstext,'limiting_solution_error':mp.nstr(error,10),'frame_wronskian_error':mp.nstr(abs(W-exact),10)})
# Gamma-regularized function at the simultaneous zero vanishes identically.
H,Hd=reg(a,-1,g,d,mp.mpf('.5'),z)
assert H==0 and Hd==0
result={'status':'PASS','z':str(z),'paths':rows,'base_regularized_solution':'identically zero by recurrence','provenance':'independent checks of classical resonance mechanisms; no novelty claim'}
Path(__file__).with_name('generic_che_resonance_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','path_cases':len(rows),'file':str(Path(__file__).with_name('generic_che_resonance_results.json'))}))
