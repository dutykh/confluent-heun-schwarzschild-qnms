"""Scalar Schwarzschild horizon resonance audit; Python 3 + SymPy + mpmath."""
import sympy as sp
import mpmath as mp
L = sp.symbols('L', nonnegative=True)
def mat(N):
    A=sp.zeros(N)
    for k in range(N):
        A[k,k]=N*N+L-k*(k+1)
        if k<N-1: A[k,k+1]=-(k+1)*(N-k-1)
        if k>0: A[k,k-1]=-N*(N-k)
    return A
for N in range(1,7):
    A=mat(N)
    for k in range(N):
        off=sum(-A[j,k] for j in range(N) if j != k)
        assert sp.simplify(A[k,k]-off-(L+N-k)) == 0
    pprevprev,pprev=sp.Integer(0),sp.Integer(1)
    for n in range(1,N+1):
        p=(n*(n-1)-N*N-L)*pprev-N*(n-1)*(n-1-N)**2*pprevprev
        pprevprev,pprev=pprev,sp.expand(p)
    det=sp.factor(A.det())
    assert sp.simplify(pprev-(-1)**N*det)==0
    if N<=4: print('N=',N,'det A =',det)
print('Column dominance and first regularized coefficient identities: PASS')
mp.mp.dps=90
r=mp.mpf('1.3'); t=1-r

def hh(lam,ell,regularized):
    ll=ell*(ell+1)
    if regularized:
        pm2,pm1=mp.mpf(0),mp.mpf(1)
        val=mp.rgamma(1-lam)
        for n in range(1,200):
            p=(n*(n-1)-lam*lam-ll)*pm1-lam*(n-1)*(n-1-lam)**2*pm2
            val += p*t**n*mp.rgamma(n+1-lam)/mp.factorial(n)
            pm2,pm1=pm1,p
        return val
    hm2,hm1=mp.mpf(0),mp.mpf(1); val=hm1
    for n in range(1,200):
        h=((n*(n-1)-lam*lam-ll)*hm1-lam*(n-1-lam)*hm2)/(n*(n-lam))
        val += h*t**n
        hm2,hm1=hm1,h
    return val

def uh(lam,ell,regularized):
    # lambda=2 i Omega, so e^(-2iOmega)=e^-lambda and (r-1)^(-iOmega)=(r-1)^(-lambda/2)
    return mp.exp(-lam)*mp.exp(lam*r/2)*r*(r-1)**(-lam/2)*hh(lam,ell,regularized)
for N in range(1,7):
    for ell in range(3):
        det=int(mat(N).det().subs(L,ell*(ell+1)))
        lhs=uh(mp.mpf(N),ell,True)
        rhs=mp.exp(-N)*det/mp.factorial(N)*uh(mp.mpf(-N),ell,False)
        err=abs(lhs-rhs)/max(abs(rhs),mp.mpf(1))
        assert err<mp.mpf('1e-65'), (N,ell,err)
print('Upper-half-plane identity: PASS for N=1..6, ell=0..2 at r=1.3 (90 digits, 200 terms)')
for eps in [mp.mpf('1e-10'),mp.mpf('1e-20'),mp.mpf('1e-30')]:
    lam=1+eps
    assert abs(hh(lam,0,True)-hh(lam,0,False)*mp.rgamma(1-lam))<mp.mpf('1e-65')
print('Regularization equality off resonance: PASS near lambda=1')
