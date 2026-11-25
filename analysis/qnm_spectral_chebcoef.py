"""Chebyshev coefficient form of same collocation, to avoid nodal D^2 conditioning."""
import numpy as np
from scipy.linalg import eig
from numpy.polynomial.chebyshev import chebder,chebmul

def derivative_integer(v):
    out=np.zeros(max(1,len(v)-1),dtype=np.int64)
    for k in range(1,len(v)):
        for j in range(k-1,-1,-2):
            out[j]+=v[k]*k*(1 if j==0 else 2)
    return out

def product_scaled32(a,b_scaled16):
    out=np.zeros(len(a)+len(b_scaled16)-1,dtype=np.int64)
    for i,ai in enumerate(a):
        for j,bj in enumerate(b_scaled16):
            out[i+j]+=ai*bj
            out[abs(i-j)]+=ai*bj
    return out

def coefficient_matrices(ell,n):
    # Integer assembly: each returned entry is exactly an integer /32.
    # Tk Tl=(T(k+l)+T|k-l|)/2, and derivatives of Tk have integer coefficients.
    out=[np.zeros((n+1,n+1)) for _ in range(3)]
    for j in range(n+1):
        v=np.zeros(j+1,dtype=np.int64); v[-1]=1
        v1=derivative_integer(v); v2=derivative_integer(v1)
        terms=[[(v2,[4,-2,-4,2]),(v1,[4,-16,12]),(v,[-16*ell*(ell+1)-8,8])],
               [(v1,[16,64,-16]),(v,[32,-32])],[(v,[96,-32])]]
        for index,tt in enumerate(terms):
            col=np.zeros(n+2,dtype=np.int64)
            for a,b in tt:
                pp=product_scaled32(a,b)
                col[:len(pp)]+=pp
            col[n-1]+=col[n+1]
            out[index][:,j]=col[:n+1]/32
    return out

def solve(ell,n):
    A0,A1,A2=coefficient_matrices(ell,n)
    # row scale (same polynomial roots)
    scale=np.maximum(np.max(abs(A0),axis=1),np.maximum(np.max(abs(A1),axis=1),np.max(abs(A2),axis=1)))
    A0=A0/scale[:,None]; A1=A1/scale[:,None]; A2=A2/scale[:,None]
    I=np.eye(n+1); Z=np.zeros_like(I)
    M=np.block([[Z,I],[A0,A1]])
    B=np.block([[I,Z],[Z,A2]])
    vals=eig(M,B,right=False,check_finite=False)
    vals=-1j*vals
    vals=[v for v in vals if 1e-5<v.real<2 and -2<v.imag<-1e-5]
    return sorted(vals,key=lambda z:-z.imag)
if __name__=='__main__':
 for ell in range(3):
  for n in [40,60,80,100,120]:
   print(ell,n,solve(ell,n)[:6])
