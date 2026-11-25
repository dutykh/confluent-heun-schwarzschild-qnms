#!/usr/bin/env python3
"""Generic CHE path-dependent resonant limit and actual finite-point matrices.

Maple convention H''+[a+(b+1)/z+(g+1)/(z-1)]H'+(A*z+B)/(z*(z-1))*H=0.
The normalization H(0)=1 gives H'(0)=B/(b+1). We approach b=-1,B=0
along B=c*(b+1) and verify the full basis and connection-matrix limits.
The regularized function Phi=H/Gamma(b+1) vanishes at this intersection.
"""
import json
from pathlib import Path
import mpmath as mp
import sympy as sp


def ab(p):
    a,b,g,d,e=p
    return d+a*(b+g+2)/2,e+b/2+(g-a)*(b+1)/2


def local_coefficients(p,N=800):
    a,b,g,_,_=p;A,B=ab(p)
    aa=[mp.mpc(1)];previous=mp.mpc(0)
    for n in range(N):
        current=aa[-1]
        aa.append(((n*(n+b+g+1-a)+B)*current+(a*(n-1)+A)*previous)/((n+1)*(n+b+1)))
        previous=current
    return aa


def evaluate_series(aa,z):
    val=aa[-1];der=mp.mpc(0)
    for a in reversed(aa[:-1]):der,val=der*z+val,val*z+a
    return val,der


def origin_basis(p,z,N=800):
    a,b,g,d,e=p
    f=evaluate_series(local_coefficients(p,N),z)
    h=evaluate_series(local_coefficients((a,-b,g,d,e),N),z)
    q=(z**(-b)*h[0],z**(-b)*(h[1]-b*h[0]/z))
    return mp.matrix([[f[0],q[0]],[f[1],q[1]]])


def one_basis(p,z,N=800):
    a,b,g,d,e=p
    Y=origin_basis((-a,g,b,-d,e+d),1-z,N)
    Y[1,0]=-Y[1,0];Y[1,1]=-Y[1,1]
    return Y


def regularized_polynomial_series(p,z,N=350):
    a,b,g,_,_=p;A,B=ab(p)
    previous=mp.mpc(0);current=mp.mpc(1)
    total=mp.rgamma(b+1);factorial=mp.mpf(1);power=mp.mpf(1)
    for n in range(N):
        nxt=(n*(n+b+g+1-a)+B)*current+n*(n+b)*(a*(n-1)+A)*previous
        factorial*=n+1;power*=z
        total+=nxt*mp.rgamma(b+n+2)*power/factorial
        previous,current=current,nxt
    return total


def ordinary_basis(a,g,A,z,N=800):
    columns=[]
    for initial in [(1,0),(0,1)]:
        aa=list(map(mp.mpf,initial))
        for n in range(N-1):
            aa.append(((n+1)*(n+g+1-a)*aa[-1]+(a*n+A)*aa[-2])/((n+2)*(n+1)))
        columns.append(evaluate_series(aa,z))
    return mp.matrix([[columns[0][0],columns[1][0]],[columns[0][1],columns[1][1]]])


def symbolic_checks():
    n,a,b,g,A,B=sp.symbols('n alpha beta gamma A B')
    am,an,ap=sp.symbols('am an ap')
    recurrence=-(n+1)*(n+b+1)*ap+(n*(n+b+g+1-a)+B)*an+(a*(n-1)+A)*am
    collected=(-n*(n+1)-(b+1)*(n+1))*ap+(n*(n-1)+(-a+b+g+2)*n+B)*an+(a*(n-1)+A)*am
    assert sp.expand(recurrence-collected)==0
    eps,c=sp.symbols('epsilon c')
    p=(sp.Rational(7,10),-1+eps,sp.Rational(2,5),sp.Rational(1,5),sp.Rational(1,2)+(c-sp.Rational(7,20))*eps)
    AA,BB=ab(p)
    assert sp.simplify(AA-(sp.Rational(69,100)+sp.Rational(7,20)*eps))==0
    assert sp.simplify(BB-c*eps)==0
    assert sp.simplify(BB/(p[1]+1)-c)==0
    return {'coefficient_recurrence_difference':'0','A_path':str(sp.expand(AA)),
            'B_path':str(sp.expand(BB)),'normalized_origin_derivative':str(sp.simplify(BB/(p[1]+1)))}


def strings(value,digits=50):
    if isinstance(value,mp.matrix):return [[strings(value[i,j],digits) for j in range(value.cols)] for i in range(value.rows)]
    return {'real':mp.nstr(mp.re(value),digits),'imag':mp.nstr(mp.im(value),digits)}


def maxabs(M):return max(abs(x) for x in M)


def finite_point_checks(a,g,d):
    eps=mp.mpf('1e-4');c=mp.mpf(1)
    p=(a,-1+eps,g,d,mp.mpf(1)/2+(c-mp.mpf(7)/20)*eps)
    points=[mp.mpf(3)/10,mp.mpf(1)/2,mp.mpf(7)/10]
    matrices=[one_basis(p,z)**-1*origin_basis(p,z) for z in points]
    discrepancy=max(maxabs(C-matrices[0]) for C in matrices)
    determinant_error=abs(mp.det(matrices[0])+p[1]*mp.exp(-a)/g)
    assert discrepancy<mp.mpf('1e-85') and determinant_error<mp.mpf('1e-85')
    return {'epsilon':str(eps),'c':str(c),'matching_points':list(map(str,points)),
            'matrix_matching_point_max_discrepancy':mp.nstr(discrepancy,20),
            'determinant_identity':'det(C01)=-beta*exp(-alpha)/gamma',
            'determinant_identity_error':mp.nstr(determinant_error,20)}


def main():
    mp.mp.dps=100
    a=mp.mpf(7)/10;g=mp.mpf(2)/5;d=mp.mpf(1)/5;z=mp.mpf(3)/10;A0=mp.mpf(69)/100
    p0=(a,mp.mpf(-1),g,d,mp.mpf(1)/2)
    U=ordinary_basis(a,g,A0,z)
    expected_w=mp.exp(-a*z)*(1-z)**(-g-1)
    assert abs(mp.det(U)-expected_w)<mp.mpf('1e-90')
    def vf(t,Y):return [Y[1],-(a+(g+1)/(t-1))*Y[1]-A0/(t-1)*Y[0]]
    cc=[mp.odefun(vf,0,initial,tol=mp.mpf('1e-80'),degree=40)(z) for initial in [(1,0),(0,1)]]
    Uode=mp.matrix([[cc[0][0],cc[1][0]],[cc[0][1],cc[1][1]]])
    ode_error=maxabs(U-Uode);assert ode_error<mp.mpf('1e-75')
    Y1lim=one_basis(p0,z);Cref=Y1lim**-1*U
    assert maxabs(Y1lim*Cref-U)<mp.mpf('1e-90')
    Phi0=regularized_polynomial_series(p0,z);assert Phi0==0
    out={'convention':'Maple alpha beta gamma delta eta; Y0=Y1*C01; derivatives in z; real powers on (0,1)',
         'working_decimal_digits':100,'z':str(z),'alpha':str(a),'gamma':str(g),'delta':str(d),
         'symbolic_checks':symbolic_checks(),'finite_point_checks':finite_point_checks(a,g,d),
         'limiting_ordinary_equation':'u\'\'+[0.7+1.4/(z-1)]u\'+0.69/(z-1)*u=0',
         'ordinary_reference_basis':strings(U),'ordinary_reference_Wronskian':strings(mp.det(U)),
         'expected_Wronskian':strings(expected_w),'independent_ode_max_difference':mp.nstr(ode_error,20),
         'limiting_normalized_basis_at_one':strings(Y1lim),'Cref_from_ordinary_initial_data_to_one':strings(Cref),
         'regularized_solution_at_degenerate_point':strings(Phi0),'paths':[]}
    for c in map(mp.mpf,[0,1,-2]):
        T=mp.matrix([[1,0],[c,1]]);target=U*T;Ctarget=Cref*T
        path={'c':str(c),'limit_basis':strings(target),'limit_connection_matrix':strings(Ctarget),
              'basis_change_Tc':strings(T),'runs':[]}
        for eps in [mp.mpf(10)**(-k) for k in [2,4,8,12,20]]:
            p=(a,-1+eps,g,d,mp.mpf(1)/2+(c-mp.mpf(7)/20)*eps)
            Y0=origin_basis(p,z);Y1=one_basis(p,z);C=Y1**-1*Y0
            Phi=regularized_polynomial_series(p,z)
            phi_identity=abs(Phi-mp.rgamma(eps)*Y0[0,0]);assert phi_identity<mp.mpf('1e-85')
            W0=(1-eps)*mp.exp(-a*z)*z**(-eps)*(1-z)**(-g-1)
            assert abs(mp.det(Y0)-W0)<mp.mpf('1e-80')
            path['runs'].append({'epsilon':str(eps),'eta':str(p[-1]),'A':strings(ab(p)[0]),'B':strings(ab(p)[1]),
                'H0_at_z':strings(Y0[0,0]),'H0prime_at_z':strings(Y0[1,0]),
                'H0_origin_derivative':strings(ab(p)[1]/(p[1]+1)),
                'basis_max_difference_from_limit':mp.nstr(maxabs(Y0-target),20),
                'connection_max_difference_from_limit':mp.nstr(maxabs(C-Ctarget),20),
                'connection_matrix_C01':strings(C),'U_inverse_Y0':strings(U**-1*Y0),
                'regularized_value':strings(Phi),'regularized_over_epsilon':strings(Phi/eps),
                'gamma_regularization_identity_error':mp.nstr(phi_identity,20)})
        assert mp.mpf(path['runs'][-1]['basis_max_difference_from_limit'])<mp.mpf('1e-18')
        assert mp.mpf(path['runs'][-1]['connection_max_difference_from_limit'])<mp.mpf('1e-18')
        out['paths'].append(path)
    out['conclusion']='H0 tends to u0+c*u1; the full origin Frobenius basis tends to U*[[1,0],[c,1]], and C01 tends to Cref*[[1,0],[c,1]]. Unit origin value does not produce a path-independent resonant continuation.'
    dest=Path(__file__).with_name('generic_che_path_verify.json')
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print('ordinary U=',mp.nstr(U,30));print('Cref=',mp.nstr(Cref,30));print(dest.resolve())

if __name__=='__main__':main()
