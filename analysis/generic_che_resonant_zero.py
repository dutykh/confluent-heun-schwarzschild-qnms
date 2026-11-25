#!/usr/bin/env python3
"""Resolve one generic CHE connection zero through an apparent-origin resonance.

At epsilon=0 the equation reduces exactly to Kummer. Its Tricomi solution gives
an independent normalization/phase check of the generic infinity continuation.
For epsilon!=0 the generic CHE solver is used, with B=epsilon*t; the limiting
basis slope is chosen by zero exponential connection amplitude.
No novelty claim is attached to this illustrative resolved parameter curve.
"""
import json,time
from pathlib import Path
import mpmath as mp
from generic_che_infinity_numeric import infinity_state,algebraic_state_from_m,coefficient_from_states,coefficient,record,asymptotic_coefficients,polynomial_state
from generic_che_path_verify import ordinary_basis


def parameters(epsilon,t):
    return (mp.mpf(7)/10,-1+epsilon,mp.mpf(2)/5,mp.mpf(1)/5,
            mp.mpf(1)/2+epsilon*(t-mp.mpf(7)/20))


def main():
    mp.mp.dps=75
    alpha=mp.mpf(7)/10;aa=mp.mpf(69)/70;bb=mp.mpf(7)/5
    z=mp.mpc('.5','.5');p0=parameters(mp.mpf(0),mp.mpf(0))
    pref=mp.exp(aa*mp.log(alpha)-mp.j*mp.pi*aa)
    U0=mp.hyperu(aa,bb,alpha)
    U1=mp.hyperu(aa+1,bb+1,alpha)
    cstar=alpha*aa*U1/U0
    exact_y=pref*mp.hyperu(aa,bb,alpha*(1-z))
    exact_yp=pref*alpha*aa*mp.hyperu(aa+1,bb+1,alpha*(1-z))
    B0_exact=mp.exp(mp.j*mp.pi*bb)*pref*aa*U1
    B1_exact=-mp.exp(mp.j*mp.pi*bb)*pref*U0/alpha
    # Positive real Tricomi integral independently verifies the slope formula.
    U0int=mp.quad(lambda t:mp.exp(-alpha*t)*t**(aa-1)*(1+t)**(bb-aa-1),[0,1,mp.inf])/mp.gamma(aa)
    U1int=mp.quad(lambda t:mp.exp(-alpha*t)*t**aa*(1+t)**(bb-aa-1),[0,1,mp.inf])/mp.gamma(aa+1)
    cref_integral=alpha*aa*U1int/U0int
    ordinary=ordinary_basis(alpha,bb-1,mp.mpf(69)/100,z)
    out={'working_decimal_digits':75,'normalization':'Omega-independent generic CHE, algebraic branch on ray c=-1 with starting log z=log r+i*pi',
         'limiting_parameters':p0,'Kummer_a':aa,'Kummer_b':bb,
         'exact_yalg':exact_y,'exact_yalg_prime':exact_yp,
         'exact_B_u0':B0_exact,'exact_B_u1':B1_exact,'exact_c_star':cstar,
         'Tricomi_integral_slope_difference':abs(cstar-cref_integral),
         'limiting_infinity_checks':[],'fixed_slope_runs':[],'zero_curve_runs':[],
         'qualification':'Analytic Kummer limit and implicit-function structure are exact; computed finite-epsilon zero locations are numerical with refinement checks, not interval-certified enclosures.'}
    configs=[dict(R=48,N=24,taylor_order=44,step_fraction='.20'),
             dict(R=80,N=40,taylor_order=56,step_fraction='.18'),
             dict(R=112,N=60,taylor_order=64,step_fraction='.16')]
    for cfg in configs:
        state,diag=infinity_state(p0,**cfg)
        y,yp=algebraic_state_from_m(p0,z,state)
        B0=coefficient_from_states(p0,z,(ordinary[0,0],ordinary[1,0]),state)
        B1=coefficient_from_states(p0,z,(ordinary[0,1],ordinary[1,1]),state)
        rr={'settings':cfg,'yalg_error':abs(y-exact_y),'yalg_prime_error':abs(yp-exact_yp),
            'B_u0':B0,'B_u1':B1,'B_u0_error':abs(B0-B0_exact),'B_u1_error':abs(B1-B1_exact),
            'c_star_from_generic_transport':-B0/B1,'c_star_error':abs(-B0/B1-cstar),**diag}
        zR=-mp.mpf(cfg['R']); logzR=mp.log(-zR)+mp.j*mp.pi
        yyR=pref*mp.hyperu(aa,bb,alpha*(1-zR))
        ypR=pref*alpha*aa*mp.hyperu(aa+1,bb+1,alpha*(1-zR))
        mR=mp.exp(aa*logzR)*yyR
        mzR=mp.exp(aa*logzR)*(ypR+aa*yyR/zR)
        PN,PNp=polynomial_state(asymptotic_coefficients(p0,cfg['N']),zR)
        rr['exact_terminal_m_error']=abs(mR-PN)
        rr['exact_terminal_mz_error']=abs(mzR-PNp)
        assert rr['exact_terminal_m_error']<=diag['terminal_tail_bounds']['terminal_m_error_bound']
        assert rr['exact_terminal_mz_error']<=diag['terminal_tail_bounds']['terminal_mz_error_bound']
        out['limiting_infinity_checks'].append(rr)
        print('limit',cfg,'yerr',mp.nstr(rr['yalg_error'],8),'cstarerr',mp.nstr(rr['c_star_error'],8),flush=True)
        Path('analysis/generic_che_resonant_zero.json').write_text(json.dumps(record(out),indent=2)+'\n')
    settings=dict(R=80,N=40,taylor_order=56,step_fraction='.18')
    for eps in [mp.mpf('1e-3'),mp.mpf('1e-5'),mp.mpf('1e-7')]:
        val=coefficient(parameters(eps,cstar),**settings)
        out['fixed_slope_runs'].append({'epsilon':eps,'t':cstar,'amplitude':val,'amplitude_abs':abs(val),'amplitude_abs_over_epsilon':abs(val)/eps,'settings':settings})
        print('fixed slope',eps,'Babs',mp.nstr(abs(val),10),flush=True)
        calls=[]
        def f(t):
            value=coefficient(parameters(eps,t),**settings)
            calls.append({'t':t,'amplitude':value})
            return value
        root=mp.findroot(f,(cstar,cstar+mp.mpf('.01')),solver='secant',tol=mp.mpf('1e-22'),maxsteps=12,verify=True)
        residual,diag=coefficient(parameters(eps,root),diagnostics=True,**settings)
        # Finer terminal and transport settings assess, rather than assume, root residual.
        higher=dict(R=112,N=60,taylor_order=64,step_fraction='.16')
        highres=coefficient(parameters(eps,root),**higher)
        out['zero_curve_runs'].append({'epsilon':eps,'root_t':root,'root_t_minus_cstar':root-cstar,
            'accessory_B':eps*root,'eta':parameters(eps,root)[-1],'amplitude_at_root':residual,
            'higher_settings':higher,'higher_amplitude_at_same_root':highres,
            'arithmetic_settings':settings,'secant_calls':calls,'diagnostics':diag})
        print('zero',eps,'t',mp.nstr(root,28),'finerB',mp.nstr(abs(highres),8),flush=True)
        Path('analysis/generic_che_resonant_zero.json').write_text(json.dumps(record(out),indent=2)+'\n')
    print('cstar=',mp.nstr(cstar,45),flush=True)
    return out

if __name__=='__main__':main()
