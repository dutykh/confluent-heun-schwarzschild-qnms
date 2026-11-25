# Generic CHE connection-amplitude numerical audit

This note covers `generic_che_infinity_numeric.py`, `generic_che_forward_numeric.py`, and `generic_che_resonant_zero.py`, with their corresponding JSON files. The path-dependence test is separately preserved as `generic_che_path_verify.py/.json` and `generic_che_path_notes.md`. The stale duplicate `generic_che_resonance_verify.json` was removed only after verifying byte-for-byte equality with the preserved path JSON. The other agent's frame/resonance files are separate and were not modified.

## Conventions, parameter choices, and continuation paths

The equation is

H''+[alpha+(beta+1)/z+(gamma+1)/(z-1)]H'+(Az+B)H/[z(z-1)]=0,
A=delta+alpha(beta+gamma+2)/2,
B=eta+beta/2+(gamma-alpha)(beta+1)/2.

The computed amplitude is denoted B_exp here to distinguish it from the accessory numerator B. It is the coefficient of a **unit-leading growing exponential branch** on the specified recessive chart for the algebraic branch. Adding an algebraic solution to the growing solution does not change B_exp.

Set rho=-A/alpha, u=beta+1+2rho, v=gamma+1, C=rho(rho+beta), D=A+B+rho(gamma+1). Then y_alg=z^rho*m and

m''+p*m'+V*m=0,
p=alpha+u/z+v/(z-1), V=C/z^2+D/[z(z-1)].

The two unrestricted complex examples, in order (alpha,beta,gamma,delta,eta), are:

| Case | alpha | beta | gamma | delta | eta |
|---|---|---|---|---|---|
| A | -1.2+0.3i | 0.3+0.2i | -0.4+0.1i | 0.2-0.15i | 0.6+0.05i |
| B | 0.8+0.6i | -0.2+0.15i | 0.4-0.25i | 0.35+0.2i | -0.1+0.4i |

Both are off Konno's factorization locus. Their (delta-delta0,eta-eta0) differences are respectively (1.835-0.495i,0.21+0.20i) and (-0.09-0.38i,-0.10375+0.845i).

For each example c=-conj(alpha)/|alpha| is fixed, so Re(alpha*c)=-|alpha|<0. The backward path starts at z=cR, moves radially to |z|=2, follows a homotopic polygon approximating the circle from arg(c) to pi/4, and moves radially to z*=0.5+0.5i. Starting log z is log R+i arg(c), with arg(c) in (0,2pi] (the positive real ray is assigned2pi). Log(z-1)-log(z) tends to zero at infinity. At z* both logs are principal, log(z-1) having argument3pi/4. All further matching checks remain in the upper half-plane. The forward method traverses the reverse path. A circle of radius3 provides a homotopic path check. These choices are part of the answer: a different winding/continuation may give different connection data.

For parameter families, a fixed ray from the parameter-chart center should be supplied; it must retain Re(alpha*c)<0. Alpha is fixed throughout every tested family here. The infinity routines allow an explicitly fixed ray. Alpha=0 is excluded.

## Backward method and checks

The finite polynomial PN=sum_(n=0)^N a_n z^(-n), a0=1, uses

a_(n+1)={[n(n+1-u-v+alpha)+C+D]a_n-[(n-1)(n-u)+C]a_(n-1)}/[alpha(n+1)].

Its exact residual is

E PN=(E_N z+F_N)/[z^(N+2)(z-1)],
E_N=[N(N+1-u-v+alpha)+C+D]a_N-[(N-1)(N-u)+C]a_(N-1),
F_N=-[N(N+1-u)+C]a_N.

No convergence of the infinite asymptotic series is assumed. The polynomial and its derivative provide finite terminal data; arbitrary-precision Taylor continuation propagates m,m_z backward. The local Taylor coefficients are computed from the rational first-order system. Steps are bounded by a fraction of the distance to the nearest finite singularity and by1/|alpha|.

The amplitude is

B_exp=exp(alpha*z)z^(beta+1)(z-1)^(gamma+1)W(H0,y_alg)/alpha,

where H0 is evaluated by its independent convergent origin recurrence. The main campaign used75 decimal digits and (R,N,Taylor order,step fraction)=(24,16,36,.22),(40,24,44,.20),(64,32,52,.18),(96,48,64,.16). At the final setting the amplitude changes between z=.5+.5i,.4+.3i,.6+.2i were below7.6e-51. Increasing Taylor order to76, reducing the step fraction to.12, and changing the circle radius to3 changed the results by2.23e-54 (A) and2.31e-51 (B). These are empirical transport checks, not interval enclosures.

## Forward method: only the regular solution is transported

Define

q=J*z^(-rho)*H0,
J=exp(alpha*z)z^u(z-1)^v.

It satisfies

q''-p*q'+W*q=0,
W=V-p'=(C+u)/z^2+D/[z(z-1)]+v/(z-1)^2.

The initial q,q' at z* are obtained directly from H0,H0'. Propagating this equation forward to cR needs no numerical infinity solution. The finite coefficient approximant is

B_N=[(p*PN+PN')q-PN*q']/alpha.

Exactly,

d/dz{(p*PN+PN')q-PN*q'}=q*E(PN),
B_exp-B_N=alpha^(-1) integral_(cR)^(infinity_c) q*E(PN) dz.

For sigma=-Re(alpha*c)-|u|/R-|v|/(R-1)>0, let

F_R=R/(R-1)[|E_N|R^(-N-1)/(N+1)+|F_N|R^(-N-2)/(N+2)],
W_R=|C+u|/R+|D|log(R/(R-1))+|v|/(R-1).

Then a valid exact-state tail bound is

|B_exp-B_N| <= F_R/|alpha| * (|q(R)|+|q'(R)|/sigma)*exp(W_R/sigma).

The numerical evaluations of that bound use computed endpoint states. Their transport and arithmetic errors have not been enclosed by interval arithmetic, so this is **not a certified total computer error bound**.

With exact transport the forward and backward PN approximants are algebraically equal. Their agreement therefore checks the distinct propagation equations/directions, gauge factors, and phases, not an independent model of the asymptotic truncation. The observed R96 agreement was about5.6e-46, limited by the45-digit stored backward values.

The independent forward truncation refinement to R128,N64 (Taylor order80,step fraction.14) is more informative:

| Case | Observed coefficient difference | Exact-state tail bound at96,48 | Ratio | Tail bound at128,64 |
|---|---:|---:|---:|---:|
| A | 2.4310043116e-42 | 2.5037704148e-42 | .97093739 | 1.1120898546e-55 |
| B | 2.3950153015e-35 | 2.544090334e-35 | .94140340 | 4.6467149091e-47 |

At fixed R128,N64, simultaneously raising precision75 to95 digits, raising Taylor order80 to92, reducing step fraction.14 to.12, and using the homotopic radius3 arc changed the coefficient by5.54e-71 (A) and3.78e-68 (B). These measured propagation changes lie far below the reported tail bounds.

Thus the observed truncation changes are consistent with the analytic bounds. These differences are not themselves exact errors. Final numerical amplitudes (extra digits retained in JSON) are

A: -0.306015253608359116222223426391789609630631748
   +0.575567504362596887600953170198668974089652972 i,
B: +0.553776247131803972838220847997246412619512908
   +0.844387692756856058900040384060636097621273035 i.

## Exact apparent-resonance check and a resolved accessory zero

Fix alpha=.7,gamma=.4,delta=.2 and beta=-1+epsilon,
B=epsilon*t, eta=1/2+epsilon*(t-7/20).
At epsilon=0 the origin becomes ordinary and A=.69. With a=69/70,b=7/5,w=alpha(1-z), the equation reduces exactly to Kummer's equation. For the negative-real starting ray c=-1 and log z=log r+i*pi,

y_alg=alpha^a exp(-i*pi*a) U(a,b,alpha(1-z)).

This supplies an independent full normalization and phase check of the generic infinity solver. For R,N=(48,24),(80,40),(112,60), the terminal m errors against this exact Tricomi expression were2.0471e-14,1.0908e-23,1.8013e-33, respectively, versus analytic terminal bounds4.3314e-14,2.1159e-23,3.4693e-33. Both value and derivative bounds passed. The errors in y_alg at the matching point were3.26e-14,1.71e-23,2.89e-33.

Let u0,u1 have ordinary-point initial data(1,0),(0,1). Write B0=B_exp[u0],B1=B_exp[u1]. Exact values are approximately

B0=.1911741690115361927004719344813542+.6927040442799373803404042862363335i,
B1=-.2557672092251640124843368294197407-.9267516691220602857604628221265342i.

The zero-amplitude limiting slope is exactly

t_star=-B0/B1=alpha*a*U(a+1,b+1,alpha)/U(a,b,alpha)
      =.747453786553367359040558990081093331866319248...

The positive Tricomi integral for a>0,alpha>0 proves U(a,b,alpha)>0, hence B1!=0. An independently evaluated positive integral reproduces the slope. The generic transport ratio at R112,N60 agrees with it to2.6e-51. Its ratio is more accurate than the separate amplitudes because common normalization error cancels; that improvement is not assumed for generic amplitudes.

In the resolved chart(epsilon,t), the normalized origin solution is jointly holomorphic, with H0(0,t)=u0+t*u1. The infinity construction is holomorphic on the fixed chart. Therefore F(epsilon,t)=B_exp[H0] is holomorphic, F(0,t)=B0+t B1, and partial_t F(0,t_star)=B1!=0. The holomorphic implicit-function theorem yields a unique local zero curve t(epsilon). This is an exact existence statement, distinct from the numerical locations below.

Using the scalar gamma-regularized solution without resolving its vanishing normalization would obscure this conclusion: B_exp[Phi]=F(epsilon,t)/Gamma(epsilon) vanishes at epsilon=0 for every t because Phi is then the zero solution. That trivial normalization zero is not a connection eigencondition. The resolved H0 frame gives the nontrivial condition B0+t B1=0.

| epsilon | computed t(epsilon) | finer-setting amplitude at the same root |
|---:|---:|---:|
| 1e-3 | .7477441452401204816551441447 | 1.0172e-41 |
| 1e-5 | .7474566901786279866189859707 | 1.0166e-41 |
| 1e-7 | .7474538155896238054048569455 | 1.0166e-41 |

Roots used R80,N40, Taylor56,step.18; the residual check used R112,N60,Taylor64,step.16. Their computed imaginary parts are about9e-42 and are numerical error at this scale. Keeping t=t_star instead gives amplitudes proportional to epsilon, with |F|/epsilon approaching .279154. The accessory numerator is B=epsilon*t and eta=1/2+epsilon*(t-.35). This illustrates a generic CHE accessory-parameter connection zero through the resolved apparent resonance. It is not a QNM identification and is not presented as a new pole-skipping phenomenon.

## Reproduction and limits

Run from the project directory:

    python3 analysis/generic_che_path_verify.py
    python3 analysis/generic_che_infinity_numeric.py --campaign
    python3 analysis/generic_che_forward_numeric.py
    python3 analysis/generic_che_resonant_zero.py

Dependencies are mpmath and sympy. Path tests use100 digits; infinity/forward/zero campaigns use75 digits. The files contain executable recurrences, Taylor propagation, all settings, branches, and diagnostics. The calculations establish internal consistency on the tested parameter charts and an exact special-limit check. They do not establish all Stokes connection matrices, global completeness, or interval-certified computer enclosures. No novelty or priority claim is made for the underlying constructions.
