# Joint resonant limit of the normalized generic CHE solution

Conventions: H''+[alpha+(beta+1)/z+(gamma+1)/(z-1)]H'+(Az+B)H/[z(z-1)]=0, with A=delta+alpha(beta+gamma+2)/2, B=eta+beta/2+(gamma-alpha)(beta+1)/2. All local real-interval powers use real logarithms. The finite connection convention is Y0=Y1 C01, with solution columns and their z derivatives.

## Recurrence and gamma regularization

For H0=sum a_n z^n, a_0=1, a_{-1}=0,

(n+1)(n+beta+1)a_(n+1)=[n(n+beta+gamma+1-alpha)+B]a_n+[alpha(n-1)+A]a_(n-1).

In particular H0'(0)=B/(beta+1). Set p_n=n!(beta+1)_n a_n. Then p_0=1, p_1=B, and

p_(n+1)=[n(n+beta+gamma+1-alpha)+B]p_n+n(n+beta)[alpha(n-1)+A]p_(n-1).

Every p_n is polynomial in the equation parameters. The regularized solution is

Phi=H0/Gamma(beta+1)=sum p_n z^n/[n! Gamma(beta+n+1)].

Its coefficients are entire in the parameters. Local uniform convergence on compact parameter sets and |z|<=r<1 follows from a polynomial majorant: beyond a parameter-uniform index, the original coefficient recurrence has denominator bounded below by a fixed multiple of n^2, its a_n multiplier is 1+O(1/n), and its a_(n-1) multiplier is O(1/n). Thus the maximum of successive coefficient magnitudes grows at most polynomially. Finitely many regularized initial coefficients are bounded on each compact parameter set. This proves a jointly holomorphic regularized local solution, without dividing by parameter factors that vanish at resonance.

## Explicit path dependence

Fix alpha=7/10, gamma=2/5, delta=1/5, and set

beta=-1+epsilon, B=kappa,
eta=1/2+kappa-7epsilon/20,
A=69/100+7epsilon/20.

Along kappa=c epsilon, H0(0)=1 and H0'(0)=c exactly. At epsilon=kappa=0 the origin becomes ordinary, since the equation reduces to

u''+[0.7+1.4/(z-1)]u'+0.69/(z-1)u=0.

Let u0,u1 be its ordinary-point solutions with initial data (1,0) and (0,1), respectively. In the normalized recurrence, a_1=c and every subsequent denominator (n+1)(n+epsilon), n>=1, is nonzero near epsilon=0. The local coefficient family is holomorphic there and converges locally uniformly, giving

H0 -> u0+c u1.

This is an analytic statement; numerical checks are supplementary. At z=.3:

u0=1.03845223794182439010363075806870056,
u1=.34702636137952558033193565806914601.

The c=0,1,-2 paths give different limits 1.03845223794182439,1.38547859932134997,.344399515182773229. All have value one at the origin.

The second Frobenius column is z^(-beta)HeunC(alpha,-beta,gamma,delta,eta;z). It tends to u1, since its limiting value at zero is zero and derivative is one. Therefore, with U the ordinary-point fundamental matrix,

Y0^(c) -> U T_c, T_c=[[1,0],[c,1]].

For every common invertible limiting target basis Z, its connection matrix consequently tends to Z^(-1)U T_c. Changing the target basis cannot remove the dependence on c. This applies in particular to a sectorial infinity basis whenever that basis has a common invertible parameter limit. No numerical infinity calculation is claimed by this finite-point test.

## Actual finite-point matrix check

The nonresonant basis at z=1 is evaluated via the reflection

(alpha,beta,gamma,delta,eta) -> (-alpha,gamma,beta,-delta,eta+delta), w=1-z,

and its second exponent -gamma. Both columns are regular in the chosen parameter paths. At the limiting parameter point,

Cref=Y1^(-1)U approximately
[[-.71249366180062113372873264018010, -2.46612476929139951087256045358403],
 [1.10243928914095154985942055923941, 2.07340732600842901838433816084007]].

The code verifies C01 -> Cref T_c for the three paths with epsilon=1e-2,1e-4,1e-8,1e-12,1e-20. At the final epsilon, the maximum matrix errors are respectively2.17e-20,6.01e-20,1.55e-19. Matching-point independence at z=.3,.5,.7 is checked to8.6e-101 for a generic finite epsilon; det C01=-beta exp(-alpha)/gamma is checked to1.0e-100. The ordinary reference solutions agree with an independent mpmath Taylor ODE solver to6.7e-86.

## What gamma regularization does and does not fix

At epsilon=kappa=0, Phi identically vanishes: its constant coefficient contains1/Gamma(0)=0, p1=B=0, the n=1 recurrence has n(n+beta)=0, and all subsequent p_n vanish. Consequently multiplying by Gamma(epsilon) does not define a path-independent joint extension of H0.

There is a constructive replacement. In the explicit parameter chart (epsilon,kappa),

Phi=epsilon u0+kappa u1+terms of total degree at least two.

Thus the two parameter derivatives partial_epsilon Phi and partial_kappa Phi at the intersection provide u0 and u1. Differentiating the equation is legitimate because Phi itself is zero at the intersection, so both derivative functions solve the limiting homogeneous equation. This recovers a useful ordinary/resonant frame, but depends on stating the chart or another normalization convention.

The example does not prohibit choosing a fixed continuation prescription or a physical one-parameter curve. It disproves a joint, path-independent normalized HeunC extension based solely on H0(0)=1 at this degenerate parameter point.

Executable: generic_che_path_verify.py. Results: generic_che_path_verify.json. Dependencies mpmath and sympy. Working precision100 decimal digits. This audit is distinct from generic_che_frame_verify.py and its companion frame notes.
