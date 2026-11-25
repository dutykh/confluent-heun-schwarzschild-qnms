# Generic CHE interval coefficient certificates

`generic.py` certifies two nonresonant complex-parameter examples already studied by floating-point continuation. It uses a single forward solution of the adjoint equation, interval rational arithmetic, and MPFI enclosures of elementary functions. No infinity solution is propagated. The final result includes all three errors: origin-series truncation, finite-path Taylor continuation (including arithmetic), and the infinite asymptotic tail. Decimal strings in the JSON are displays; the exact rational lower and upper endpoints define the certificate.

## Equation and normalization

The original equation is

\[
H''+\left(\alpha+\frac{\beta+1}{z}+\frac{\gamma+1}{z-1}\right)H'
+\frac{Az+B}{z(z-1)}H=0,
\quad H(0)=1,
\]

with the Maple parameter dictionary

\[
A=\delta+\tfrac\alpha2(\beta+\gamma+2),\quad
B=\eta+\tfrac\beta2+\tfrac12(\gamma-\alpha)(\beta+1).
\]

Put rho=-A/alpha, u=beta+1+2rho, v=gamma+1,
C=rho(rho+beta), D=A+B+rho(gamma+1). The original algebraic infinity solution is z^rho m(z), with m asymptotic to 1. Its equation has p=alpha+u/z+v/(z-1) and V=C/z²+D/[z(z-1)]. The forward unknown is

\[
q=e^{\alpha z}z^{u-\rho}(z-1)^v H,
\qquad q''-p q'+Qq=0,
\quad Q=V-p'=\frac{C+u}{z^2}+\frac{D}{z(z-1)}+\frac{v}{(z-1)^2}.
\]

For the canonical exponential connection coefficient b, q tends to b along the selected damping ray. The coefficient is obtained from the finite concomitant and the residual integral as in the parent derivation.

## Certified seed

All input decimal parameters are parsed as exact rational numbers and then enclosed. The regular Frobenius coefficients obey

\[
a_{n+1}=\frac{[n(n+\beta+\gamma+1-\alpha)+B]a_n+
[\alpha(n-1)+A]a_{n-1}}{(n+1)(n+\beta+1)},
\quad a_{-1}=0,\ a_0=1.
\]

Let c0=1+Re(beta)>0, K=|beta+gamma+1-alpha|, L=|B|,
U=|alpha| and T=|A-alpha|. The two recurrence multipliers are bounded by
S=max(1+K+L/c0,U+T/c0). Taking M=S+1 gives |a_n|<=M^n by induction, since M²-SM-S=1. In actual code upper enclosures are used for S,M and a lower enclosure for c0.

For r=M|z0|<1 and truncation through degree J,

\[
|H-H_J|\le \frac{r^{J+1}}{1-r},\qquad
|H'-H_J'|\le M r^J\frac{J+1-Jr}{(1-r)^2}.
\]

Here z0=(1+i)/64, J=220, and the code verifies r<1/2. Both errors are added as complex rectangles. MPFI principal logs at z0 and z0-1 determine the gauge and its derivative; both inputs lie strictly above the branch cut. The complete seed enclosure is then propagated.

## Finite transport and Cauchy bound

At an exact rational center zc choose

\[
d=\tfrac12\min(\|z_c\|_\infty,\|z_c-1\|_\infty,2),
\quad r_j=\|z_c-j\|_\infty-d>0.
\]

These rational infinity norms underestimate Euclidean distances, so the closed disk of radius d contains neither singularity. Throughout that disk,

\[
|p|\le P=|\alpha|+|u|/r_0+|v|/r_1,
\quad |Q|\le Q_*=|C+u|/r_0^2+|D|/(r_0r_1)+|v|/r_1^2.
\]

For the scaled first-order state Y=(q,dq'), the matrix infinity norm is bounded by max(1/d,dQ_*+P). Thus radial Grönwall bounds its norm throughout the disk by

\[
Z e^{K_d},\quad
Z=\max(|q(z_c)|,d|q'(z_c)|),\quad
K_d=\max(1,dP+d^2Q_*).
\]

Initial-state uncertainty is included when computing Z. Cauchy's coefficient estimate gives, for either component's Taylor polynomial of degree L and tau>=|step|/d<1, the bound

\[
E=Z e^{K_d}\frac{\tau^{L+1}}{1-\tau}.
\]

The q remainder is at most E; the q' remainder is at most E/d. The code uses the exact rational l1 norm of the step as the upper bound and verifies tau<=1/4. It computes q coefficients through degree L+1, evaluates q through degree L and q' through degree L, and adds these remainder boxes. Roundoff and propagation of prior uncertainty are enclosed by interval arithmetic.

The efficient Taylor recurrence comes from

\[
A_0(z)q''+B_0(z)q'+C_0(z)q=0,
\]

where A0=z²(z-1)²,
B0=-alpha A0-u z(z-1)²-v z²(z-1),
C0=(C+u)(z-1)²+D z(z-1)+v z².
It has at most four backward lags. Division by the shifted constant A0(zc) is valid because zc is different from 0 and 1.

## Paths and infinity tail

The seed first travels through (1+i)/2 and 2+2i. Case A then passes through 4+2i, 4, and +R; case B passes through -2+2i, -4, and -R. Exact dyadic subdivisions lie along these segments. The paths are in the same homotopy classes as the previous upper-half-plane continuations. Case A's ray angle moves from arg(4+i) to 0; case B's from arg(-4+3i) to pi. These rotations stay inside Re(alpha z)<0 and do not change the canonical infinity solution. There is no extra winding around either regular singularity; in particular case A uses argument 0 at positive infinity, not 2pi.

For the degree-N inverse-power polynomial PN, let RN=LPN=(E_N z+F_N)/[z^(N+2)(z-1)]. The code evaluates exactly the two residual coefficients via the recurrence. For |z|=R>1 and unit ray c=+1 or -1, define

\[
\sigma=-\Re(\alpha c)-|u|/R-|v|/(R-1)>0,
\]

\[
F_R=\frac{R}{R-1}\left[\frac{|E_N|}{(N+1)R^{N+1}}+
\frac{|F_N|}{(N+2)R^{N+2}}\right],
\]

\[
W_R=\frac{|C+u|}{R}+|D|\log\frac{R}{R-1}+\frac{|v|}{R-1}.
\]

The parent forward-coefficient theorem gives

\[
\left|b-\frac{(pP_N+P_N')q-P_Nq'}{\alpha}\right|
\le\frac{F_R}{|\alpha|}
\left(|q(Rc)|+\frac{|q'(Rc)|}{\sigma}\right)e^{W_R/\sigma}.
\]

All quantities on the right use outward enclosures of the already certified endpoint state. This disk is added to the finite concomitant rectangle. This proof applies to the two specified parameter tuples, with the explicit positive sigma test; it is not a certificate at resonance, alpha=0, a Stokes boundary, or an arbitrary frequency rectangle.

Run: `python cert/generic.py --R 64 --N 32 --degree 192 --dps 160`.

The trusted numerical base is mpmath's outward rational interval operations, MPFI/MPFR elementary operations, and the exact-endpoint bridge in elementary.py. No `iv.findroot` result is involved.

## Completed numerical certificates

At R=64, N=32, Taylor degree 192 and 160 decimal working digits:

| Case | Certified real interval (outward decimal endpoints) | Certified imaginary interval (outward decimal endpoints) |
|---|---|---|
| A | [-0.3060152536083591162222234264723686128426, -0.3060152536083591162222234263470880920160] | [0.5755675043625968876009531700788881294812, 0.5755675043625968876009531702041686503078] |
| B | [0.5537762471318039728382188972524773924480, 0.5537762471318039728382471409990674875430] | [0.8443876927568560589000227480816169483819, 0.8443876927568560589000509918282070434769] |

Case A used 471 Taylor steps; its finite-concomitant component widths were approximately 4.24e-62, and its final component half-widths are below 6.265e-29. Case B used 491 steps; finite-concomitant component widths were approximately 2.48e-46 and final component half-widths are below 1.413e-23. The infinite-tail estimate dominates both final rectangles. The prior 45-digit floating-point values lie inside both certificates; this comparison checks normalization and paths and is not used to prove the enclosures.

The two interval runs took about 38 and 41 seconds on the current machine. A separate mathematical audit checked the polynomial recurrence, scaled Cauchy majorant, residual integral and normalization. Unsupported N=0 is explicitly rejected. These are certificates for two finite connection coefficients, not root certificates or a general symbolic connection matrix.

Implementation arithmetic reference: N. Revol and F. Rouillier, Motivations for an Arbitrary Precision Interval Arithmetic and the MPFI Library, Reliable Computing 11 (2005), 275–290, DOI 10.1007/s11155-005-6891-y. The elementary bridge is separately exercised by elementary.py self-tests.
