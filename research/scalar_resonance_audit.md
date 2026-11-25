# Scalar Schwarzschild horizon resonances and safe connection coefficients

Independent audit for the root agent. No main report was edited.

## Fixed conventions

The equation supplied by the root is

\[
u_{xx}+[\Omega^2-U(z)]u=0,\qquad
x=z+\log(z-1),\qquad f=(z-1)/z,
\quad U=f[L/z^2+1/z^3],\quad L=\ell(\ell+1).
\]

The normalized horizon solution is

\[
u_H(\Omega,z)=e^{-2i\Omega}e^{i\Omega z}z(z-1)^{-i\Omega}
 H_H(\Omega,1-z),
\]

with

\[
H_H=\operatorname{HeunC}(-2i\Omega,-2i\Omega,0,-2\Omega^2,2\Omega^2-L;t),\quad t=1-z.
\]

Thus u_H/e^{-iΩx}→1 at the horizon. Set λ=2iΩ, c=1−λ. The exceptional nonzero negative imaginary frequencies are Ω_N=−iN/2, N=1,2,…, at which λ=N and c=1−N. Frequency zero is a separate problem, not covered by the following formulas involving 1/(2iΩ).

## 1. Exact horizon recurrence

Writing H_H=Σ h_n t^n, h_0=1, h_{−1}=0, direct substitution into the canonical CHE gives

\[
n(n-\lambda)h_n=
[n(n-1)-\lambda^2-L]h_{n-1}
-\lambda(n-1-\lambda)h_{n-2}. \tag{1}
\]

Equivalently in Ω,

\[
n(n-2i\Omega)h_n=
[n(n-1)+4\Omega^2-L]h_{n-1}
-[2i\Omega(n-1)+4\Omega^2]h_{n-2}.
\]

For Ω=Ω_N, coefficients h_n for n<N are well defined. At the Nth step the compatibility obstruction is

\[
R_N=-(N+L)h_{N-1}+N h_{N-2}. \tag{2}
\]

If R_N≠0, the normalized horizon family has a simple pole (an anomalous point). If R_N=0, its Nth recurrence is formally 0/0; one must calculate the parameter limit, not choose h_N arbitrarily. This is the miraculous case. For a general one-frequency family, the correctly continued coefficient there is R_N'(Ω_N)/(−2iN), where the derivative includes the dependence of earlier coefficients and equation parameters on Ω.

## 2. For this scalar problem there are NO miraculous horizon frequencies Ω_N

This conclusion can be proved for every real L≥0, so in particular for every physical scalar angular momentum. Define the N×N real tridiagonal matrix A_N(L), indexed k=0,…,N−1, by

\[
(A_N)_{kk}=N^2+L-k(k+1),
\]
\[
(A_N)_{k,k+1}=-(k+1)(N-k-1),\qquad
(A_N)_{k,k-1}=-N(N-k).
\]

A compatible recurrence through order N would give a nonzero nullvector (h_0,…,h_{N−1}) of A_N. But its diagonal entries are positive and it is strictly COLUMN diagonally dominant:

\[
(A_N)_{kk}-\sum_{j\ne k}|(A_N)_{jk}|=L+N-k>0. \tag{3}
\]

Hence A_N is nonsingular. It is also a nonsingular M-matrix; alternatively, the off-diagonal signs and strict column Gershgorin dominance put all eigenvalues in the open right half-plane. In either argument

\[
\Delta_N(L):=\det A_N(L)>0,\qquad
R_N=-\frac{\Delta_N(L)}{[(N-1)!]^2}\ne0. \tag{4}
\]

The determinant formula follows by the usual tridiagonal determinant recurrence, or from solving the first N−1 rows successively. Examples are

\[
\Delta_1=L+1,
\]
\[
\Delta_2=L^2+6L+6,
\]
\[
\Delta_3=L^3+19L^2+93L+99,
\]
\[
\Delta_4=L^4+44L^3+604L^2+2832L+3216.
\]

Therefore **all Ω_N are genuine simple poles of the unit-amplitude scalar horizon solution**. It would be incorrect to claim that its gamma regularization automatically introduces a spurious zero at each Ω_N. The gamma zero cancels an actual pole here.

This is a scalar-specific proof. The gravitational Regge–Wheeler equation has a miraculous algebraically special point, so the general warnings from gravitational literature cannot be transplanted without checking the scalar recurrence.

## 3. Division-free gamma-regularized series

Set

\[
p_n=n!(1-\lambda)_n h_n.
\]

The p_n can be calculated without frequency-dependent divisions:

\[
p_0=1,\quad p_{-1}=0,\qquad
p_n=[n(n-1)-\lambda^2-L]p_{n-1}
-\lambda(n-1)(n-1-\lambda)^2p_{n-2}. \tag{5}
\]

The regularized local solution is

\[
\widehat H_H(\Omega,t)=\frac{H_H(\Omega,t)}{\Gamma(1-2i\Omega)}
=\sum_{n=0}^{\infty}\frac{p_n(2i\Omega,L)}{n!\,\Gamma(n+1-2i\Omega)}t^n. \tag{6}
\]

Use reciprocal gamma directly. Do NOT evaluate a diverging HeunC and multiply its numerical result by zero. The expression converges for |t|<1 and extends locally analytically in Ω across the exceptional points; its coefficients and Ω derivatives can be computed using the polynomial recurrence and reciprocal gamma.

At λ=N all terms n<N vanish and

\[
p_N(N,L)=(-1)^N\Delta_N(L),\quad
\widehat H_H(\Omega_N,t)=\frac{(-1)^N\Delta_N(L)}{N!}t^N+O(t^{N+1}). \tag{7}
\]

In particular the regularized solution is NONZERO. A convenient matching point is z in (1,2), so |t|<1; continuation beyond that disk uses an ordinary nonsingular ODE initial-value step.

## 4. Exact upper-half-plane identity at the exceptional frequencies

Let \(\widehat u_H=u_H/\Gamma(1-2i\Omega)\). The scalar differential equation is even in Ω, and at Ω_N the regularized solution has the small horizon exponent +N/2. Comparing its leading coefficient to the normalized horizon solution at −Ω_N=+iN/2 gives

\[
\boxed{
\widehat u_H(\Omega_N,z)
=\frac{e^{-N}\Delta_N(L)}{N!}\,u_H(-\Omega_N,z).
} \tag{8}
\]

This holds identically in z on the selected spatial sheet by uniqueness of the small Frobenius solution. It supplies a particularly safe way to evaluate the regularized horizon vector exactly on the resonance grid, using an ordinary upper-half-plane solution with nonsingular recurrence. For N=1, ell=0, the factor is e^{-1}.

The fact that the regularized continued boundary solution has the *other* leading exponent at the exceptional frequency is not a change of physical prescription: outgoing/ingoing boundary families are defined by analytic continuation from their initial domain. Pointwise exponent tests at a resonance can misclassify them.

## 5. Which determinant zeros are physical?

Let g=u_out be the selected outgoing infinity solution, continued to a specified local frequency sheet, and define

\[
D=W_x(u_H,g)=2i\Omega A_{\mathrm{in}},\quad
\widehat D=W_x(\widehat u_H,g)=D/\Gamma(1-2i\Omega).
\]

The Green kernel (up to the Wronskian sign convention) is

\[
G(\Omega;x,y)=\frac{u_H(x_<)g(x_>)}{D}
=\frac{\widehat u_H(x_<)g(x_>)}{\widehat D}. \tag{9}
\]

**For the present scalar family, at every Ω_N the horizon vector \(\widehat u_H\) is nonzero. Provided the selected infinity vector is analytic and nonzero there, zeros of \(\widehat D\) are legitimate poles of the continued Green kernel, with their order determined by the zero order of \(\widehat D\).** They are not automatic gamma artifacts. This does not assert that such zeros actually occur; the local horizon analysis does not settle that spectral question.

At such a point the unregularized D can be finite and nonzero while u_H has a pole. Therefore a criterion based only on zeros of the raw meromorphic A_in can miss an exceptional Green-function pole. Use the determinant of nonzero analytic endpoint sections, not an arbitrary singular normalization.

For a GENERAL CHE or a gravitational equation with a miracle, u_H may remain finite and nonzero at Ω_N. Then multiplication by 1/Γ(1−2iΩ) makes \(\widehat u_H\) the zero solution at that frequency, and \(\widehat D\) has an automatic zero. Divide out the common vanishing factor of the endpoint vector before diagnosing spectral zeros.

A rigorous implementation rule is to represent each endpoint by its analytic initial-data vector (u,u_x) at the matching point. If both components vanish to common order m, divide the entire vector by (Ω−Ω_*)^m. Do this for both endpoints. The determinant of the resulting nonzero analytic vectors has the correct spectral zero set. The same cancellation follows directly from the numerator and denominator of (9). Dividing by a nonvanishing analytic component on a local patch is an additional harmless conditioning choice; dividing by a Euclidean norm is nonholomorphic and is unsuitable for holomorphic derivative or argument-principle calculations.

## 6. A warning about Laurent finite parts

Suppose L(Ω)u_H(Ω)=0 and

\[
u_H(\Omega)=a_{-1}/\varepsilon+a_0+O(\varepsilon),\quad
L(\Omega)=L_0+\varepsilon L_1+\cdots,\quad\varepsilon=\Omega-\Omega_*.
\]

Then

\[
L_0a_{-1}=0,\qquad L_0a_0=-L_1a_{-1}. \tag{10}
\]

**The Laurent finite part a_0 is generally NOT a solution of the homogeneous equation at Ω_*.** Thus it cannot simply replace the horizon solution in a Wronskian. If a homogeneous logarithmic companion is needed, first subtract the pole times a suitable analytic *solution family at the same frequency*, such as the opposite-horizon family, and only then take the limit. This produces a homogeneous solution but changes the endpoint boundary family unless the subtracted term is tracked explicitly. At a miracle the correct free resonant Frobenius coefficient is determined by the frequency limit (l'Hopital/parameter derivative), not by an arbitrary choice.

## 7. Infinity Volterra construction and the negative imaginary axis

For the outgoing normalization write g=e^{iΩx}v. On an admissible complex contour to infinity,

\[
v(x)=1+\int_x^{\infty_\Gamma}
\frac{e^{2i\Omega(t-x)}-1}{2i\Omega}\,U(t)v(t)\,dt. \tag{11}
\]

The sign shown follows directly from v''+2iΩv'=Uv. This is a coefficient-only construction: compute (g,g_x) at a matching point and evaluate its Wronskian with the horizon vector; no full connection matrix is needed.

A valid contour must stay on the chosen analytic spatial sheet, avoid the singularities of the radial equation and tortoise map, and have a tail on which the exponential is bounded/decaying and \(\int|U\,dx|<\infty\). For fixed Ω away from zero, a tail with Im(Ω e^{iθ})>0 supplies exponential decay. On a compact frequency patch and an admissible tail, the Volterra series converges with standard ordered-integral/factorial bounds; sufficiently far out one can also use a contraction bound. This establishes an analytic nonzero outgoing solution on that specified patch.

At Ω=−iN/2, a contour in the positive-real spatial direction does NOT satisfy that decay condition. A rotated tail must point into a suitable left half-plane sector, with upper/lower deformations giving specified lateral continuations. For Schwarzschild the infinity tail produces a frequency branch cut customarily placed on the negative imaginary axis. Hence write \(g^+\), \(g^-\), and \(\widehat D^+\), \(\widehat D^-\) on that axis. There is no unqualified single outgoing value there without a sheet/continuation prescription. Rotating a contour is not merely numerical stabilization: the permitted deformation and its homotopy encode the continuation branch. Working in r=z with x=z+log(z−1) explicitly tracked is safer than assuming the inverse tortoise map is single valued in an arbitrary complex-x domain.

## Primary literature verified

1. A. Maassen van den Brink, Phys. Rev. D 62, 064009 (2000), DOI 10.1103/PhysRevD.62.064009, arXiv:gr-qc/0001032, especially §II A (Frobenius/anomalous/miraculous distinction and the derivative method for the resonant coefficient), §III (gravitational algebraically special miracle), and §VI (uniqueness of that RWE miracle). https://arxiv.org/pdf/gr-qc/0001032
2. P. T. Leung, A. Maassen van den Brink, K. W. Mak, K. Young, Classical Quantum Gravity 20, L217–L223 (2003), DOI 10.1088/0264-9381/20/16/101, arXiv:gr-qc/0301018. Page 2 explicitly describes f's poles at −in/2, cancellation by rescaling without changing the Green kernel, miraculous exceptions, and the two lateral infinity solutions. https://arxiv.org/pdf/gr-qc/0301018
3. G. B. Cook and M. Zalutskiy, Phys. Rev. D 94, 104074 (2016), arXiv:1607.07406, provides later Kerr examples distinguishing anomalous and miraculous polynomial modes. Abstract verified; full PDF was unavailable through the search tool, so no detailed theorem in the above audit relies on it. https://arxiv.org/abs/1607.07406

The scalar no-miracle proof, polynomial recurrence, and identity (8) above were derived independently in this audit; do not attribute those derivations to a paper unless a matching source is checked.

## Verification

`scalar_resonance_verify.py` checks strict column dominance symbolically, verifies p_N=(−1)^N det A_N for N=1,…,6, and tests (8) at 90-digit precision for N=1,…,6, ell=0,1,2, z=1.3 using 200 series terms. Relative errors are below 10^{-65}. It also checks (6) against H/Γ(c) off resonance at λ=1+10^{-10},1+10^{-20},1+10^{-30}. These checks supplement the proofs; they do not establish any QNM frequency.
