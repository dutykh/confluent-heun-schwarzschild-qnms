# Generic CHE resonance: a holomorphic frame exists, but the distinguished Frobenius line need not extend

## Main conclusion and limits of novelty

There are two different questions.

1. **Can one construct a nondegenerate holomorphic two-solution frame through a simultaneous resonance and vanishing Frobenius obstruction? Yes, locally on a fixed punctured spatial domain.** An explicit formula is given below. It necessarily mixes the two off-resonance Frobenius solutions and is closely related to Motygin's local subtraction regularization.
2. **Can one preserve the canonically normalized exponent-zero Frobenius solution line and extend it holomorphically, without extra data, through every such parameter intersection? No.** Already N=1 gives a rigorous counterexample with every projective solution direction as a possible limiting line. A blow-up records these directions, but is not needed merely to construct an arbitrary holomorphic full frame.

The no-go phenomenon is not new to Heun theory: the classical identities M(a,a,z)=exp(z) and M(0,b,z)=1 show the same path dependence at the Kummer parameter point (a,b)=(0,0). Thus these results are useful mathematical structure and safeguards for a connection algorithm, not by themselves a defensible breakthrough or established novelty claim.

## 1. Equation and recurrence

Use the Maple canonical CHE

\[
H''+\left(\alpha+\frac{\beta+1}{z}+\frac{\gamma+1}{z-1}\right)H'
+\frac{Az+B}{z(z-1)}H=0,
\]

where

\[
A=\delta+\frac{\alpha(\beta+\gamma+2)}2,
\qquad
B=\eta+\frac\beta2+\frac{(\gamma-\alpha)(\beta+1)}2.
\]

For the normalized exponent-zero solution H_0=Σ h_n z^n, h_0=1, h_{−1}=0,

\[
n(n+\beta)h_n=a_nh_{n-1}+b_nh_{n-2},\tag{1}
\]
\[
a_n=(n-1)(n+\beta+\gamma-\alpha)+B,
\qquad b_n=\alpha(n-2)+A.
\]

Define p_n=n!(β+1)_n h_n. Then

\[
p_0=1,\quad p_{-1}=0,\qquad
p_n=a_np_{n-1}+(n-1)(n-1+\beta)b_np_{n-2}.\tag{2}
\]

The gamma-regularized solution is

\[
\widehat H=\frac{H_0}{\Gamma(\beta+1)}
=\sum_{n\ge0}\frac{p_n}{n!\Gamma(n+\beta+1)}z^n.\tag{3}
\]

All coefficients in (3) are holomorphic in the five CHE parameters. The series is normally convergent on compact parameter sets and |z|≤r<1. One proof is to use (2) for the finite potentially resonant indices, and then use (1) for the regularized coefficients for sufficiently large n: uniformly on a compact parameter set, the two recurrence multipliers are bounded by 1+C/n and C/n. A polynomial majorant n^K with K sufficiently large gives convergence on each r<1, including z derivatives. Therefore parameter differentiation of (3) is legitimate locally.

Throughout the frame discussion, fix a simply connected spatial domain inside 0<|z|<1 and a branch of log z. No assertion is made that both frame elements are holomorphic through z=0 for noninteger β.

## 2. Resonant obstruction and collapse of the regularized solution

Fix N≥1 and put

\[
\varepsilon=\beta+N,\qquad
\kappa=\frac{p_N(\alpha,\beta,\gamma,\delta,\eta)}{N!}.\tag{4}
\]

Only the restriction of κ to ε=0 is intrinsic to the resonance; (4) chooses one convenient analytic extension off the resonance divisor. Other extensions differing by ε times a holomorphic function change the complementary basis but not the obstruction locus.

Let

\[
S(z)=z^{-\beta}\operatorname{HeunC}(\alpha,-\beta,\gamma,\delta,\eta;z).\tag{5}
\]

Near β=−N, S is holomorphic in all parameters on the fixed spatial branch because its transformed exponent parameter −β is near the positive integer N.

At β=−N, the first N coefficients (indices 0,…,N−1) of (3) vanish. The coefficient at z^N is p_N/N!=κ. The larger-exponent Frobenius solution S starts with z^N and is unique with that leading coefficient. Hence the exact identity on the resonance divisor is

\[
\widehat H\big|_{\varepsilon=0}
=\kappa\big|_{\varepsilon=0}\,S\big|_{\varepsilon=0}.\tag{6}
\]

If κ=0 as well, the entire regularized solution is zero. Directly from (2), p_N=0 implies p_{N+1}=0 at β=−N because the second recurrence multiplier then vanishes, and all later p_n vanish by induction.

The restriction p_N|_{β=−N} is a monic polynomial of degree N in η, since B is affine in η and the first term in (2) raises its degree by one. Its smooth zeros are the ordinary points of the miraculous/apparent locus. Multiple zeros require separate local geometry; the N=1 counterexample below has no such complication.

## 3. A uniform holomorphic full frame: explicit construction

Because the numerator vanishes identically on ε=0, the function

\[
\boxed{F=\frac{\widehat H-\kappa S}{\varepsilon}}\tag{7}
\]

has a holomorphic removable extension across ε=0. This is analytic division by a coordinate, not a numerical subtraction prescription. For ε≠0 it is a linear combination of two solutions of the same CHE; by holomorphic extension it also solves the equation at ε=0.

Use W(f,g)=fg'−f'g, and set

\[
E(z)=e^{-\alpha z}z^{-\beta-1}(1-z)^{-\gamma-1}.
\]

From the leading Frobenius coefficients and Abel's identity,

\[
W(H_0,S)=-\beta E(z),\qquad
W(F,S)=\frac{-\beta}{(\beta+N)\Gamma(\beta+1)}E(z).\tag{8}
\]

Since 1/Γ(1−N+ε)=(-1)^{N−1}(N−1)!ε+O(ε²),

\[
\boxed{W(F,S)\big|_{\beta=-N}=(-1)^{N-1}N!\,
 e^{-\alpha z}z^{N-1}(1-z)^{-\gamma-1}\ne0.}\tag{9}
\]

Thus (F,S) is a nondegenerate holomorphic frame near every point of β=−N, including κ=0. This construction works even at multiple zeros of the obstruction polynomial. It does not require a blow-up.

What it does **not** do is preserve the selected H_0 line. The exact relation is

\[
\widehat H=\varepsilon F+\kappa S.\tag{10}
\]

That distinction prevents a mistaken claim that arbitrary basis regularization automatically fixes a physical endpoint boundary condition.

## 4. No-go theorem for the distinguished line at a smooth simultaneous zero

Assume ε and κ are independent local parameter coordinates at a point P where ε=κ=0, and hold the other three coordinates fixed. The frame (F,S) is holomorphic and nondegenerate there. On the punctured parameter neighborhood where the distinguished solution is nonzero, its projective coordinates in this frame are exactly

\[
[\widehat H]=[\varepsilon:\kappa].\tag{11}
\]

Along κ=cε with ε→0,

\[
\frac{\widehat H}{\varepsilon}\longrightarrow F_P+cS_P.\tag{12}
\]

Since F_P and S_P are independent, different c give different solution lines. Along ε=0, κ→0 the limiting line is [S_P], supplying the point at infinity. Consequently the distinguished line has no continuous, and hence no holomorphic, extension to P. In particular no meromorphic scalar rescaling of H_0 or \(\widehat H\) can produce a nonzero holomorphic section through P while preserving its off-locus line.

This is a failure of a selected line, **not** a failure of holomorphic dependence of the differential equation's full solution space.

## 5. Concrete N=1 counterexample with exact initial jets

Let

\[
\varepsilon=\beta+1,\qquad \kappa=B,
\]

and hold α,γ,δ fixed. The transformation between (β,η) and (ε,κ) is nonsingular. Equation (3) gives the exact jet

\[
\boxed{\big(\widehat H(0),\widehat H'(0)\big)
=\frac{(\varepsilon,\kappa)}{\Gamma(1+\varepsilon)}.}\tag{13}
\]

At ε=κ=0 the singularity at zero is removable. The limiting equation is

\[
y''+\left(\alpha+\frac{\gamma+1}{z-1}\right)y'
+\frac{\delta+\alpha(\gamma+1)/2}{z-1}y=0.\tag{14}
\]

Let U,V solve it with initial jets (1,0) and (0,1), respectively. The derivatives ∂_ε\(\widehat H\) and ∂_κ\(\widehat H\) at P equal U and V: parameter differentiation of L\(\widehat H\)=0 produces no inhomogeneous term because \(\widehat H(P)\equiv0\). Therefore

\[
H_0\big|_{\kappa=c\varepsilon}\longrightarrow U+cV.\tag{15}
\]

Their Wronskian is e^{-αz}(1−z)^{-γ−1}, so the limiting lines are distinct at any ordinary matching point, not only as formal series.

An explicit genuinely confluent family (α≠0, and noninteger γ) is

\[
\alpha=1,\quad\gamma=\tfrac12,\quad\delta=-\tfrac34,\quad
\beta=-1+\varepsilon,\quad
\eta=\tfrac12+\kappa-\tfrac14\varepsilon.
\]

Then A=ε/2 and B=κ. At P,

\[
U(z)=1,\qquad V(z)=\int_0^z e^{-t}(1-t)^{-3/2}\,dt.
\]

The paths κ=0 and κ=ε yield U and U+V. Both have value one at the origin and solve the same limiting equation, but they are not proportional. This is a concrete counterexample to a direction-independent extension of the canonical normalized solution.

## 6. Monodromy normal form

Let m=exp(−2πiβ)=exp(−2πiε). Around zero, H_0 is single valued off resonance and S has multiplier m. In the frame (F,S), with column-vector coefficients, the exact monodromy matrix is

\[
M_0=
\begin{pmatrix}
1&0\\
\kappa\,\dfrac{1-m}{\varepsilon}&m
\end{pmatrix}.\tag{16}
\]

Its apparent quotient is removable, and at ε=0,

\[
M_0=\begin{pmatrix}1&0\\2\pi i\kappa&1\end{pmatrix}.\tag{17}
\]

For κ≠0 it is a nontrivial unipotent matrix whose eigenline is S. At κ=0 it is the identity and every solution line is an eigenline. The off-resonance eigenline with eigenvalue one is (ε,κ), explaining its direction-dependent limit in elementary matrix terms.

## 7. Blow-up and connection determinants

At a smooth point of the codimension-two base locus, blowing up the ideal (ε,κ) resolves the map [ε:κ]. In the chart κ=cε, the pulled-back selected solution divided by ε is F+cS and remains nonzero at ε=0. In the other chart ε=dκ, division by κ gives dF+S. The exceptional projective line records the missing limiting boundary datum. These are local generators of a line bundle; a single global nonvanishing generator over the whole exceptional projective line should not be presumed.

For any other analytically normalized nonzero endpoint solution G, write

\[
\widehat D=W(\widehat H,G)=\varepsilon A_*+\kappa B_*,\quad
A_*=W(F,G),\quad B_*=W(S,G).\tag{18}
\]

At P the scalar \(\widehat D\) automatically vanishes together with the selected horizon vector. There is no single holomorphic common scalar factor to cancel in the original multi-parameter neighborhood. In the blow-up chart its meaningful transform is

\[
\widehat D/\varepsilon=A_*+cB_*.
\]

Thus the limiting spectral condition depends on the direction c. This is why multivariate normalization-zero removal requires a boundary prescription, rather than blindly dividing an isolated scalar zero.

For a prescribed analytic one-parameter physical family through P, the situation is simpler: pull back (ε,κ) to that curve, remove its common lowest power, and retain the resulting limiting direction. If ε=e_1t+… and κ=k_1t+… with e_1≠0, the limit is [F_P+(k_1/e_1)S_P]. One-variable removable-singularity reasoning must not be promoted to a direction-independent statement in five-parameter space.

## 8. Prior work and restrained assessment

- **Motygin 2020:** *On regularization of the Heun functions*, arXiv:2010.09053, §3. It explicitly describes H/Γ, collapse onto the other Frobenius solution, and identically zero regularized limits at vanishing obstruction. Its Eq. (12) subtracts the singular multiple of the second solution locally; the later cutoff is introduced for a globally C∞ definition. Thus local holomorphic mixing is already present. The present CHE frame (7) is an explicit gamma-scaled implementation of the same local mechanism, not a fundamentally new idea. https://arxiv.org/pdf/2010.09053
- **Motygin 2018:** *On evaluation of the confluent Heun functions*, arXiv:1804.01007, is cited by the 2020 paper as the earlier confluent-Heun regularization source. Its PDF was not retrieved in this targeted check, so detailed claims here rely on the verified 2020 discussion. https://arxiv.org/abs/1804.01007
- **Classical Kummer counterexample:** DLMF 13.6.1 and 13.6.3 give M(a,a,z)=exp(z) and M(0,b,z)=1. Both paths tend to (a,b)=(0,0), where the limiting equation is y''−y'=0 away from zero. They produce distinct normalized solution lines. DLMF 13.2.3–13.2.4 also defines the entire gamma-regularized Kummer function. Thus multi-parameter path dependence is already implicit in elementary classical formulas. https://dlmf.nist.gov/13.6 and https://dlmf.nist.gov/13.2
- **Blow-up:** resolution of the projective map [ε:κ] by blowing up its base ideal is standard geometry, not a new special-function operation. See the Stacks Project, §31.33, Definition 31.33.1 and affine charts in Lemma 31.33.2. https://stacks.math.columbia.edu/tag/01OF

The targeted search did not locate this exact CHE parameter-coordinate theorem stated in the same form. That does not establish originality. The most defensible contribution would be a careful, explicit synthesis applied to connection coefficients, with an implemented and validated treatment of a genuinely difficult application. The proposed universal canonical nonzero holomorphic boundary solution is mathematically impossible without retaining extra direction/path data; a universal holomorphic **full frame** is locally possible and already aligned with known regularization methods.
