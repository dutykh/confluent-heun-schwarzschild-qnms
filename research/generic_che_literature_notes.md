# Generic CHE certification and resonant normalization: targeted literature audit

Research date: 19 November 2026. This is an internal source audit, not a proof of novelty. Primary papers were checked wherever accessible. Details available only from abstracts are explicitly marked. Earlier connection/QNM sources are recorded separately in `qnm_connection_literature_notes.md`.

## Assessment

Neither outgoing Volterra construction, explicit residual remainder estimates, parameter-uniform asymptotics, nor the failure of a single normalized horizon line at simultaneous resonance/apparency is a new general principle. There are exceptionally close predecessors for each. The strongest defensible project is a practical, explicit five-parameter CHE specialization with evaluated constants, derivative control, a correct parameter-domain statement, and an independently validated implementation. Its novelty still needs comparison at the level of actual formulas and performance, not merely methodology.

## 1. Generic CHE numerical evaluation already exists

O. V. Motygin, **On evaluation of the confluent Heun functions**, arXiv:1804.01007 (2018), [primary PDF](https://arxiv.org/pdf/1804.01007), [author implementation](https://github.com/motygin/confluent_Heun_functions).

The convention is

\[
H''+\left(\frac\gamma z+\frac\delta{z-1}+\epsilon\right)H'
+\frac{\alpha z-q}{z(z-1)}H=0.
\]

This is the generic five-parameter CHE. Sections 3–5 provide origin Frobenius recurrences, exceptional logarithmic cases, infinity formal expansions, and ordinary-point Taylor continuation. Section 7 performs regular-to-regular and regular-to-irregular connection matching, distinguishing sectors at infinity. The infinity initialization uses optimal truncation, with a heuristic radius chosen from the least factorial term. Section 6, equations (21)–(22), supplies residual/last-term error proxies, not rigorous enclosures. Thus a certified replacement for this initialization and continuation is a meaningful algorithmic improvement, although its mathematical ingredients have the predecessors below. The 2015 Motygin paper concerns general Heun; the directly relevant confluent paper is 2018.

## 2. Rank-one irregular residual and derivative bounds: Olver

F. W. J. Olver, **Asymptotic solutions of linear ordinary differential equations at an irregular singularity of rank unity**, *Methods and Applications of Analysis* 4 (1997), 375–403, [publisher full text](https://intlpress.com/site/pub/files/_fulltext/journals/maa/1997/0004/0004/MAA-1997-0004-0004-a002.pdf).

Arbitrary-order equations with coefficients analytic outside a disk and distinct leading characteristic values are treated. Theorem 3.1 constructs normalized sectorial solutions, uniform closed-subsector asymptotics, and differentiated asymptotics through order n−1. Sections 5–9 use the residual of a finite formal truncation to construct a singular Volterra equation, solve it by Liouville–Neumann iteration, and estimate the remainder and its differentiated forms. Equation (5.7) is the key residual integral equation. Generic CHE outside a disk is within this class when the leading exponential rates are distinct. Therefore a CHE-specific residual certificate is a specialization of a well-developed theory.

Related verified primary references: Olver (1965), DOI [10.1137/0702017](https://doi.org/10.1137/0702017), second-order rank-one analysis; Olver–Stenger (1965), DOI [10.1137/0702018](https://doi.org/10.1137/0702018), arbitrary-rank errors. Their publisher metadata was checked, but the readable 1997 full text is the better theorem audit source.

[DLMF §2.7(ii)](https://dlmf.nist.gov/2.7) gives the scalar formal recurrence and rank-one sectorial solutions, equations 2.7.10–2.7.16, and points to Olver's *Asymptotics and Special Functions*, Chapter 7, for errors. Coincident characteristic roots require a different treatment; the ε→0 CHE degeneration must not be included silently in a fixed-rank certificate.

## 3. Uniform holomorphic parameter dependence is also established

G. Cotti, D. Guzzetti, D. Masoero, *Journal of Differential Equations* 428 (2025), 1–58, DOI [10.1016/j.jde.2025.01.085](https://doi.org/10.1016/j.jde.2025.01.085), [arXiv:2310.19739v2](https://arxiv.org/html/2310.19739v2).

For systems Y′=(Λ+R)Y, Theorem 2.21 establishes a jointly holomorphic normalized fundamental matrix with uniform parameter asymptotics under uniform Levinson, integrability, and decay hypotheses. Corollary 2.22 yields relative O(z^−δ) from R=O(z^−1−δ), uniformly in parameters. Theorem 2.27 and Corollary 2.28 address an individual subdominant column on enlarged domains. Equation (2.35) is an infinity-endpoint Volterra operator. This exceeds generic CHE in scope, but is not a ready-to-evaluate CHE numerical certificate.

G. Degano, D. Masoero, **A primer of the complex WKB method, with application to the ODE/IM correspondence**, [arXiv:2501.05957v2](https://arxiv.org/pdf/2501.05957), Lecture I, Proposition 1.4 and Corollary 1.5: general scalar equations with analytic parameters, admissible curve families, factorial iterated-Volterra estimates, convergent Neumann construction, and holomorphic dependence. Do not import constants uncritically: the uniqueness statement of Theorem 1.7 appears to omit an endpoint derivative condition retained in Lemma 1.1; finite ordinary endpoints need care. Independently prove the specialized kernel bounds and uniqueness actually used.

The supplied Sibuya book is the 1975 polynomial-potential monograph. Its §6, Theorem 6.1, printed pp.15–17, gives entire parameter dependence and compact-parameter uniform function/derivative asymptotics for polynomial Schrödinger potentials; §7, Theorem 7.2, gives uniqueness. It is valuable background but not literally a rational-CHE theorem.

## 4. Certified finite continuation and logarithmic Frobenius series

M. Mezzarobba, **Rigorous Multiple-Precision Evaluation of D-Finite Functions in SageMath**, [arXiv:1607.01967](https://arxiv.org/pdf/1607.01967), §§3–4: ball-arithmetic Taylor continuation and connection matrices between ordinary or regular-singular endpoints, using canonical power-logarithm data. Irregular endpoint normalization is outside that implementation's described scope.

M. Mezzarobba, **Truncation Bounds for Differentially Finite Series**, *Annales Henri Lebesgue* 2 (2019), 99–148, DOI [10.5802/ahl.17](https://doi.org/10.5802/ahl.17), [author full text](https://www.marc.mezzarobba.net/papers/Mezzarobba_AdaptiveBounds_2019.pdf).

The algorithm combines residuals with majorants, covers regular-singular generalized powers and logarithms, and bounds derivatives as well as values. Section 1.2 explicitly excludes irregular expansions, even individually convergent solutions there. Algorithms 7.1 and 7.4 rigorously bound rational recurrence sequences using interval arithmetic; their interval polynomial operations are natural tools for parameter boxes, although this paper does not state our desired CHE-specific compact-parameter theorem. Section 8 explains numerical remainder extraction. Therefore certified finite CHE continuation is established in greater generality.

## 5. Elementary CHE integral series

P.-L. Giscard, A. Tamar, **Elementary integral series for Heun functions: Application to black-hole perturbation theory**, *J. Math. Phys.* 63 (2022), 063501, DOI [10.1063/5.0071081](https://doi.org/10.1063/5.0071081), [arXiv:2010.03919v1](https://arxiv.org/pdf/2010.03919).

Corollary 2 supplies generic CHE integral-series solutions based at an ordinary point. For any finite endpoint on a singularity-free real interval, bounded kernels yield convergent Volterra series. This is not, by itself, an outgoing infinity normalization or a certified horizon-to-infinity coefficient. In the checked preprint appendix after equation (28), the displayed tail estimate omits interval-length factors and is labelled relative although displayed as absolute. Do not reuse it as a numerical certificate without correcting it and checking the publisher version.

Independent elementary estimate: if |f|≤M on a path segment of length L, then

\[
|f^{*n}|\le \frac{M^nL^{n-1}}{(n-1)!},\qquad
\left|\sum_{n=m+1}^{\infty}f^{*n}\right|
\le M e^{ML}\frac{(ML)^m}{m!}.
\]

The estimate printed in that preprint fails already for f≡1, L=2, m=1. This observation concerns that displayed bound, not the integral-series convergence theorem.

## 6. Gamma regularization can vanish: explicit prior art

O. V. Motygin, **On regularization of the Heun functions**, [arXiv:2010.09053](https://arxiv.org/pdf/2010.09053), especially §3.

This paper concerns general Heun. It explicitly observes that dividing the origin-normalized solution by Γ(γ) can produce the identically zero solution at γ a nonpositive integer when the resonant recurrence obstruction vanishes. The proposed subtraction in equation (12) mixes the two Frobenius branches to remove the pole locally. The global patch in equation (14) is smooth in parameters rather than holomorphic, owing to cutoff functions. It is not a scalar regularization preserving the original selected Frobenius line. Thus the vanishing gamma-regularized section is already documented; a CHE version is not a discovery of the phenomenon.

## 7. Directional resonance/apparency obstruction: pole-skipping literature

M. Blake, R. A. Davison, D. Vegh, **Horizon constraints on holographic Green's functions**, *JHEP* 01 (2020), 077, [arXiv:1904.12883v3](https://arxiv.org/pdf/1904.12883), DOI [10.1007/JHEP01(2020)077](https://doi.org/10.1007/JHEP01%282020%29077).

Section 2.2, equation (2.17), and §3.2, equation (3.8), explicitly demonstrate that the horizon Frobenius coefficient selected at resonance/apparency depends on the direction of approach in frequency–momentum space. Appendix A identifies the logarithmic obstruction. Equation (3.9) discusses multiple-root cases in which the generic first-order slope picture fails. Thus directional failure of a unique ingoing line has long been part of pole-skipping theory.

Y.-S. Choun, **Fuchsian Resonance and a Horizon-to-Boundary Dictionary for Pole-Skipping**, arXiv:2609.00942v1, 1 September 2026. [Primary abstract](https://arxiv.org/abs/2609.00942); [full author manuscript on mirror](https://arxiv.cnoctave.top/pdf/2609.00942). The official abstract was search-verified; direct official PDF fetch failed, so theorem inspection used the manuscript mirror.

Theorem 4.2 (p.30) proves that a pole-removed ingoing family restricts at resonance to a nonvanishing normalization factor times the Frobenius obstruction times the higher-root solution. Theorem 4.3 (p.31) proves the simultaneous vanishing of both regularized connection coefficients precisely at obstruction zero, assuming an analytic independent boundary basis and a fixed continuation path. Sections 5 and 6.5 show how approach curves select different projective connection data. A blow-up is not explicitly named, but the projective nonextension mechanism is already present. One must not claim it as a new negative theorem merely by expressing it geometrically.

K.-i. Kubota, H. Motohashi, **Pole skipping, avoided crossing, and resonant excitation in Kerr quasinormal modes near algebraically special frequencies**, *Phys. Rev. D* 114 (2026), L081502, published 28 October 2026, DOI [10.1103/29n2-dsmr](https://doi.org/10.1103/29n2-dsmr), [arXiv:2605.17840](https://arxiv.org/pdf/2605.17840). Primary publisher metadata and manuscript checked. This is directly relevant to QNM exceptions: zeros associated with Matsubara frequencies can cancel poles near algebraically special Kerr/Schwarzschild frequencies. A zero of a normalization-dependent determinant must not automatically be reported as a physical Green-function pole.

## 8. Elementary CHE illustration, independently derived

For the convention in §1, H(0)=1 implies H′(0)=−q/γ. Fix α,δ,ε, and approach (γ,q)=(0,0) through q=cγ. The recurrence then has limits H(0)=1 and H′(0)=−c; at the limiting parameters the origin becomes ordinary and these initial data determine distinct solutions. Hence the original selected Frobenius line cannot have a path-independent continuous extension through that two-parameter point. A nonzero jointly holomorphic scalar-renormalized section preserving this line is impossible there. This does not obstruct an analytic basis of the solution bundle at a fixed ordinary point, and a prescribed one-parameter physical path may still select a limit.

The first coefficients of H/Γ(γ) have leading forms γ and −q, suggesting the ideal (γ,q) and its standard blow-up. For higher resonances the natural candidate is (γ+N, obstruction_N). Proving an ideal identity and a lifted line bundle, with multiplicities and nontransverse apparent loci treated explicitly, is more precise than the basic slope argument; this search does not establish its novelty. The 2020 and 2026 sources above must be cited before making such a claim.

## 9. What a useful theorem/implementation would have to deliver

1. Explicit CHE residual coefficients and computable value-plus-derivative error constants, uniform on a stated compact parameter chart.
2. A fixed admissible infinity ray/sector with quantitative exponential separation, ε bounded away from the rank-changing locus, and a fixed singularity-avoiding continuation homotopy.
3. Certified Frobenius and Taylor transport, connected to known D-finite algorithms, followed by an explicit propagated bound for one Wronskian/connection entry.
4. A normalization-domain statement excluding or resolving the gamma-regularized zero locus. On a physical parameter curve, remove only proved common vanishing orders; do not discard poles or zeros by convention.
5. Verified numerical comparison with Leaver/MST or another independent implementation, preferably certified QNM root counts or isolations rather than agreement of floating-point approximations alone.

The recent extremal-Kerr/DCHE paper Chen et al., arXiv:2610.39095v3, Appendix B, already supplies endpoint residual Volterra bounds, derivative bounds, finite Taylor certificates, and an interval winding-count example; see the separate QNM notes. This prevents any claim of the first such certified black-hole pipeline. A generic CHE theorem can be a useful extension if its exact domain, constants, and practical advantage are demonstrated.

## Leads not fully audited

- R. Schäfke (1984), SIAM J. Math. Anal. 15, 253–271, DOI 10.1137/0515021: primary abstract describes regular-to-rank-one-irregular connection coefficient limits and a generalized-Heun application; full quantitative theorem not retrieved.
- Olde Daalhuis–Olver (1994), Proc. R. Soc. A445,39–56, [institutional full text](https://ir.cwi.nl/pub/1486/1486D.pdf): exponentially improved rank-one remainders with Stokes constants as input, not an input-only connection solver.
- Jeong–Pedroni, arXiv:2607.24920v2: resonant gauge-theoretic Floquet Heun solutions and NS blow-up equations; primary abstract checked, full text not retrieved. This is a distinct resonance variable, not automatically the origin Frobenius base ideal.
- Janowicz, arXiv:2609.16969: degenerate double-confluent Heun, indexed primary abstract claims Olver-type infinity solutions with explicit errors; full text inaccessible. No precise theorem claim made.

No exhaustive absence or impossibility claim follows from this audit.
