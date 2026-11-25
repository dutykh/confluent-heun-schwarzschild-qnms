# Rigorous arithmetic and scalar QNM prior: 21 November 2026

This is a targeted primary-source audit, not an exhaustive novelty determination.

## Validated arithmetic

1. Fredrik Johansson, *Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic*, IEEE Transactions on Computers **66**(8), 1281–1292 (2017), DOI [10.1109/TC.2017.2690633](https://doi.org/10.1109/TC.2017.2690633), [author manuscript](https://arxiv.org/pdf/1611.02831). Section 2 describes `arb_t` real balls and `acb_t` complex rectangles; Section 5 states enclosure semantics. A mathematical error bound can be added to the enclosure, but arithmetic alone does not prove the analytic approximation theorem or the existence of a zero. Arb is now incorporated into FLINT. [Current official acb documentation](https://flintlib.org/doc/acb/index.html).

2. Nathalie Revol and Fabrice Rouillier, *Motivations for an Arbitrary Precision Interval Arithmetic and the MPFI Library*, Reliable Computing **11** (2005), 275–290, DOI [10.1007/s11155-005-6891-y](https://doi.org/10.1007/s11155-005-6891-y), [author PDF](https://perso.ens-lyon.fr/nathalie.revol/publis/Revol-Rouillier2005_Article_MotivationsForAnArbitraryPreci.pdf). MPFI provides real interval arithmetic using MPFR endpoints and directed rounding. Complex rectangles may be implemented by pairs of real intervals; do not describe MPFI as inherently the complex ball library Arb/acb. MPC is correctly rounded complex floating-point arithmetic, not interval arithmetic.

3. [Official mpmath interval context documentation](https://mpmath.readthedocs.io/en/stable/contexts.html) describes `iv.mpc` rectangular complex intervals, but calls interval support experimental and says many functions do not properly support intervals. Crucially, [official root-finding documentation](https://mpmath.org/doc/current/calculus/optimization.html), immediately after the `findroot` description, explicitly says an interval returned by `iv.findroot` is **not guaranteed to contain a zero**. It verifies a small residual, not an existence/uniqueness enclosure. A separate interval-Newton/Rouché/argument-principle proof is required. Exact integers/rationals or decimal strings should initialize enclosures; arbitrary Python binary floats are different exact inputs.

## Closest recent direct scalar amplitude computation

Danilo S. Ferreira and Luis B. Castro, *Scalar quasinormal and quasibound states of weakly magnetized Ernst–Schwarzschild spacetime: a confluent-Heun Wronskian approach*, arXiv:2609.17727v1, 15 September 2026. [Canonical identifier](https://arxiv.org/abs/2609.17727); [full primary manuscript copy read](https://borowiecki.dev/pdf/2609.17727). Official arXiv fetch was disabled in this session, but title/authors/date were independently indexed by SciRate and alphaXiv and the full PDF carries the arXiv identifier.

The paper treats the separated weak-field Ernst scalar equation, equivalent to a massive scalar Schwarzschild radial equation. The massless Schwarzschild case is included when the effective magnetic mass vanishes. Its Section 3 gives the incoming amplitude and spectral Wronskian through a large-order Jaffé coefficient limit (Eqs. 39–41), explicitly citing Schäfke (1984) and Lutz–Schäfke (1997). For their normalization,

\[
S_k=e^{-2\sqrt{-\alpha}\sqrt{k}}k^{3/4-\chi/2}d_k,
\qquad
W=-4h\sqrt\pi e^{\alpha/2}(-\alpha)^{\chi/2-1/4}\lim_{k\to\infty}S_k.
\]

Branches are fixed to the selected Stokes sector. Equation 42 identifies the QNM condition with the zero incoming coefficient on the radiative sheet. Threshold alpha=2h=0 and resonant parameters are outside the stated generic construction or require a limiting prescription. Section 4 approximates the limit with N=500 and checks stability with increasing N and arithmetic precision; Section 5 uses N=700 for quasibound states. Schwarzschild comparison is explicit. No finite-N rigorous tail bound or interval root certificate was identified in the full nine-page manuscript. Therefore a single coefficient without outgoing back-propagation plus Schwarzschild validation is established prior. A quantitatively proved forward-state residual certificate is a narrower possible distinction.

## Schäfke theorem access limitation

The full 1984 and 1997 original articles remained unavailable through accessible endpoints. Publisher abstracts support the broad scope only:

- Reinhard Schäfke, *A Connection Problem for a Regular and an Irregular Singular Point*, SIAM Journal on Mathematical Analysis **15** (1984), 253–271, DOI [10.1137/0515021](https://doi.org/10.1137/0515021).
- Donald A. Lutz and Reinhard Schäfke, *Calculating connection coefficients for meromorphic differential equations*, Complex Variables **34** (1997), 145–170, DOI [10.1080/17476939708815044](https://doi.org/10.1080/17476939708815044).

An accessible explicit reproduction is G. V. Kraniotis, *The massive Dirac equation in the Kerr–Newman–de Sitter and Kerr–Newman black hole spacetimes*, [arXiv:1801.03157v5](https://arxiv.org/pdf/1801.03157), Theorem 2 / Eq. 167, printed page 31. It attributes the generalized-Heun connection-coefficient limit to Schäfke (1984), imposes mu_1 not a positive integer, alpha outside the nonnegative real axis, and Re(a)<1/2, and states a k^(-1/2) asymptotic expansion of the scaled sequence. This is not itself a verified explicit finite-N remainder certificate. Do not silently replace the original theorem's hypotheses with a generic-sector assertion.

## What the current search does and does not establish

The strongest previously verified full interval QNM certificate located is Chen–Zhang–Cao–Jing, *Complete quasinormal modes of Type-D black holes. II. Extremal Kerr black holes*, [arXiv:2610.39095v3](https://arxiv.org/abs/2610.39095), Appendix B: explicit asymptotic Volterra and finite-continuation error bounds, interval arithmetic, and argument-principle count for a specified extremal Kerr gravitational l=m=3 rectangle. This is not a scalar Schwarzschild result, and it does not certify the whole numerical survey. The paper was fully inspected in the preceding audit; new direct fetches were disabled during this update.

Targeted searches for scalar Schwarzschild QNM interval/certified/computer-assisted results did not locate a verified direct prior root enclosure. This limited search does **not** establish a first-ever claim. The current computation gives reproducible rigorous enclosures for nine specified scalar modes and a fixed-frequency comparison of the actual errors before and after the residual correction.
