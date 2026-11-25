# Independent scalar Schwarzschild Chebyshev benchmark

All frequencies in this note are **Omega = 2 M omega**, with time factor exp(-i omega t), a massless scalar field, and L=ell(ell+1). No Leaver frequencies or externally tabulated values were supplied to the spectral scripts.

## Exact reduction

Put z=r/(2M), f=1-1/z, x_*=z+log(z-1). The scalar radial function u satisfies

u_(x_*x_*) + [Omega^2 - f(L/z^2+1/z^3)] u = 0.

Use y=(z-1)/z in [0,1], h=y(1-y)^2, so d/dx_*=h d/dy. Factor

u=e^{i Omega z}(z-1)^{-i Omega}z^{2 i Omega}v(y)
 =e^{i Omega/(1-y)} y^{-i Omega}(1-y)^{-i Omega}v(y).

This gives the ingoing horizon factor and the outgoing infinity formal factor. Set lambda=i Omega, k=1-2(1-y)^2. The exact reduced ODE is

h v'' + [h'+2 lambda k]v'
 + [4 lambda(1-y)-4 lambda^2(2-y)-L-1+y]v = 0,

with h'=(1-y)(1-3y). The two limiting endpoint equations used by collocation are

(1-2lambda) v'(0) + (4lambda-8lambda^2-L-1)v(0)=0,
2lambda v'(1) - (4lambda^2+L)v(1)=0.

These are limiting ODE rows; they are not independent rigorous definitions of the continued scattering boundary conditions. The derivation and coefficient representation are checked symbolically in `qnm_spectral.py`.

## Collocation implementation and conditioning

The preliminary nodal formulation used N+1 Lobatto nodes, y_j=(1-cos(pi*j/N))/2, and the quadratic pencil

A0 + lambda A1 - lambda^2 A2,
A0=diag(h) D^2+diag(h')D-diag(L+1-y),
A1=2diag(k)D+4diag(1-y), A2=4diag(2-y).

A generalized companion eigensolve gave reasonable low roots at small N, but at larger N the nodal D^2 matrices became badly conditioned. In double precision, increasing N from 80 to120 could spoil even the ell=0 fundamental; small residuals did not prevent false accuracy. This preliminary implementation is preserved as `qnm_spectral_nodal.py` and its results as `qnm_spectral_nodal_results.json`, not used for final frequencies.

The final implementation uses **the identical collocation equations in Chebyshev coefficient space**. Put x=2y-1, v=sum_{j=0}^N c_j T_j(x). The differential pencil is

Q0=[(1-x-x^2+x^3)/2] D_x^2 + [(-1-2x+3x^2)/2]D_x + x/2-L-1/2,
Q1=2(1+2x-x^2)D_x+2(1-x),
Q2=2(3-x),

and Q0+lambda Q1-lambda^2 Q2=0. The residual polynomial has degree at most N+1. On the Lobatto nodes T_(N+1)=T_(N-1), so adding the N+1 residual coefficient to the N-1 row gives the exact equivalent collocation pencil. This is not a different physical boundary condition.

The coefficient matrices are assembled using integer Chebyshev product and differentiation identities. Every tested entry is exactly an integer divided by32, with no irrational-grid or floating differentiation error. They are upper Hessenberg. The generalized double-precision eigensolve is used only to generate initial seeds at N40. The finite-matrix roots are then refined at **50 decimal working digits** by differentiated Hessenberg elimination, evaluating det'(Q)/det(Q) for Newton steps. No-pivot elimination requires nonzero intermediate pivots; it succeeded for the reported cases, but is not claimed universally robust. The script stops if Newton fails.

The following checks distinguish discretization error from arithmetic error:

- Exact symbolic operator identities and symbolic action on six Chebyshev basis polynomials.
- Grids N=40,60,80,100,160,240,320 for ell=0,1,2 and three low-damping positive-frequency branches.
- Additional ell=0 n1,n2 grids N=400,480, and n2 N=640.
- Fixed-grid repeats at 70 versus50 working digits, ell0 n2 N160 and ell2 n2 N100. Differences are about3.3e-45 and3.5e-45, dominated by the stored45-digit output strings; hence arithmetic is not responsible for the much larger grid drift.
- Finite pencil residuals, Chebyshev coefficient tails, endpoint residuals, and an exact continuous polynomial-residual bound are provided for N320. These are **not** frequency error estimates.

## Filtering and labels

Initial candidates are restricted to 0<Re Omega<2 and -2<Im Omega<0. The first three positive-frequency branches at N40 are continued through the grids using their previous-grid frequencies as Newton seeds. The labels n=0,1,2 are provisional damping-order labels for those continued branches, not a completeness theorem. New grid-dependent roots are not automatically counted as physical overtones. In particular, a new spurious-looking double-precision root can move ahead of an existing branch in damping order at large N; sorting afresh at every grid is unsafe. An independent scattering/Leaver comparison is needed before final physical classification.

## Essential boundary-condition limitation

For Im Omega<0, an incoming solution at infinity divided by the outgoing prefactor has the asymptotic form

exp[-2i Omega/(1-y)] (1-y)^{2i Omega} [1+...].

Its exponential magnitude decays like exp[-2(-Im Omega)/(1-y)]. It is flat at y=1 along the real interval: adding it can preserve all algebraic endpoint asymptotics and C-infinity regularity. Also, the outgoing formal inverse-radius series generally diverges; claiming ordinary analyticity at the compactified endpoint would not fix this issue correctly.

Therefore endpoint boundedness, endpoint smoothness, a small ODE residual, or a rapidly decaying Chebyshev tail alone does not rigorously exclude incoming contamination. The continued outgoing solution is defined by analytic continuation from its defining domain or by an equivalent sectorial prescription. This spectral calculation supplies independent numerical cross-validation, not a proof of that continuation or of global spectral completeness.

The phenomenon is visible numerically: ell0 n2 at N320 has a continuous polynomial ODE residual bound about2e-44, yet its frequency still moves by roughly3.3e-7 at N400. A tiny residual is not a tiny QNM frequency error in this non-self-adjoint problem.

## Reproduction

From the project directory, start a new working directory, regenerate the full numerical campaign, then compare with the archived records:

```sh
che_replay=$(mktemp -d)
OPENBLAS_NUM_THREADS=1 python -B analysis/qnm_spectral_release.py --recompute "$che_replay/spectral.json"
python -B analysis/qnm_spectral_release.py --check-result "$che_replay/spectral.json"
python -B analysis/qnm_compare.py --spectral "$che_replay/spectral.json" --output-dir "$che_replay"
cmp "$che_replay/qnm_comparison_results.json" analysis/qnm_comparison_results.json
```

The generator uses NumPy, SciPy, SymPy and mpmath. For the two 70-digit repeats the Newton step stopping threshold is 10^-50; the 50-digit grid study uses 10^-30. This distinction reproduces the archived extended-precision final Newton steps. It reconstructs the 68 roots through degree 640, the two 70-digit working-precision repeats, and all nine degree-320 residual checks. With no grid or precision overrides, `qnm_spectral.py --output PATH` runs this same reconstruction. To perform a custom grid study, pass explicit resolution or mode arguments to `qnm_spectral.py`; this defaults to a separate base-grid JSON file. The archived supplementary grid records and the integrated final campaign furnish the comparisons in `--check-result`.

The exact continuous residual bound used in the diagnostics is worth recording. When max|c_j|=1, the unaliased residual's top coefficient is

r_(N+1) = ((N+1)/2-lambda)^2 c_N.

For an exact finite collocation root its entire residual is

r_(N+1)[T_(N+1)(x)-T_(N-1)(x)],

hence its absolute supremum on [-1,1] is at most2|r_(N+1)|. Small finite-matrix residuals are added to that bound when using rounded frequency strings. This validates the residual computation; it does not provide a resolvent estimate or a frequency enclosure.

## Final tracked frequencies and grid differences

The long decimal strings in JSON identify finite-matrix roots; the last-grid differences determine what can responsibly be claimed about convergence.

| ell | n | final N | Omega (selected digits) | absolute change from previous grid |
|---:|---:|---:|---|---:|
| 0 | 0 | 320 | 0.22090987816083937 - 0.20979143417376192 i | 5.0668232007073385393e-21 |
| 0 | 1 | 480 | 0.17223383667284729 - 0.69610489361286667 i | 2.0966810803608126368e-12 |
| 0 | 2 | 640 | 0.15148387116784554 - 1.2021571801069683 i | 4.6164424115582524844e-9 |
| 1 | 0 | 320 | 0.58587226653456542 - 0.19531997782715644 i | 1.265923711494910739e-38 |
| 1 | 1 | 320 | 0.52889730120966508 - 0.61251478311809425 i | 1.4839839347117451693e-24 |
| 1 | 2 | 320 | 0.45907866986260334 - 1.0802668500382144 i | 6.2844995365345246998e-16 |
| 2 | 0 | 320 | 0.96728774442142597 - 0.19351755195657573 i | 0.0 |
| 2 | 1 | 320 | 0.92770115803953113 - 0.59120787397592505 i | 5.0198226187185118603e-37 |
| 2 | 2 | 320 | 0.86108810875315154 - 1.017116804308549 i | 1.1957712655899063906e-27 |
