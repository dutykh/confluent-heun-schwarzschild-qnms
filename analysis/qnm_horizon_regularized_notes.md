# Regularized scalar horizon evaluator: scope, validation, and prior art

`qnm_horizon_regularized.py` evaluates

\[
\widehat u_H=\frac{u_H}{\Gamma(1-2i\Omega)},\qquad \partial_x\widehat u_H,
\]

at the default matching point z=1.5 (or any supplied ordinary point with 0<|1−z|<1). The equation, normalization, and time convention are documented at the top of the executable module. It has only the `mpmath` dependency.

## Usage

```python
from qnm_horizon_regularized import horizon_regularized, evaluate_horizon_regularized
u, ux = horizon_regularized('0.7-0.2j', ell=2, dps=80)
data = evaluate_horizon_regularized('-0.5j', ell=0, dps=80)
print(data.u, data.ux, data.terms, data.ode_residual_scaled)
```

```bash
python3 qnm_horizon_regularized.py --omega-re 0 --omega-im -0.5 --ell 0 --dps 80
python3 qnm_horizon_regularized.py --test --output qnm_horizon_regularized_results.json
```

Decimal string input preserves the supplied precision. A Python float/complex input has already been rounded to binary-double precision. Nearby frequencies are not snapped to the resonance grid.

## Recurrence and exact-grid identity

Set λ=2iΩ, L=ell(ell+1), t=1−z. The implemented polynomial recurrence is

\[
p_0=1,\quad p_{-1}=0,\qquad
p_n=[n(n-1)-\lambda^2-L]p_{n-1}
-\lambda(n-1)(n-1-\lambda)^2p_{n-2},
\]

and

\[
\widehat H_H=\sum_{n\ge0}\frac{p_n t^n}{n!\Gamma(n+1-\lambda)}.
\]

Reciprocal gamma is evaluated directly, so the method does not divide by zero or multiply a divergent computed HeunC by zero. H, H_t, and H_tt are summed together; the returned x derivative follows from the explicit prefactor and chain rule. Series stopping requires 12 consecutive small terms after a frequency-dependent minimum index. A cancellation indicator and the independent differential-equation residual trigger increased working precision if needed. These are numerical diagnostics, **not rigorous interval enclosures or a certified tail bound**.

For Ω_N=−iN/2, N≥1, let Δ_N(L)=det A_N(L), where

\[
(A_N)_{kk}=N^2+L-k(k+1),\quad
(A_N)_{k,k+1}=-(k+1)(N-k-1),\quad
(A_N)_{k,k-1}=-N(N-k).
\]

All indices run from 0 to N−1. The strict column-dominance gap is L+N−k>0. Thus Δ_N(L)>0, the normalized horizon pole is genuine, and the regularized vector is nonzero. With this module's normalization,

\[
\widehat u_H(\Omega_N,z)=\frac{e^{-N}\Delta_N(L)}{N!}u_H(-\Omega_N,z).
\]

The module exposes `resonance_determinant(N, ell)`, evaluated by an exact integer tridiagonal determinant recurrence.

## Tests actually run

The JSON output records all cases. At requested 80 digits, with 110 working digits:

- 36 exceptional-frequency cases: N=1,2,3,4,8,12,16,24,32 and ell=0,1,2,5. Direct regularized series compared to the ordinary upper-half-plane recurrence using the exact identity. Maximum relative/absolute-scaled discrepancy: 1.951×10^−93. Maximum normalized ODE residual: 2.618×10^−92.
- Eight generic complex-frequency cases compared to the ordinary recurrence times reciprocal gamma.
- Three near-resonant cases, N=1,4,12, displaced by a real 10^−30, checked against the independently computed unregularized recurrence before gamma scaling.
- Six static cases, ell=0,…,5, checked against the exact Legendre solution u=z P_ell(2z−1) and its x derivative.

No test here asserts a QNM frequency, constructs the infinity solution, selects a frequency sheet there, or certifies a global spectral count.

## Targeted novelty assessment

The main mechanisms are established prior art and should be described as independent derivations and implementation checks.

1. **Scalar poles at every negative half-integer imaginary frequency, with no miracle:** explicitly stated for the general-spin Schwarzschild problem in Casals and Ottewill, *Analytic Investigation of the Branch Cut of the Green Function in Schwarzschild Space-time*, Phys. Rev. D **87**, 064010 (2013), arXiv:1210.0519, §§IV A–V. The sole stated exception is the gravitational s=2 algebraically special frequency; the scalar s=0 case has no such exception. https://arxiv.org/pdf/1210.0519
2. **Regularized horizon solution proportional to the opposite-frequency solution:** explicitly stated in the same paper, §V, immediately after Eq. (5.1). Their regularizing factor is a sine; the present gamma factor and e^−N Δ_N/N! give a specific normalization of the already known mechanism.
3. **Gamma regularization, collapse to the larger-exponent solution, and accidental identically zero limits:** explicitly discussed in Motygin, *On regularization of the Heun functions*, arXiv:2010.09053, §3, PDF p. 4. That paper studies general Heun and explicitly traces the confluent-Heun regularization to Motygin's earlier 2018 paper, arXiv:1804.01007. https://arxiv.org/pdf/2010.09053
4. **Primitive endpoint vectors before determinant formation:** explicitly required in Chen et al., *Complete quasinormal modes of Type-D black holes. II. Extremal Kerr black holes*, arXiv:2610.39095v3, Appendix B.1. The extremal endpoint construction differs from the present nonextremal Frobenius series, but the normalization principle is already explicit. https://arxiv.org/pdf/2610.39095

The exact column-dominance proof, the displayed p_n scalar recurrence, and the explicit Δ_N normalization were not located in the targeted search. This is **not evidence of novelty**: the recurrence is an algebraic rescaling of a standard recurrence, and the other two formulas provide an independent proof and explicit normalization of known properties. No novelty claim should be based on them without a substantially deeper historical review.

Suggested phrasing: “We give a direct matrix proof of the known absence of miraculous horizon cancellations in the scalar Schwarzschild problem and implement a division-free evaluator with explicit exceptional-frequency normalization.”
