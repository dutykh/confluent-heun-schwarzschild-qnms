# Forward-only incoming-coefficient calculation

Files: `qnm_forward_coefficient.py`, `qnm_forward_coefficient_results.json`, and `qnm_forward_coefficient_comparisons.py`. The first script implements and runs the forward calculations; the second script adds same-endpoint bound and backward-transport controls. Earlier deliverables are unchanged.

## Exact construction

Use Ω=2Mω, exp(−iωt), z=r/(2M), x=z+log(z−1), f=(z−1)/z, and V(z)=ℓ(ℓ+1)/z²+1/z³. Define h=exp(iΩx)u_H for the unit ingoing horizon solution. Then

 h_xx−2iΩh_x−fVh=0.

The pair (h,h_x) therefore obeys the existing rational Taylor system with Ω replaced by−Ω. Initialize at z=1.5 from the convergent horizon series, propagate along the real interval to the chosen anchor c if necessary, and then FORWARD along z=c+is to z_R=c+iR. No outgoing backward ODE is solved.

Let P_N be the outgoing finite inverse-power polynomial and r_N=fP_N''+(f'+2iΩ)P_N'−VP_N its two-term exact residual. Define

 A_N=[(2iΩh−h_x)P_N+h fP_N']/(2iΩ) at z_R.

Direct differentiation, using the h-system, gives

 d(2iΩ A_N(z))/dz=h(z)r_N(z).

On the sector and ray considered here, h is bounded, the forward Volterra equation implies a finite limit h(∞), and h_x→0. That limit is the incoming amplitude in the chosen continued scattering normalization. Consequently

 A_in−A_N=∫_(z_R)^(∞ on ray) h(z)r_N(z)dz/(2iΩ).

The ray, branch, and frequency-domain qualifications from the contour construction remain essential. This is not asserted to be the first use of forward Jost functions or finite-radius matching.

## Conditional tail estimates

Let F bound∫|r_N|ds, V_R bound∫|V|ds, and k bound the forward kernel on the tail. For exact terminal h_R,hx_R,

 sup_tail|h|≤(|h_R|+k|hx_R|)exp(kV_R),
 E_direct=F exp(kV_R)(|h_R|+k|hx_R|)/(2|Ω|).

These formulas are exact-arithmetic statements. In this implementation h_R,hx_R are numerically propagated, so the displayed tail bounds do not include their arithmetic/Taylor errors and are not full certificates.

At the SAME z_R and N, the conventional outgoing-state error estimates E_m=kFexp(kV_R), E_mx=Fexp(kV_R) imply

 E_state=Fexp(kV_R)[|h_R|+k|2iΩh_R−hx_R|]/(2|Ω|).

Neither bound uniformly dominates. Near a QNM, hx_R≈2iΩh_R, so the conventional terminal-state projection can be sharper. Off-root, after the decaying component becomes small, hx_R≈0 and the direct bound can be sharper. Matching-point relocation must not be confused with a new uniformly sharper theorem.

## A useful weighted refinement

Write

 a=h_R−hx_R/(2iΩ), b=hx_R/(2iΩ),
 h_free(s)=a+b exp(2iΩ[x(s)−x(R)]), Q=|a|+|b|.

The exact forward Volterra series yields

 |h(s)−h_free(s)|≤Q[exp(k V(R,s))−1]≤Q expm1(kV_R).

If Fw bounds∫exp(−μ(s−R))|r_N(s)|ds, where the ray gives |exp(2iΩΔx)|≤exp(−μ(s−R)), then

 E_weighted=[|a|F+|b|Fw+Q expm1(kV_R)F]/(2|Ω|).

This uses the actual a,b without assuming that either vanishes. The existing positive weighted-power bounds evaluate Fw. Take the minimum of independently valid bounds, including the free-component sup bound QFexp(kV_R)/(2|Ω|).

In identical-endpoint comparisons at the CF frequencies, E_state/E_direct was about0.46–0.50 for the fundamentals and0.35 forℓ=0,n=2: the unweighted direct bound was WORSE near the roots. The weighted bound reversed that comparison modestly:

| Mode | E_state/E_weighted at R20 | R30 | R40 |
|---|---:|---:|---:|
|0,0|1.732|1.735|1.768|
|2,0|1.754|1.832|1.851|
|0,2|1.053|1.084|1.097|

At the rounded off-root fundamental seeds, the unweighted direct bound was about2.0–2.2 times better than the terminal-state projection. For the more damped scalar overtone, the decaying component remains substantial at these radii, and that off-root advantage did not occur. These quantitative comparisons make no claim of dramatic or universal improvement.

## Numerical campaign and independence

Arithmetic40 decimal digits; Taylor order48; step fraction0.20;180 horizon terms; same rounded seeds and secant tolerance1e−29 as the earlier calculations. For each fixed R, N minimizes the uniform residual majorant over a stated finite scan on a frequency disk of radius0.01 centered at the rounded seed. The settings stay fixed during each root search. CF frequencies are used only for diagnostic evaluations and comparisons, not to seed roots. Radius is increased between completed independent searches until the conditional unweighted amplitude-tail target1e−30 is reached.

| Mode | R,N | Root seconds | Difference from finite-depth CF | Conditional unweighted tail bound |
|---|---:|---:|---:|---:|
|0,0|20,11|4.65|2.06e−11|6.13e−10|
|0,0|30,17|5.11|4.60e−16|1.36e−14|
|0,0|40,22|5.29|1.15e−20|3.45e−19|
|0,0|50,28|6.17|2.77e−25|8.29e−24|
|0,0|60,34|6.43|6.91e−30|2.04e−28|
|0,0|80,46|7.64|9.21e−32|1.30e−37|
|2,0|20,37|6.22|3.25e−36|7.15e−35|
|0,2|20,47|6.14|1.45e−21|5.14e−20|
|0,2|30,70|8.44|1.51e−28|9.06e−32|

The best weighted bounds at the selected final settings were3.42e−38,2.03e−35, and2.95e−32. In low-radius scalar fundamental tests, the weighted bound was close to the observed discrepancy: R20 gives approximately1.62e−10 versus observed |A_N|1.49e−10; R40 gives8.99e−20 versus8.29e−20. These are empirical checks, not proofs of numerical transport accuracy.

The final ℓ=0 comparison differences are limited largely by the finite-depth CF references: their double-depth frequency-correction estimates are9.21e−32 and1.51e−28. Thus agreement at that level does not measure the forward solver's exact error. Stored decimal frequencies are display precision; error-comparison fields were evaluated at the working precision before serialization.

A55-digit repeat and a Taylor56/step0.18 check were performed at the selected roots. The finer transport changed the approximate coefficient by about2.58e−36,6.16e−37, and3.37e−40, respectively. These checks are empirical. They do not turn the conditional terminal-tail formulas into a full certified root enclosure.

## Same-endpoint transport control and honest speed assessment

Forward and backward finite-P_N coefficients are mathematically identical at the SAME R,N by Wronskian conservation. The implementation comparison confirms that differences arise from numerical transport/conditioning, not a different exact truncated spectral function.

| Mode, same R,N | Forward root seconds | Backward root seconds |
|---|---:|---:|
|0,0 at80,46|7.64|7.08|
|2,0 at20,37|6.22|5.81|
|0,2 at30,70|8.44|8.46|

No intrinsic speed advantage from transport direction alone is established. Relative to the earlier conservative endpoint choices, smaller radii can reduce work, and the forward formulation provides h_R at no extra solve, making coefficient-specific a posteriori tail assessment convenient. It may also improve numerical conditioning in individual cases. Those are narrower, defensible advantages. The weighted coefficient bound produces a modest quantitative improvement at identical endpoints; neither forward matching itself nor the unweighted direct bound justifies a broad novelty claim.

Reproduce:
```
python3 analysis/qnm_forward_coefficient.py
python3 analysis/qnm_forward_coefficient_comparisons.py
```
