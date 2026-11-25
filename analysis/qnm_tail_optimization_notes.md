# Sharper outgoing-tail bounds and controlled cost comparisons

This extension changes none of the earlier contour or Leaver deliverables. The new files are `qnm_tail_optimization.py` and `qnm_tail_optimization_results.json`; the follow-up scripts reproduce the additional legacy-adaptive timing and global-residual diagnostics. Calculations below concern exact-arithmetic tail errors. Taylor truncation and floating-point transport remain empirical, and no frequency-error certificate is claimed.

## Exact geometry and residual cancellation

Write z=c+is, c>1, s≥R>0. For p>1,

 J_p(c,R)=∫_R^∞(c²+s²)^(−p/2)ds
 =c^(1−p) B_{c²/(c²+R²)}((p−1)/2,1/2)/2.

Thus the bound based on |z|≥R can be replaced by exact positive power-tail integrals. The finite inverse-power polynomial P_N has exact residual

 r_N(z)=C_N z^(−N−2)+D_N z^(−N−3),
 C_N=[N(N+1)−L]a_N−N²a_(N−1), D_N=−(N+1)²a_N.

The triangle bound is F_R≤|C_N|J_(N+2)+|D_N|J_(N+3). An additional cancellation-aware majorant follows without assuming cancellation persists under integration. Write C_N z+D_N=iC_N s+K_N, K_N=cC_N+D_N. Convexity of the norm on1/s∈[0,1/R] gives

 |iC_N+K_N/s|≤C*=max(|C_N|,|iC_N+K_N/R|),
 F_R≤C*(c²+R²)^(−(N+1)/2)/(N+1).

Take the smaller of these two independently valid bounds. This cancellation refinement uses the actual complex coefficients at a specified frequency. The frequency-disk selector instead uses positive coefficient majorants and makes no unsupported phase-cancellation assumption.

## A smaller uniform kernel constant

Let Ω=p−iq, p>0,q≥0. Along an upward vertical ray, Δx=U+iV has0≤U/V≤ρ, where ρ=1/[2√(c(c−1))]. Under qρ<p, set

 κ=(q/p+ρ)/(1−(q/p)ρ), S=√(1+κ²), k=S/[(1+S)|Ω|].

For w=2iΩΔx, Re w≤0 and |Im w|≤κ(−Re w). Writing a=−Re w,

 |1−e^w|≤min[S(1−e^(−a)),1+e^(−a)]≤2S/(1+S).

Therefore |K|≤k<1/|Ω|. On s≥R use the smaller slope bound
ρ_R=R/[R²+c(c−1)] when R≥√(c(c−1)), otherwise useρ. The code distinguishes k_tail from k_segment.

Crucially, differentiation does not acquire the smaller constant: |∂_xK|≤1. Thus the residual estimates are E0=k_tail F_R exp(k_tail V_R) and E1=F_R exp(k_tail V_R), NOT E1=|Ω|E0. With A=E0+k_segment E1, finite-segment propagation gives

 E_m≤A exp(k_segment V_segment),
 E_mx≤E1+(A/k_segment)[exp(k_segment V_segment)−1].

The scalar phase-cone derivation was independently checked by another agent.

## Decaying exponential in derivative estimates

The phase derivative obeys

 d Im[Ωx(c+is)]/ds=p+[p(c−1)−qs]/[(c−1)²+s²].

Its positive minimum gives |exp(2iΩΔx)|≤exp(−μ(t−s)). The minimum over a tail is evaluated at the larger ofR and (c−1)(p+√(p²+q²))/q, with q=0 handled by the limit. Put r=√(c²+R²), λ=R/r². Convexity gives |c+i(R+u)|≥r(1+λu), hence

 ∫_R^∞e^(−μ(s−R))|c+is|^(−p)ds
 ≤r^(−p)e^(μ/λ)E_p(μ/λ)/λ,

where E_p(a)=∫_1^∞e^(−at)t^(−p)dt. This is also bounded byJ_p. Let Fw and Vw denote the resulting weighted residual and potential bounds. A uniform tail value bound is

 E0≤min[k_tail F_R,(F_R+Fw)/(2|Ω|)] exp(k_tail V_R),
 E1≤min[F_R exp(k_tail V_R),Fw+E0 Vw].

The E0 used here bounds the entire tail, using decreasing power majorants, not just its endpoint. The segment derivative estimate is also compared with its exponentially weighted counterpart.

Numerical caution: at40 digits, mpmath's high-order expint evaluation lost positivity in one large-argument case. The implementation now evaluates the positive scaled expression with guarded extra precision, checks positivity, and was independently compared with positive normalized quadrature. Three upper-bound/quadrature ratios were1.0000000672,1.0000067061, and1.0004938343. This repairs numerical evaluation of the formula; it is not outward-rounded interval certification.

## Frequency domains and order selection

Each search uses the disk of radius0.01 around the rounded seed, wholly in ReΩ>0,ImΩ<0. Set p_min=ReΩ_c−r, q_max=−ImΩ_c+r, Ω_min=|Ω_c|−r. The phase bounds use these conservative rectangle extrema. Uniform coefficient majorants satisfy

 A0=1,
 A_(k+1)=[|k(k+1)−L|A_k+k²A_(k−1)]/[2Ω_min(k+1)].

The common requested bound is on |δm|+|δm_x|/Ω_min at the anchor:10^(−25). It is NOT a common certified amplitude or frequency error. The selector scans radii20,25,… and asymptotic orders, takes the first feasible radius, and then the smallest feasible order. It also records the best order at that radius. This is a conservative grid selector, not a proof of globally minimal wall cost. Weighted refinements are verified after the conservative selection. At the preceding grid radius, additional refinements and the global-residual bound still failed the target in all three examples.

| Mode ℓ,n | Anchor | Legacy adaptive R,N | Improved adaptive R,N | Best searched N at improved R | Improved uniform state bound |
|---|---:|---:|---:|---:|---:|
|0,0|1.5|105,49|105,47|61|5.95e−26|
|2,0|1.5|40,47|35,62|66|6.54e−26|
|0,2|5.0|30,40|30,39|70|6.09e−26|

At exactly the improved selected settings, the pointwise state bounds improved by factors2.65,41.9, and2.58 relative to the old formulas. Exact geometry/combined residual alone was a modest part of the improvement; the phase-cone constant was particularly useful forℓ=2.

## Controlled numerical comparisons

All transport/root comparisons use the same40-digit arithmetic, Taylor order48, step fraction0.20,180 horizon terms, rounded seeds, and secant stopping tolerance. Endpoint settings stay fixed during each root search. The continued-fraction frequencies are used only for independent comparisons. Times are illustrative wall measurements from this runtime, not statistical performance guarantees.

| Mode | Outgoing steps old140/65 → legacy adaptive → improved adaptive | Root seconds old → legacy → improved |
|---|---|---|
|0,0|57 →47 →47|10.39 →8.63 →8.79|
|2,0|147 →49 →44|29.72 →9.87 →8.74|
|0,2|170 →37 →37|34.83 →8.10 →8.21|

Most of the speedup comes from removing the old fixed endpoint's overresolution. Beyond that ordinary adaptation, the sharper bound gave a modest additional root-time gain forℓ=2 (about13%). No extra measurable speed gain is established in the other two cases. Mean single-evaluation timing forℓ=2 was noisy enough to favor the legacy adaptive setting despite fewer improved steps; both raw repetitions are retained. Selector setup costs were0.15–0.66 seconds, separately recorded.

| Mode | CF frequency difference legacy adaptive | CF difference improved adaptive | Actual-root improved amplitude-tail bound |
|---|---:|---:|---:|
|0,0|9.21e−32|9.21e−32|4.34e−28|
|2,0|2.76e−35|2.87e−35|9.49e−29|
|0,2|2.04e−28|1.32e−28|6.78e−26|

These amplitude-tail bounds insert numerically computed horizon factors into exact analytic bounds and are therefore diagnostics, not certificates covering horizon/transport arithmetic. Forℓ=0,n=2 the state-to-amplitude norm amplification is about2.50e6, showing why a uniform state tolerance alone cannot establish equal amplitude accuracy. Frequencies agree empirically much better than10^(−20) in these selected examples, but the CF references and contour propagation also have their own truncation errors. The compared methods are not claimed to deliver identical achieved precision merely because they meet the same state-tail target.

## Additional single-global-residual and direct-coefficient checks

Join the exact finite-segment propagation toP_N on the tail; the joined approximation is C1 and its residual is supported only pastR. With one full-ray constantk,

 E_m≤kF_R exp(kV0), E_mx≤F_R exp(kV0).

For h=exp(iΩx)u_H, the identity2iΩK−exp(2iΩΔx)=−1 yields the direct determinant estimate

 |δD|≤(|h|+k|h_x|)F_R exp(kV0).

The code records this and the alternative obtained by applying the determinant's linear functional to the global component bounds. The better global forms mildly improve the fundamental examples; they can be weaker than the separate tail/segment treatment for the damped overtone because the tail admits a much better phase cone than the full ray. Keep the minimum of independently valid bounds, not a blanket claim that the global form always wins. These added diagnostics do not change the timed searches or their settings.
