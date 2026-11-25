# Explicit free-tail correction: sharper coefficients at added cost

Files: `qnm_free_correction.py` and `qnm_free_correction_results.json`. Earlier deliverables are unchanged. Run `python3 analysis/qnm_free_correction.py` to reproduce the nine fixed-endpoint comparisons and two independent corrected-root searches.

## Exact correction and remainder

For the forward horizon solution at z_R=c+iR, set

 a_H=h_R−hx_R/(2iΩ), b_H=hx_R/(2iΩ), Q=|a_H|+|b_H|.

Let the exact residual of P_N be

 r_N(z)=C_N z^(−N−2)+D_N z^(−N−3).

The free tail is h_free(z)=a_H+b_H exp(2iΩ[x(z)−x(z_R)]). Its explicit defect integral gives

 A_corr=A_N+[a_H I0+b_H Iexp]/(2iΩ),
 I0=C_N z_R^(−N−1)/(N+1)+D_N z_R^(−N−2)/(N+2),
 Iexp=∫_(z_R)^(∞ along ray) exp(2iΩ[x(z)−x(z_R)])r_N(z)dz.

The sign in I0 is positive: the primitive−z^(1−p)/(p−1) vanishes at the upper endpoint. The exact forward Volterra series gives |h−h_free|≤Q expm1(kV_R), hence

 |A_in−A_corr|≤Q expm1(kV_R)F_R/(2|Ω|).

This proof and the orientation/normalization were independently checked. The estimate assumes exact transported h_R,hx_R and exact I0,Iexp. State, Taylor, arithmetic, and quadrature errors must be added separately. No full numerical certificate is asserted.

## Quadrature and independent reference

Iexp is integrated upward with dz=i du. The exponential uses the stable phase

 Δx=iu+log1p(iu/(z_R−1)).

The integrand is normalized by the magnitudes of its residual terms before quadrature, preventing tiny residuals from being mistaken for numerical zero. A phase/power decay scale rescales the integration variable. The baseline uses40 digits and breakpoints0,1,4,12,∞; the repeat uses55 digits and0,0.5,2,8,24,∞. A separate no-exponential quadrature reproduced the elementary I0 with relative discrepancies below5e−57 in the tested R20 cases. The largest estimated coefficient change induced by the Iexp precision/breakpoint repeat was below1e−50. These are convergence checks, not certified quadrature-error bounds.

The comparison does NOT treat A_in at a finite-depth CF frequency as exactly zero. At that same frequency, a separate55-digit backward contour evaluation with Taylor56, step0.18, R140,N65 supplies a higher-resolution numerical reference. This reference is itself empirical. All corrected and uncorrected fixed-endpoint evaluations otherwise use the same40-digit transport, Taylor48, step0.20,180 horizon terms, and identical R,N.

## Fixed-endpoint improvement

For the scalar fundamental:

| R,N | Uncorrected coefficient discrepancy | Corrected discrepancy | Conditional corrected remainder | Bound improvement over BEST uncorrected weighted bound |
|---|---:|---:|---:|---:|
|20,11|1.48942e−10|2.69691e−14|6.78401e−13|239×|
|30,17|3.32356e−15|1.79778e−19|6.67041e−18|539×|
|40,22|8.28975e−20|2.02167e−24|9.53253e−23|943×|

The observed discrepancy reductions were about5523×,18487×, and41004×. These are meaningful improvements at identical endpoints. The comparison uses the best earlier weighted bound, not the weaker unweighted forward bound.

Other cases:

| Mode,R,N | Uncorrected discrepancy | Corrected discrepancy | Conditional corrected remainder |
|---|---:|---:|---:|
|ℓ2,n0;20,37|1.50870e−35|6.24667e−37|5.11596e−36|
|ℓ2,n0;30,56|6.15638e−37|6.15638e−37|2.83901e−53|
|ℓ2,n0;40,76|6.15638e−37|6.15638e−37|2.23900e−70|
|ℓ0,n2;20,47|1.30167e−20|2.03996e−25|1.64571e−23|
|ℓ0,n2;30,70|2.21943e−32|4.38939e−34|1.31500e−35|
|ℓ0,n2;40,94|4.38169e−34|4.38169e−34|1.09303e−47|

The floors are essential to interpretation. The ℓ2 R30/R40 and overtone R40 errors do not follow their much smaller tail bounds: improving the analytic tail cannot remove finite-transport/arithmetic/reference error. The overtone R30 correction reaches a similar transport floor. Quadrature refinements are far smaller, supporting this diagnosis, but a fully certified separation of errors has not been performed.

## Roots and cost: no speed-superiority claim

Both corrected-root searches start from the same rounded seeds as before and use the same secant stopping tolerance. They do not use the reference frequency as the root seed.

| Mode and formula | R,N | Difference from CF reference | Root wall seconds |
|---|---:|---:|---:|
|ℓ0,n0 corrected|40,22|2.79689e−25|9.33|
|ℓ0,n0 plain, same endpoint|40,22|1.14685e−20|5.29|
|ℓ0,n0 plain, larger endpoint|50,28|2.77382e−25|6.17|
|ℓ0,n2 corrected|20,47|2.28116e−26|18.32|
|ℓ0,n2 plain, same endpoint|20,47|1.44592e−21|6.14|
|ℓ0,n2 plain, larger endpoint|30,70|1.51468e−28|8.44|

Thus the correction sharply improves fixed-endpoint coefficients and their conditional remainder bounds, but the present quadrature implementation is NOT cheaper than simply extending the plain integration in these examples. The larger-endpoint plain scalar-fundamental calculation reaches comparable accuracy faster; the plain overtone R30 calculation is both faster and more accurate. No statement of equal achieved accuracy is inferred merely from a common arithmetic setting or root stopping tolerance.

The strongest demonstrated result is consequently a sharper explicit corrected coefficient with a controlled exact-arithmetic residual, obtained at additional numerical cost. This is a standard Volterra-type defect correction, not a claimed priority result. Its relevance to a new research contribution depends on broader novelty analysis and on a future complete error budget or a genuinely advantageous evaluation strategy.
