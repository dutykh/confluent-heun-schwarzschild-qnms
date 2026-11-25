# Independent multiprecision contour diagnostic

Files:
- `qnm_contour_mp.py`: mpmath implementation using local Taylor propagation; no SciPy import and no import of the double-precision contour solver.
- `qnm_contour_mp_overtone_run.py`: additional scaled-horizon ℓ=0,n=2 campaign.
- `qnm_contour_mp_results.json`: full settings, states, amplitudes, analytic tail estimates, roots, and precision checks.

Reproduce in order:
```
python3 analysis/qnm_contour_mp.py --campaign
python3 analysis/qnm_contour_mp_overtone_run.py
```
The first command recreates the two fundamental-mode records; the second appends/replaces the overtone record.

## Independent propagation

For m=exp(−iΩx)u_out and v=m_x, the z-system is

 m_z=bv,  v_z=V m−2iΩ b v,
 b=z/(z−1),  V=ℓ(ℓ+1)/z²+1/z³.

The exact Taylor coefficients about z0 are b0=z0/(z0−1), b_k=(−1)^k/(z0−1)^(k+1) for k≥1, and

 V_k=(−1)^k[(k+1)ℓ(ℓ+1)/z0^(k+2)+(k+1)(k+2)/(2z0^(k+3))].

If A(z)=ΣA_k(z−z0)^k, then y_(n+1)=Σ_(k=0)^n A_k y_(n−k)/(n+1). Horner evaluation advances the solution along the vertical ray. Each step has |h|≤step_fraction min(|z0|,|z0−1|) and |h|≤1/|Ω|. This is distinct from the double-precision DOP853 propagation.

The horizon is initialized from its convergent Heun recurrence at z0=1.5 with 180 terms. For the higher overtone, use h=exp(iΩx)u_H=z exp(2iΩ(z−1))H_H. Then (h,h_x) obey the same rational system with Ω replaced by −Ω, permitting real-segment Taylor propagation to c=5. The determinant is D=2iΩ h m+h m_x−h_x m. This scaling avoids carrying unnecessary large endpoint exponentials separately, but does not eliminate all cancellation.

## Calculations

Frequency convention Ω=2Mω; time dependence exp(−iωt). Baseline arithmetic is 40 decimal digits, with a 55-digit repeat at the final root. The final root solves the order48, step_fraction0.20 propagation with outgoing inverse-power polynomial order65 at R=140. Analytic bounds refer to exact outgoing terminal data and exact propagation; Taylor truncation and floating-point errors are NOT interval-certified.

At the independently computed continued-fraction frequencies, the absolute incoming amplitudes were:

| ℓ,n | c | Taylor order32, fraction0.25, R100,N50 | Order40, fraction0.25, R140,N65 | Order48, fraction0.20, R140,N65 |
|---|---:|---:|---:|---:|
| 0,0 |1.5|5.07305e−21|5.07573e−26|6.66066e−31|
| 2,0 |1.5|5.28775e−20|4.87206e−25|1.40401e−34|
| 0,2 |5.0|9.61856e−18|1.22712e−23|1.33151e−27|

The final matching-amplitude tail-only bounds are respectively 5.57254e−37, 5.57386e−67, and 1.11926e−73. These very small endpoint bounds do not bound the full implemented numerical error; propagation and arithmetic must be assessed separately.

The secant roots were obtained from explicitly rounded seeds (.22−.21i), (.97−.19i), and (.15−1.20i), with a second seed displaced by .001+.001i. The continued-fraction values were used only for diagnostic comparisons, never as root seeds.

| ℓ,n | Computed contour Ω (digits shown are output, not a claim of certification) | Difference from finite-depth CF reference |
|---|---|---:|
|0,0|0.22090987816083937175092301233587754 −0.20979143417376191756347813371075195 i|9.21472e−32|
|2,0|0.96728774442142597345555887030351719 −0.19351755195657572531730746670174373 i|3.02591e−35|
|0,2|0.15148387107035152594044328823929967 −1.2021571800716070281873120374184572 i|1.47908e−28|

The ℓ=0 fundamental CF reference itself has an estimated double-depth correction9.21456e−32; the overtone reference has correction1.51468e−28. Thus their comparison differences are largely limited by the finite-depth reference. These are local Newton diagnostics, not certified CF error bounds.

At the fixed final contour roots, 40-to55-digit amplitude changes were3.03e−41,3.60e−41, and8.19e−35. The overtone's greater arithmetic sensitivity is visible, but remains far below its double-precision discrepancy. As an additional overtone check, at the same root change c from5 to5.4, Taylor order48 to56, and step_fraction0.20 to0.18, using55 digits. The resulting |A_in| is3.57155e−29. This checks a different matching ray and propagation accuracy; it is not a rigorous frequency-error bound.

The evidence supports substantially more than20-digit agreement of these selected contour computations with an independent recurrence method. It does not prove a uniform performance guarantee, an interval-certified determinant error, completeness of the spectrum, or validity of every proposed contour deformation. The mathematical hypotheses of the analytic tail estimates remain separate requirements.
