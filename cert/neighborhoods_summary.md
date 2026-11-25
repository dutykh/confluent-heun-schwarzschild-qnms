# Uniform enclosures of the Schwarzschild incoming amplitude and its frequency derivative

For each of the nine fixed exact centers in the campaign, put r=10⁻²⁸, s=10⁻⁴ and D=10⁻², in Ω=2Mω. All nine outer disks of radius D lie in the holomorphy chart specified in the original scalar certificate. Their respective upper bounds M on |F|, where F=A_in, are obtained independently from the horizon recurrence, a real-segment Gronwall estimate and a forward Volterra estimate. The new certificate reuses these bounds, and the exact coefficient enclosures at the nine centers; it adds 72 outward-rounded point evaluations. Eight new samples per mode use the same certified finite transport and uncorrected residual functional as the earlier frequency-shifted evaluations. The corrections and truncation errors of these enclosures enter all subsequent coefficient intervals. All interpolation, error comparisons and decimal publication bounds are obtained from rational arithmetic.

## The interpolation construction

In t=(Ω−c)/s let G(t)=F(c+st). Its known analytic disk has radius T=D/s=100. Use nodes u_j=j/4 (j=−4,…,4), write Q(t)=∏(t−u_j)=t∏_{j=1}⁴(t²−(j/4)²), and let P(t)=Σ_{k=0}⁸b_k t^k interpolate the nine exact values of G. To enclose every b_k, form the nine Lagrange cardinal polynomials with rational coefficients and combine the sample rectangles with the corresponding signed rational weights. Denote the resulting power-coefficient rectangles by B_k. The values b_k belong to a common exact interpolant. Upper bounds using individual B_k are valid even though those boxes can be correlated.

The analytic interpolation formula, also valid at the nodes by continuation, is

    G(t)−P(t)=Q(t)H(t),
    H(t)=(1/(2πi))∮_{|ξ|=T}G(ξ)/[Q(ξ)(ξ−t)] dξ.

This follows either from residues or Cauchy's integral representation; the poles at the interpolation nodes recover their cardinal-polynomial terms. On the outer circle |Q(ξ)|≥m_Q=T∏_{j=1}⁴(T²−(j/4)²). On |t|≤1 put H₀=TM/[m_Q(T−1)] and H₁=H₀/(T−1). These bound |H| and |H′|. For Q=Σ q_k t^k set Q₀=Σ|q_k| and Q₁=Σ k|q_k|, so the interpolation errors for G and G′ are bounded by Q₀H₀ and Q₁H₀+Q₀H₁. The latter is divided by s for the physical frequency derivative F′.

Write β=mid(B₁)/s. Use the farthest-corner bound on the circumscribed radius rad(B₁) and the maximum modulus of each B_k. Then, uniformly throughout the closed disk |w|≤s,

    |F′(c+w)−β| ≤ e,
    e = [rad(B₁)+Σ_{k=2}⁸ k sup|B_k| + Q₁H₀+Q₀H₁]/s,
    |F(c+w)−F(c)−βw| ≤ e|w|.

The last inequality integrates the derivative along the segment from 0 to w. The same estimate holds for any two points of the disk, with |w₁−w₂| on the right. The available certified center amplitude satisfies |F(c)|≤U₀. Thus |F|≤U₀+(|β|+e)s on the entire fixed disk, |F′|≥|β|−e there, and at every point in the disk |F(c+w)|≥(|β|−e)|w|−U₀. Upper and lower bounds for |β| in this argument use directed integer square roots.

All nine strict inequalities (|β|−e)r>U₀ have passed. Consequently every closed annulus r≤|w|≤s is free of zeros and each inner boundary has a positive lower modulus. The absolute inner boundary bound in the table is the larger of that bound and the strictly positive original Rouché gap. Every circle of radius u in [r,s] satisfies |F(c+w)−βw|<|βw|, and Rouché gives exactly one zero counted with multiplicity inside. Simplicity also follows from |β|−e>0. The anchored Newton map z↦z−F(z)/β maps either closed disk (radius r or s) strictly into itself and has Lipschitz constant at most e/|β| on the closed s disk. In particular, this excludes every additional zero between the two radii. The argument concerns the nine individual disks and their intervening annuli. It makes no count in the regions between distinct mode neighborhoods.

## Published bounds

The columns below give the panel center β* (displayed rounded), e* covering F′ about β* after enclosing the center-rounding displacement, derivative modulus lower bound, coefficient modulus upper bound multiplied by 10⁴, inner boundary lower bound multiplied by 10²⁸, s-boundary lower bound multiplied by 10⁴, and contraction upper bound. Quantities following e* use the exact β archived in JSON.

| ℓn | Re β* | Im β* | e* | min modulus of F′ | 10⁴ max modulus of F | 10²⁸ inner lower | 10⁴ outer lower | contraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 00 | 5.7536 | 4.3754 | 0.00617 | 7.222 | 7.235 | 7.208 | 7.222 | 0.0009 |
| 01 | 8.0071 | 1.6586 | 0.00812 | 8.168 | 8.186 | 7.658 | 8.168 | 0.0010 |
| 02 | 8.9291 | 1.1457 | 0.01672 | 8.985 | 9.020 | 8.589 | 8.985 | 0.0019 |
| 10 | -2.6137 | 5.1346 | 0.00594 | 5.755 | 5.768 | 5.734 | 5.755 | 0.0011 |
| 11 | 2.0535 | 3.6338 | 0.00445 | 4.169 | 4.179 | 4.160 | 4.169 | 0.0011 |
| 12 | 4.2703 | 2.3931 | 0.00571 | 4.889 | 4.901 | 4.843 | 4.889 | 0.0012 |
| 20 | -4.0517 | -2.2612 | 0.00521 | 4.634 | 4.646 | 4.614 | 4.634 | 0.0012 |
| 21 | -1.7665 | 1.4168 | 0.00232 | 2.262 | 2.267 | 2.252 | 2.262 | 0.0011 |
| 22 | -0.0016 | 2.0510 | 0.00225 | 2.048 | 2.054 | 2.046 | 2.048 | 0.0011 |


## Verification

- Recompute all new samples and inequalities using `python -m cert.neighborhoods --workers 5`.
- Recheck all inequalities from archived exact point enclosures using `python -m cert.neighborhoods --check-saved`. The latter reconstructs the outer majorant and the previous Rouché test from its original samples.
- Run exact polynomial coefficient recovery and an analytic rational-function remainder test using `python -m cert.test_neighborhoods`.
- The files `cert/neighborhood_samples.json` and `cert/results_neighborhoods.json` store the 72 additional samples and all nine exact rational enclosures. Decimal publication bounds use outward integer division.
