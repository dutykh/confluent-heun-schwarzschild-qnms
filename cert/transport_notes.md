# Validated horizon initialization and finite transport

`transport.py` implements finite transport of the scalar Schwarzschild state
(h,v)=(exp(iΩx)u_H,h_x), with x=z+log(z−1). It uses directed `mpmath.iv`
arithmetic for algebra and the MPFI/MPFR elementary-function bridge in
`cert.elementary`. The resulting enclosures include Frobenius truncation,
local Taylor truncation, arithmetic rounding, input-rectangle uncertainty,
and enclosure propagation. They do not by themselves certify a QNM: an
outgoing-tail estimate and root-isolation argument are separate layers.

## API

```
transport_horizon(omega, ell, anchor='1.5', R=80,
                  dps=100, order=96, step_fraction='0.18',
                  horizon_terms=220, return_steps=False)
```

The return value is a dictionary containing `h`, `hx` as rectangular interval
complex numbers; `centers` and Euclidean complex-disk `radii`; the target,
precision/settings, step count, and timing. The method accepts a complex
frequency rectangle, not merely its midpoint. The path is 1.5→anchor→anchor+iR.
`horizon_enclosure` returns the initial enclosure at1.5 and its explicit
Frobenius-tail metadata. `exact_box` serializes interval endpoints losslessly
as dyadic rational strings.

## Horizon tail

Write H_H(z)=Σc_n(1−z)^n. The exact recurrence is

 c_n={[n(n−1)+4Ω²−L]c_(n−1)−[2iΩ(n−1)+4Ω²]c_(n−2)}/[n(n−2iΩ)].

The interval recurrence rejects a denominator whose enclosure can contain
zero. Let O≥|Ω| on the input rectangle, λ=2O, C=4O²+L, n0=N+1, and

 Abar=(1+C/n0²)/(1−λ/n0),
 Bbar=(λ/n0+C/n0²)/(1−λ/n0).

These bound the absolute recurrence coefficients for every n≥n0. Choose
ρ∈(1/2,1) and verify Abarρ+Bbarρ²<1 using outward interval arithmetic.
Then, with

 CH=max(|c_(N−1)|ρ^(N−1),|c_N|ρ^N),

induction gives |c_n|≤CHρ^(−n) for all subsequent n. At |w|=1/2,
q=|w|/ρ,

 |tail H|≤CH q^(N+1)/(1−q),
 |tail H_z|≤(CH/|w|)q^(N+1)[(N+1)/(1−q)+q/(1−q)²].

The exact horizon gauge is
h=z exp(2iΩ(z−1))H_H and
v=f z exp(2iΩ(z−1))[H_H'+(2iΩ+1/z)H_H].
The corresponding tail errors and interval arithmetic widths are enclosed
before transport begins.

## Fast fundamental-matrix Taylor recurrence

The first-order system is

 h_z=bv, v_z=Vh+2iΩbv,
 b=z/(z−1), V=L/z²+1/z³.

Elimination gives the polynomial equation

 z²(z−1)h''+(z−2iΩz³)h'−(Lz+1)h=0.

Its Taylor coefficient recurrence has only five adjacent terms. Both
fundamental columns are therefore computed in O(N) operations per step,
instead of the O(N²) rational-coefficient convolution. Initial columns are
(h,v)=(1,0) and(0,1); h_z=v/f. After evaluating the polynomial and derivative,
v=f h_z is evaluated with interval arithmetic.

## Cauchy remainder in the correct operator norm

At a center z0, take a certified disk radius r less than the lower bounds
on |z0| and |z0−1|. On that disk,

 |b|≤1+1/(d1−r),
 |V|≤L/(d0−r)²+1/(d0−r)³.

An upper bound M for the induced infinity norm of the full system matrix is

 M=max(bmax,Vmax+2|Ω|bmax).

The fundamental matrix satisfies ||Φ||∞≤S=exp(Mr) on the circle. Cauchy's
formula in this matrix norm bounds the first row in its induced l1 norm;
there is no extra factor two from treating its columns separately. For
q=|step|/r<1 and a degree-N polynomial,

 E_h,row≤S q^(N+1)/(1−q),
 E_hz,row≤(S/r)q^N[(N+1)/(1−q)+q/(1−q)²],
 E_v,row≤|f(z1)|E_hz,row.

These row-operator errors multiply the infinity norm of the incoming state
once. All constants, norms, exponentials, and final radii are rounded
outward. Ordinary point arithmetic merely chooses exactly representable
centers and steps; interval evaluation then encloses those choices.

## Wrapping control

Each state component is represented as a complex center plus a Euclidean
disk radius. After a local interval polynomial matrix T,

 new_radius_i ≤ radius(T_i center)
              +Σ_j sup|T_ij| radius_j
              +E_i,row max_j(|center_j|+radius_j).

Thus harmless complex phase rotations do not repeatedly inflate Cartesian
rectangles. The final disk is enclosed by a rectangular interval for the
remaining certificate layers.

## Verification

`python3 -m cert.test_transport` passes:

- Exact static polynomial solutions h=zP_l(2z−1), l=0,1,2, along a complex path.
- Three physical frequency points, with a higher-order/higher-precision
  enclosure contained inside the default enclosure.
- A genuinely nonpoint frequency rectangle in the horizon initialization.

At default100 digits/order96/step fraction0.18, measured complex-disk radii are:

| Mode parameters | Endpoint | Radius h | Radius h_x | Approximate time |
|---|---|---:|---:|---:|
|Scalar l0 fundamental|1.5+80i|2.44e−56|1.14e−61|2.9s|
|Scalar l2 fundamental|1.5+20i|4.73e−49|3.85e−51|2.6s|
|Scalar l0 second overtone|5+30i|7.15e−45|1.45e−47|2.5s|

The exact figures/settings are in `transport_test_results.json`. Refining to
120 digits/order112/fraction0.16 reduced the radii by roughly fifteen orders
of magnitude, with the refined boxes contained in the default boxes. These
tests support implementation correctness; the enclosure justification is
the analytic bounds and directed arithmetic above. The recurrence,
Cauchy-row norm, horizon majorant, and rounding pathway were independently
reviewed by another agent.
