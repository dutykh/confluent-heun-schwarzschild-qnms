# Three-point Cauchy–Rouché root isolation

This note gives the analytic proof implemented in `cert/isolate.py`, the independent uniform bound used for the Schwarzschild incoming amplitude, and an independent review of the interval evaluator. The certificate establishes a local zero count. It does not establish a global ordering of quasinormal overtones.

## 1. Exact rational root test

Let F be holomorphic on a neighborhood of the closed disk |z−c|≤R, and assume |F|≤M on its outer circle. Let 0<h<R be real. Write F(c+w)=Σ a_n w^n. Cauchy's coefficient inequality gives |a_n|≤M/R^n. Consequently

D_h = [F(c+h)−F(c−h)]/(2h),

|D_h−a_1| ≤ Σ_{j≥1} M h^(2j)/R^(2j+1)
             = M h²/[R(R²−h²)].

Suppose three rigorous complex rectangles enclose F(c), F(c+h), F(c−h). Interval subtraction and division give a rectangle enclosing D_h. Choose its rational midpoint d, and let ε_s bound its circumradius. Let U_0 bound the modulus of the center-value rectangle. Define

δ = ε_s + M h²/[R(R²−h²)].

For 0<r<R and |w|=r,

|F(c+w)−d w| ≤ U_0 + δ r + M r²/[R(R−r)].

Thus the strict inequality

U_0 + δ r + M r²/[R(R−r)] < |d| r

proves by Rouché that F has exactly one zero, counted with multiplicity, inside |z−c|<r. The comparison function d(z−c) has one zero and d≠0 follows from the strict inequality. The enclosed zero is therefore simple. The open circle contains no zero of F on its boundary.

`make_certificate` performs this test with `fractions.Fraction`. It replaces |d| by a rational lower bound and all other square roots by rational upper bounds. These bounds come from integer square root and are explicitly checked by squaring. Decimal hints in JSON are for display only. The center and all radii/sample offsets are exact rationals. The three input rectangles must refer to precisely these three points; `check.py` verifies that requirement from the saved exact frequency metadata.

No numerical derivative is silently assumed exact. The finite difference error is included explicitly. No small-circle sampling or argument tracking is used.

## 2. The analytic Schwarzschild amplitude

Use dimensionless Ω=2Mω and x=z+Log(z−1), with z>1 initially and continuation into the upper half of the z plane. For the massless scalar equation, L=ℓ(ℓ+1),

f=(z−1)/z,  v=L/z²+1/z³,  U=fv,

u_xx+(Ω²−U)u=0.

The unit ingoing horizon solution is normalized by u_H e^(iΩx)→1 at the horizon. Set h=e^(iΩx)u_H and j=h_x. Then

h_xx−2iΩh_x−U h=0,

d/dz (h,j)^T = [[0,z/(z−1)],[v,2iΩ z/(z−1)]](h,j)^T.

The incoming coefficient is A_in=lim_{s→∞}h(a+is), where a>1 and Re Ω>0. This is the continuation of the scattering normalization u_H=A_out u_out+A_in u_in with unit asymptotic phases e^(±iΩx). Along a vertical ray the outgoing contribution to h decays, while the incoming contribution tends to its coefficient. The normalization contains no extra constant phase in this definition of x.

For Ω=p−iq, p>0 and q≥0,

d/ds Im[Ωx(a+is)]
 = [p s²−q s+p a(a−1)]/[(a−1)²+s²].

If q²<4p²a(a−1), this is strictly positive on the whole ray. The forward free kernel

K(s,t) = [exp(2iΩ(x(a+is)−x(a+it)))−1]/(2iΩ)

therefore obeys |K|≤1/|Ω| for s≥t≥0. The forward Volterra equation gives

sup_{s≥0}|h(a+is)| ≤ (|h(a)|+k|h_x(a)|) exp(k V_a),

where k≥1/|Ω| and

V_a = ∫_0^∞ |v(a+is)|ds ≤ πL/(2a)+1/a².

The same bound applies to |A_in|. The limit exists by dominated convergence in the Volterra representation; the exponentially decaying homogeneous derivative and the integrable potential give the required limits. On a compact frequency domain with the strict phase inequality, all estimates are uniform, so the Volterra construction and its limit are holomorphic in Ω. Horizon recurrence denominators are n(n−2iΩ); their zeros Ω=−in/2 lie outside Re Ω>0. The closed outer disks used here lie strictly in Re Ω>0 and admit open neighborhoods with the same exclusions and strict damping. This supplies the holomorphy required in Section 1.

These statements identify local zeros of the continued incoming amplitude with the stated scalar QNM boundary conditions. The conventional labels ℓ,n attached to candidate disks are benchmark labels; the local certificate does not prove the global spectral ordering implied by n.

## 3. Rational uniform outer bound

`schwarzschild_outer_bound` is independent of the numerical transport and the infinity asymptotic polynomial. Let the outer frequency disk have center p_0−iq_0 and radius R_Ω. Its enclosing rectangle has p≥p_min=p_0−R_Ω and q_min=q_0−R_Ω≤q≤q_max=q_0+R_Ω. The routine checks p_min>0, q_min≥0, and

q_max²<4p_min²a(a−1).

Use Ω_min≤|center|−R_Ω and Ω_max≥|center|+R_Ω from rational square-root bounds. Set k≥1/Ω_min and Λ≥2Ω_max.

### Uniform horizon coefficients

At z_0=3/2 put w=1−z and λ=2iΩ. Write

h=z exp(λ(z−1)) H(w),   H(w)=Σ h_n w^n,   h_0=1,

n(n−λ)h_n = [n(n−1)−λ²−L]h_{n−1}−λ(n−1−λ)h_{n−2},

with h_{−1}=0. For integer n, set d_n=distance(n,[2q_min,2q_max]). Then

|n−λ|≥sqrt(d_n²+4p_min²).

The code recursively majorizes |h_n| by nonnegative rational b_n, replacing the denominator square root by a rational lower bound and using numerator majorants

[n(n−1)+Λ²+L]b_{n−1}+Λ(n−1+Λ)b_{n−2}.

For J=N+1>Λ and all n≥J, the recurrence coefficients obey

A=[1+(Λ²+L)/J²]/(1−Λ/J),
B=[(Λ/J)(1+Λ/J)]/(1−Λ/J).

Taking ρ≥A+B is sufficient: A≥1 and (A+B)²−A(A+B)−B=B(A+B−1)≥0. With

C≥max(b_N/ρ^N,b_{N−1}/ρ^(N−1)),

induction gives |h_n|≤Cρ^n for all n≥N−1. For t=ρ/2<1, the omitted value and z-derivative tails at |w|=1/2 are bounded respectively by

C t^(N+1)/(1−t),
Cρ t^N[(N+1)−Nt]/(1−t)².

Finite coefficient sums plus these tails give H_* and H'_*. The gauge and its derivative yield

G_*=(3/2)exp(2q_max(1/2)),
|h(z_0)|≤G_*H_*,
|h_x(z_0)|≤(G_*/3)[H'_*+(Λ+2/3)H_*].

### Real segment to the anchor

If a>z_0, the x-system has matrix [[0,1],[U,2iΩ]]. Since U≥0 on this real segment, its induced infinity norm is at most max(1,2Ω_max)+U. Hence both components at a are bounded by

max(|h(z_0)|,|h_x(z_0)|)
× exp{max(1,2Ω_max)[a−z_0+log((a−1)/(z_0−1))]
       +L(1/z_0−1/a)+(1/z_0²−1/a²)/2}.

No propagation is used when a=z_0, preserving separate component bounds. Combining these bounds with Section 2 gives a rational upper bound M, rounded up to an integer.

All exponential and logarithmic upper bounds in this routine are evaluated with positive rational series and explicit geometric tails. Exponentials are range reduced to [0,1/2] and squared back with upward dyadic rounding. Logarithms use log x=2 atanh[(x−1)/(x+1)] after reduction to [1,2]. The default bound π<22/7 follows, for example, from the strictly positive integral ∫_0^1 x⁴(1−x)⁴/(1+x²)dx=22/7−π. No experimental transcendental interval routine is needed in this uniform outer bound.

## 4. Independent review of the point evaluator

### Horizon and finite transport

The implemented recurrence in `transport.py` is equivalent to the recurrence above. Its separate tail scheme uses |h_n|≤Cρ^(−n), 1/2<ρ<1, with Abar·ρ+Bbar·ρ²<1. It is a different but equally valid parametrization of a geometric majorant. Its value and derivative tail formulas are correct.

For the state (h,h_x), the fundamental matrix Φ of a finite z step satisfies ||Φ||_∞≤exp(Mr) on an analytic disk of radius r, where

M=max(b_max,V_max+2|Ω|b_max),
b_max=1+1/(d_1−r),
V_max=L/(d_0−r)²+1/(d_0−r)³.

The scalar polynomial ODE used to calculate its columns is exactly

z²(z−1)h''+(z−2iΩz³)h'−(Lz+1)h=0.

If q=|Δz|/r<1 and the Taylor polynomial has degree n, the first-row remainder in its induced row norm is at most S q^(n+1)/(1−q), S=exp(Mr). The derivative-row remainder is at most

S/r · q^n[(n+1)/(1−q)+q/(1−q)²],

multiplied by |f(z+Δz)| to recover h_x. There is no missing factor two: the Cauchy estimate applies to the induced matrix norm, hence directly to each row's operator norm. It must not instead be justified by summing two independent scalar entry bounds.

Point arithmetic in that module chooses exactly representable centers and step sizes only. Each center-rounding displacement is included when a rectangle is converted to a complex disk. All majorants and propagated disk radii use outward interval arithmetic. Taylor centers and steps may be interval enclosures of exact geometry; dependency overestimation is safe.

### Infinity coefficient and residual correction

For P_N=Σ_{n=0}^N a_n z^(−n), a_0=1,a_{−1}=0,

a_{k+1}=([k(k+1)−L]a_k−k²a_{k−1})/[2iΩ(k+1)].

Its exact residual is R_N=c_N z^(−N−2)+d_N z^(−N−3), with the two coefficients implemented in `amplitude.py`. The exact identity

d/dz{(2iΩh−h_x)P_N+h fP'_N}=h R_N

implies

A_in−A_N=(2iΩ)^(−1)∫_{z_R}^{∞}h R_N dz.

Let b=h_x(z_R)/(2iΩ), a=h(z_R)−b. Subtracting the free forward solution a+b exp(2iΩΔx) gives the corrected approximation

A_corr=A_N+[a I_0+b I_e]/(2iΩ),

with |A_in−A_corr|≤(|a|+|b|)F_R[exp(kV_R)−1]/(2|Ω|).

The implemented collection of h and h_x coefficients is algebraically exact. On the tail, pR>q is a sufficient phase condition, and k=1/max(p_min,q_min) conservatively bounds 1/|Ω|. The replacement exp(x)≤1/(1−x), 0≤x<1, and exp(x)−1≤x/(1−x) is valid. Cartesian one-norms safely bound all complex moduli. Endpoint-state enclosures and integral enclosures are both included before adding the feedback-error rectangle.

### Quadrature

The residual integral without an exponential has an exact rational antiderivative. The exponential integral is evaluated by analytic Taylor quadrature. Its logarithmic phase is on the continuation consistent with the vertical ray; the ratio (z−1)/(z_R−1) has positive real part there. Local analytic disks avoid z=0 and z=1. The implemented coefficient recurrence, symmetric-panel cancellation of odd terms, first omitted even degree, and integrated Cauchy remainder are correct. The uniform exponential damping μ=2p_min−2q_max/R>0 and the infinite-tail bound are conservative. Real transcendental enclosures and principal logarithms are delegated to the MPFI bridge.

### Serialization and precision

`campaign.py` sets the interval precision before constructing the exact rational sample frequencies. Its transport precision context restores that same precision. A point frequency is enclosed, rather than rounded and treated as exact. The output serializes each finite binary endpoint tuple to an exact integer-over-power-of-two rational string; the checker parses those strings as exact Fractions. It never reads displayed decimal approximations as certificate data. Using different valid residual formulas at the center and the two shifted points causes no problem, since every returned rectangle encloses the same analytic incoming amplitude at its stated point.

`cert/elementary.py` imports and exports finite binary endpoints exactly and uses MPFI/MPFR for elementary functions. Certification remains conditional on the standard correctness of Python integer/Fraction operations, mpmath's directed rational interval operations, and the installed MPFI/MPFR/GMP implementation. This is a computer-assisted proof with explicit arithmetic enclosures, not a formal proof verified by a proof assistant.

## 5. Reproducibility and scope

The exact-rational checker tests in `cert/test_isolate.py` cover successful one-zero examples, a two-zero example that must fail, widened sample rectangles that must fail, exact endpoint transfer, and the outer-bound phase rejection. Eight tests passed during this independent audit. The saved mode certificates are generated from `python3 -m cert.check MODE --radius 1e-28`; point enclosures are generated separately by `python3 -m cert.campaign MODE POSITION`.

A successful mode certificate proves exactly one simple zero in its stated dimensionless Ω disk. It does not alone prove absence of other zeros outside that disk, completeness of a spectrum, global overtone numbering, or a corresponding statement for a different spin or background geometry.
