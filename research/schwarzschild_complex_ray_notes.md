# Outgoing Schwarzschild scalar solution on a complex ray

This note derives the estimates directly; it makes no priority claim for complex contour continuation or Jost/Volterra methods.

Write
\[
f=(z-1)/z,\quad x=z+\log(z-1),\quad
v(z)=L/z^2+1/z^3,\quad U=fv,\quad L=\ell(\ell+1),\ \ell\in\mathbb Z_{\ge0}.
\]
Use the logarithm continued from real z>1 into the upper half-plane. Set m=e^{-i\Omega x}u_out. Then
\[
m_{xx}+2i\Omega m_x-U m=0,
\qquad f m_{zz}+(f'+2i\Omega)m_z-vm=0.
\]
Its outgoing integral equation is
\[
m(z)=1+\int_z^{\infty_{\rm ray}}K(z,t)v(t)m(t)\,dt,
\qquad K(z,t)=\frac{e^{2i\Omega(x(t)-x(z))}-1}{2i\Omega}.
\]
All derivatives called m_x below are with respect to x, not z or the ray parameter.

## 1. Admissible rays and exact monotonicity condition

Take z(s)=r+e^{i\theta}s, r>1, s≥0, and a=r−1. For \Omega=p−iq with p>0,q≥0 let A=Im(\Omega e^{i\theta})=p sin\theta−q cos\theta. Then
\[
\frac{d}{ds}\operatorname{Im}(\Omega x(z(s)))
=A+\frac{aA-qs}{s^2+2a\cos\theta\,s+a^2}.
\]
For A>0 this is nonnegative for every s≥0 precisely when
\[
q\le 2A\left(a\cos\theta+\sqrt{a(a+1)}\right).
\]
The algebra is just nonnegativity of
A s²+(2aA cos\theta−q)s+A a(a+1) on s≥0. Restricting 0<\theta≤π/2 keeps the contour in the upper-right half-plane and avoids both finite singularities and the logarithm cut.

The simplest choice is vertical, \theta=π/2:
\[
\frac{d}{ds}\operatorname{Im}(\Omega x(r+is))
=\frac{p s^2-q s+p a(a+1)}{a^2+s^2}.
\]
Thus monotonicity holds iff
\[
q/p\le2\sqrt{r(r-1)},\qquad
r\ge\frac{1+\sqrt{1+(q/p)^2}}2.
\]
A strict inequality gives an open frequency region with margin. Equality is still sufficient for the pointwise integral estimates. Monotonicity implies |exp(2i\OmegaΔx)|≤1 for every ordered pair on the ray, hence |K|≤1/|\Omega|.

For a fixed r the vertical construction with this simple bound is holomorphic at least on
\[
\mathcal W_r=\{\Re\Omega>0,\ \Im\Omega>-2\sqrt{r(r-1)}\Re\Omega\}.
\]
The upper-right quadrant is included. Different r can be used pointwise, but choose a common fixed r over a compact frequency region when proving bounds or holomorphy of a numerical approximant.

## 2. Explicit potential-integral majorants

On the vertical ray, |dt|=ds and |v(r+is)|≤L/(r²+s²)+(r²+s²)^{-3/2}. Set
\[
V(s)=\frac Lr\left(\frac\pi2-\arctan\frac sr\right)
+\frac1{r^2}\left(1-\frac{s}{\sqrt{r^2+s^2}}\right).
\]
Then \int_s^∞|v(z(t))|dt≤V(s), and V(0)=πL/(2r)+1/r².
For general 0<\theta≤π/2 the analogous exact majorant is
\[
V_\theta(s)=\frac L{r\sin\theta}\left[\frac\pi2-
\arctan\frac{s+r\cos\theta}{r\sin\theta}\right]
+\frac1{r^2\sin^2\theta}\left[1-
\frac{s+r\cos\theta}{\sqrt{s^2+2rs\cos\theta+r^2}}\right].
\]
At s=0 it becomes L\theta/(r sin\theta)+1/[r²(1+cos\theta)]. The θ→0 limit recovers L/(r+s)+1/[2(r+s)²].

## 3. Entire right-frequency-half-plane continuation at fixed r

Monotonicity is useful for efficient bounds but not necessary for existence. On the same vertical ray let p=Re\Omega>0, q=max(0,−Im\Omega), and d=t−s≥0. The explicit logarithm gives
\[
|e^{2i\Omega(x(z(t))-x(z(s)))}|
\le e^{-2pd}(1+d/a)^{2q}.
\]
Indeed, the argument of a+is increases, and the modulus ratio is at most 1+d/a. The supremum M of the right side is
\[
M=\begin{cases}
1,&q\le pa,\\
\exp(-2q+2pa)[q/(pa)]^{2q},&q>pa.
\end{cases}
\]
Consequently |K|≤\kappa=(1+M)/(2|\Omega|). On every compact subset of Re\Omega>0 this bound is uniform. Neumann iterates obey (\kappa V(s))^n/n!, so the integral series converges locally uniformly in \Omega and is holomorphic on the entire right half-plane. No small-frequency or small-potential expansion is involved.

For \Omega in the first quadrant, deform the real outgoing contour through upper-right z rays. At every rotation angle the exponential term decays, while v=O(z^{-2}); the nonexponential arc contribution tends to zero. This can be applied inductively to the integral iterates: the previous iterate is holomorphic and has a uniform algebraic tail, so both terms of the split kernel have vanishing arc integrals. Uniform factorial bounds then pass equality to the sum. Thus the vertical construction agrees with the real-ray outgoing Jost solution in their first-quadrant overlap. The identity theorem identifies its lower-right values as the unique analytic continuation from that overlap.

This argument does not cross Re\Omega=0. The zero-frequency point and negative-imaginary-axis lateral limits/branch cut require a separate analysis. Matching radii and rays give the same continuation on connected overlaps by agreement in the first quadrant and the identity theorem.

## 4. Direct series and finite-contour truncation bounds

Under the practical monotonicity condition, put b_s=V(s)/|\Omega|. The exact normalized solution obeys
\[
|m(s)|\le e^{b_s},\qquad
m_x(s)=-\int_s^∞e^{2i\Omega(x(t)-x(s))}v(z(t))m(t)z'(t)dt,
\]
\[
|m_x(s)|\le|\Omega|(e^{b_s}-1).
\]
Truncating the Neumann series after terms n=0,...,N gives
\[
E_N(s)=\sum_{n>N}\frac{b_s^n}{n!}
\le e^{b_s}\frac{b_s^{N+1}}{(N+1)!},
\qquad |\delta m_x(s)|\le|\Omega|E_N(s).
\]
The derivative estimate follows because differentiating the first kernel removes exactly one factor 1/|\Omega|, and its lower-end contribution vanishes (K(z,z)=0).

If one also truncates the ray at S>s, put b_seg=[V(s)−V(S)]/|\Omega|. Let m_[N,S] be the corresponding finite-interval Neumann sum, with terminal value 1 and derivative 0. A combined bound is
\[
|m-m_{[N,S]}|\le e^{b_s}-e^{b_{seg}}+
\sum_{n>N}\frac{b_{seg}^n}{n!},
\]
with derivative error at most |\Omega| times the same expression. The first difference bounds all ordered integration regions with at least one variable beyond S. This is rigorous but has a slow algebraic tail, motivating the next accelerator.

## 5. Asymptotic-polynomial residual accelerator

Let p_N(z)=\sum_{n=0}^N a_nz^{-n}, a_0=1,a_{-1}=0, with
\[
a_{k+1}=\frac{[k(k+1)-L]a_k-k^2a_{k-1}}
{2i\Omega(k+1)}.
\]
Direct substitution, with no asymptotic limit, gives the exact residual
\[
\mathcal E p_N=R_N(z)=c_Nz^{-N-2}+d_Nz^{-N-3},
\]
where \mathcal E=fD_z²+(f'+2i\Omega)D_z−v and
\[
c_N=[N(N+1)-L]a_N-N^2a_{N-1},\qquad d_N=-(N+1)^2a_N.
\]
For e=m−p_N,
\[
e=T e-KR_N.
\]
Here T is the v-kernel integral above and KR_N denotes integration with the same K against R_N(t)dt. This sign is correct.

At z_R=r+iR, R>0, let
\[
F_R=\int_R^∞|R_N(r+is)|ds
\le\frac{|c_N|R^{-N-1}}{N+1}
+\frac{|d_N|R^{-N-2}}{N+2}.
\]
Let b_R=V(R)/|\Omega|. The residual Neumann series gives the simultaneous rigorous enclosures
\[
\boxed{|e(z_R)|\le E_R=\frac{F_R}{|\Omega|}e^{b_R},
\qquad |e_x(z_R)|\le F_Re^{b_R}=|\Omega|E_R.}
\]
The derivative estimate is slightly sharper than the also-valid F_R+V(R)E_R: differentiate the resolvent series termwise, use the unit exponential bound, and apply the same ordered-simplex estimate to the remaining v factors. This is valid for every finite N. It does not assume convergence of the usually divergent formal inverse-power series. The initial derivative supplied to a z-based solver must be p_N'(z_R), or, if using m_x, f(z_R)p_N'(z_R).

## 6. Rigorous propagation of the tail enclosure

Suppose an exact solution \widetilde m is propagated from terminal data p_N(z_R), (p_N)_x(z_R). Let E0,E1 bound the exact terminal value and x-derivative discrepancies. On the segment s∈[0,R], the error satisfies
\[
\delta m(s)=\delta m(R)-K(s,R)\delta m_x(R)
+\int_s^R K(s,t)v(z(t))\delta m(t)z'(t)dt.
\]
Set V_seg=V(0)−V(R), b_seg=V_seg/|\Omega|, A=E0+E1/|\Omega|. Then at the real matching point:
\[
|\delta m(0)|\le A e^{b_{seg}},
\qquad
|\delta m_x(0)|\le E1+|\Omega|A(e^{b_{seg}}-1).
\]
These are bounds for exact propagation of imperfect terminal data. Numerical integration, quadrature, and rounding errors must be enclosed separately before calling a computed result certified.

## 7. Error in the QNM determinant or incoming coefficient

Let h=e^{i\Omega x}u_H^in denote the normalized horizon solution evaluated at the same real matching point; let h_x be its x-derivative. Then
\[
D=W_x(u_H^{in},u_\infty^{out})
=2i\Omega h m+h m_x-h_xm.
\]
No exponentially large matching factors remain in this normalized expression. If h,h_x are exact and E_m,E_mx bound the outgoing errors,
\[
|\delta D|\le |2i\Omega h-h_x|E_m+|h|E_{mx}.
\]
If h is also approximate, writing tildes for computed values and E_h,E_hx for horizon errors yields
\[
\begin{split}
|D-\widetilde D|\le{}&|2i\Omega\widetilde h-\widetilde h_x|E_m
+|\widetilde h|E_{mx}
+|2i\Omega\widetilde m+\widetilde m_x|E_h
+|\widetilde m|E_{hx}\\
&+2|\Omega|E_hE_m+E_hE_{mx}+E_{hx}E_m.
\end{split}
\]
When the infinity basis has been normalized/continued from W(u_inf^in,u_inf^out)=2i\Omega, A_in=D/(2i\Omega), so its error is bounded by E_D/(2|\Omega|). A small computed D is not a certified root: zero counting or root enclosure additionally needs analyticity on a fixed region and a boundary test such as Rouché/argument principle, including all numerical errors.
