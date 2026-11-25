# Generic CHE: a canonical infinity solution and one connection amplitude

This is an independently derived Volterra/reduction-of-order construction on explicit parameter charts. It is not a priority claim, a finite Gamma evaluation, or a global solution of all Stokes connection data.

## 1. Exact reduction

Start with
\[
y''+P y'+Qy=0,\quad P=\alpha+\frac{\beta+1}{z}+\frac{\gamma+1}{z-1},\quad Q=\frac{Az+B}{z(z-1)},\quad \alpha\ne0.
\]
Use the canonical A,B from the main report. Put
\[
\rho=-A/\alpha,\qquad y_{\rm alg}=z^\rho m,
\]
\[
u=\beta+1+2\rho,\quad v=\gamma+1,\quad C=\rho(\rho+\beta),\quad D=A+B+\rho(\gamma+1).
\]
Then exactly
\[
\mathcal E m=m''+\left(\alpha+\frac uz+\frac v{z-1}\right)m'
+\left(\frac C{z^2}+\frac D{z(z-1)}\right)m=0.
\]
Here the symbol u in the coefficient dictionary is u=beta+1+2rho (not the original dependent variable). The residual potential V=C/z²+D/[z(z−1)] is O(z^−2). The cancellation of its O(z^−1) term uses alpha*rho=−A.

Fix branches and define
\[
J(z)=e^{\alpha z}z^u(z-1)^v,\qquad
K(z,t)=J(t)\int_z^tJ(\zeta)^{-1}\,d\zeta,
\quad (Tg)(z)=\int_z^{\infty_c}K(z,t)V(t)g(t)\,dt.
\]
Then
\[
\boxed{m=1-Tm=\sum_{n=0}^\infty(-T)^n1},\qquad
m'(z)=\int_z^{\infty_c}\frac{J(t)}{J(z)}V(t)m(t)\,dt.
\]
The signs follow from (Jm')'=−JVm. Also K(z,z)=0 and ∂_zK(z,t)=−J(t)/J(z).

**Damping sign:** the appropriate recessive rays for the algebraic branch satisfy Re(alpha*c)<0. On such rays e^(−alpha*z) grows, so a normalized algebraic solution is distinguished. Choosing Re(alpha*c)>0 makes the other branch exponentially small and m→1 alone nonunique.

## 2. Uniform parameter charts covering every alpha≠0

For complete angular coverage use a radial ray z=c r, |c|=1, r≥R0>1, instead of requiring its start to be real. Fix c on a parameter chart Re(alpha*c)<0. On a compact subset let
\[
\sigma_0=\inf[-\Re(\alpha c)]>0,\quad U_0=\sup|u|,\quad V_0=\sup|v|.
\]
Choose R0 sufficiently large that
\[
\boxed{\sigma=\sigma_0-U_0/R_0-V_0/(R_0-1)>0.}
\]
The logarithms may be chosen as log(z)=log(c)+log(r) and log(z−1)=log(c)+Log(r−c^−1). Since Re(r−c^−1)>0, these are consistent and log(z−1)−log(z)→0.

Along the ray,
\[
\frac d{dr}\log|J(cr)|
=\Re\left(\alpha c+\frac ur+\frac v{r-c^{-1}}\right)\le-\sigma.
\]
Consequently, for t≥s≥R0,
\[
|J(ct)/J(cs)|\le e^{-\sigma(t-s)},\qquad |K(cs,ct)|\le1/\sigma.
\]
Every alpha0≠0 has such a chart: choose c=−conj(alpha0)/|alpha0| and then shrink the parameter neighborhood. No nonresonance condition on beta or gamma is needed for this infinity construction. The alpha=0 confluence is excluded.

A potential tail majorant is
\[
\mathcal B(r)=\frac{|C|}{r}+|D|\log\frac r{r-1}.
\]
An often sharper alternative, using V=(C+D)/z²+D/[z²(z−1)], is
\[
\mathcal B_*(r)=\frac{|C+D|}{r}
+|D|\left(\log\frac r{r-1}-\frac1r\right).
\]
Either bounds the absolute V integral along the radial tail. Use parameter suprema for a compact-uniform bound.

With b=mathcal B(r)/sigma, ordered-simplex estimates give
\[
|T^n1|\le b^n/n!,\quad |m|\le e^b,\quad
|m'|\le\sigma(e^b-1).
\]
Parameter-local uniform convergence proves holomorphy on these charts. The series also solves the differential equation. Its successive inverse-power remainder bounds below establish the normalized algebraic asymptotic expansion. A second solution with the same full normalized recessive behavior would have bounded derivatives; its Wronskian with the first cannot equal a nonzero constant times J^−1, which grows exponentially. Therefore the normalization is unique in this recessive class.

For rays whose recessive angular domains overlap, contour deformation inside a fixed branch domain identifies the solutions: J(t)/J(s) decays and V=O(t^−2), so large arc contributions vanish, first for each iterate and then for the uniformly convergent sum. Across a Stokes boundary this argument does not apply without continuation data. The theorem defines canonical solutions on overlapping charts with specified paths; it does not compute all Stokes multipliers.

## 3. Real-start rays, vertical specialization, and sharper decay estimates

For z=a+c s, a>1, |c|=1, Re c≥0, s≥0, the simpler sufficient condition is
\[
\sigma=-\Re(\alpha c)-|u|/a-|v|/(a-1)>0.
\]
It gives the same exponential-ratio and kernel bounds. These real-start outward rays do not cover positive real alpha; the radial construction above does.

On a vertical ray c=i the exact logarithmic derivative is
\[
\frac d{ds}\log|J(a+is)|=-\Im\alpha
+\frac{(\Re u)s-a\Im u}{a^2+s^2}
+\frac{(\Re v)s-(a-1)\Im v}{(a-1)^2+s^2}.
\]
A sharper sufficient constant is
\[
\sigma_{\rm vert}=\Im\alpha-
\frac{|u|-\Im u}{2a}-\frac{|v|-\Im v}{2(a-1)}>0.
\]
It follows by maximizing each rational term over the full real line, an upper bound on its maximum for s≥0. One may improve it by proving the actual derivative is ≤−sigma on the required tail, but sampled positivity alone is not a proof uniform in s or in a frequency/parameter patch.

For the vertical tail s≥R>0, |a+is|≥s and |a−1+is|≥s, so
\[
\int_R^\infty|V(a+is)|ds\le (|C|+|D|)/R.
\]
A less coarse explicit bound integrates |C|/(a²+s²)+|D|/[sqrt(a²+s²)sqrt((a−1)²+s²)]; the first term is elementary and the second may be bounded using (a−1)²+s².

## 4. Finite inverse-power polynomial: exact recurrence and residual

Let P_N=∑_(n=0)^N a_n z^−n, a0=1,a_−1=0. Multiplication by z²(z−1) gives
\[
\boxed{a_{k+1}=\frac{[k(k+1-u-v+\alpha)+C+D]a_k
-[(k-1)(k-u)+C]a_{k-1}}{\alpha(k+1)}.}
\]
Define
\[
E_N=[N(N+1-u-v+\alpha)+C+D]a_N
-[(N-1)(N-u)+C]a_{N-1},
\]
\[
F_N=-[N(N+1-u)+C]a_N.
\]
Then the exact rational residual is
\[
\boxed{\mathcal EP_N=\frac{E_Nz+F_N}{z^{N+2}(z-1)}.}
\]
The identity was independently checked in SymPy for N=0,1,2 with generic independent alpha,u,v,C,D. There are no resonant denominator exceptions besides alpha=0; denominators are alpha(k+1).

On a vertical tail a+is, s≥R>0, its absolute integral is bounded by
\[
\boxed{\mathcal F_R=\frac{|E_N|R^{-N-1}}{N+1}
+\frac{|F_N|R^{-N-2}}{N+2}.}
\]
This bound uses the vertical geometry; do not apply it unchanged to an arbitrary radial ray.

On a radial tail z=c r, r≥R>1, one valid bound is
\[
\mathcal F_R\le\frac R{R-1}\left[
\frac{|E_N|R^{-N-1}}{N+1}
+\frac{|F_N|R^{-N-2}}{N+2}\right].
\]
A sharper alternative is
\[
\frac{|E_N|R^{-N-1}}{N+1}+|E_N+F_N|\,I_N(R),
\quad I_N(R)=\int_R^\infty\frac{r^{-N-2}}{r-1}dr
=-\log(1-R^{-1})-\sum_{k=1}^{N+1}\frac{R^{-k}}k.
\]
For numerical evaluation of I_N avoid subtracting nearly equal quantities; its positive integral or positive tail series is preferable.

## 5. Explicit truncation enclosures and finite quadratures

Since P_N→1 and P_N'→0,
\[
e=m-P_N=-Te-K(\mathcal EP_N),
\]
where KR denotes integration against K. If the J-ratio bound is at most one, |K|≤1/sigma, and mathcal B_R bounds the V tail, then
\[
\boxed{|e(z_R)|\le\epsilon_R=
\frac{\mathcal F_R}{\sigma}e^{\mathcal B_R/\sigma},
\qquad |e'(z_R)|\le\sigma\epsilon_R.}
\]
The derivative bound follows by differentiating the first kernel in the resolvent series. These are finite-polynomial residual bounds; convergence of the infinite formal asymptotic series is not assumed.

On a finite ray segment from a chosen base z0 to zR, let T_R integrate only to zR and define
\[
q_R(z)=P_N(z_R)-K(z,z_R)P_N'(z_R),\qquad
\widetilde m_K=\sum_{k=0}^K(-T_R)^kq_R.
\]
Thus only finite quadratures of explicit kernels and coefficients occur. With
\[
b=\mathcal B_{\rm seg}/\sigma,\quad
Q=|P_N(z_R)|+|P_N'(z_R)|/\sigma,\quad
\tau_K=Qe^b b^{K+1}/(K+1)!,
\]
finite-iteration errors are ≤tau_K and ≤sigma*tau_K. Terminal residual errors propagated to z0 give the total analytic bounds
\[
E_m\le2\epsilon_Re^b+\tau_K,\qquad
E_{m'}\le\sigma\epsilon_R(2e^b-1)+\sigma\tau_K.
\]
Here primes always mean z derivatives; unit-speed ray derivatives have equal moduli. If the finite segment is outside the region where the stated sigma bound holds, use a separately controlled ordinary-point propagation, rather than reusing these estimates.

For y_alg=z^rho m,
\[
E_y\le|z^\rho|E_m,\qquad
E_{y'}\le|z^\rho|\left(E_{m'}+|\rho/z|E_m\right).
\]
Numerical quadrature, finite-precision arithmetic, and any subsequent ordinary-point transport require additional bounds.

## 6. A generic exponential connection coefficient without fixing a dominant solution

Let
\[
\mathcal W(z)=e^{-\alpha z}z^{-\beta-1}(z-1)^{-\gamma-1}.
\]
On a sufficiently distant ray tail, m is nonzero (enlarge the starting radius until e^(mathcal B/sigma)−1<1). Define by reduction of order
\[
y_{\rm exp}(z)=-\alpha y_{\rm alg}(z)
\int_{z_*}^z\frac{\mathcal W(t)}{y_{\rm alg}(t)^2}\,dt.
\]
Then exactly W(y_alg,y_exp)=−alpha*mathcal W. Moreover
\[
\frac{\mathcal W}{y_{\rm alg}^2}\sim
 e^{-\alpha z}z^{-2\rho-\beta-\gamma-2},
\]
and growing-direction integration yields
\[
y_{\rm exp}\sim e^{-\alpha z}z^{-\rho-\beta-\gamma-2}
=e^{-\alpha z}z^{A/\alpha-\beta-\gamma-2},
\]
with leading coefficient one. The coefficient −alpha in the reduction formula is essential for this normalization. The leading primitive follows by integration by parts (or the standard Laplace-endpoint estimate), justified by the strictly positive exponential growth of e^(−alpha*z) and the differentiable inverse-power estimates already obtained.

Changing the integration constant adds a multiple of y_alg. It therefore changes the algebraic coefficient in a decomposition but leaves the exponential coefficient unchanged. For any solution h of the original CHE,
\[
\boxed{\mathcal B[h]=\frac{e^{\alpha z}z^{\beta+1}(z-1)^{\gamma+1}}
{\alpha}\,W(h,y_{\rm alg}).}
\]
This is constant in z and is precisely the coefficient of the unit-leading growing exponential branch on the stated recessive chart. The construction establishes one meaningful generic connection amplitude without falsely asserting that a dominant asymptotic expansion uniquely determines its solution. Its parameter holomorphy additionally requires h to be chosen holomorphically; regular-singularity resonance may obstruct naive h=HeunC normalization and must be handled separately.

For approximate local h,h' and y,y', the numerator W=h y'−h' y has the transparent error bound
\[
E_W\le|\widetilde h|E_{y'}+|\widetilde y'|E_h+
|\widetilde h'|E_y+|\widetilde y|E_{h'}+E_hE_{y'}+E_{h'}E_y.
\]
Multiply by the modulus of the explicit prefactor to bound mathcal B[h]. No statement about zeros or root counts follows without complete arithmetic/transport bounds and a contour argument.

## 7. Schwarzschild normalization check and scope

For the scalar Schwarzschild dictionary alpha=2iOmega,beta=0,gamma=−2iOmega,delta=2Omega²,eta=−L, one finds
rho=−1+2iOmega and exponential power −1. Under the original gauge
u=e^(iOmega*z) z(z−1)^(−iOmega) H,
the unit-leading algebraic H branch becomes the unit outgoing u branch, and the unit-leading exponential H branch becomes the unit incoming u branch. The horizon solution H_H in the report has an additional physical normalization e^(−2iOmega). Therefore
\[
\boxed{A_{\rm in}=e^{-2i\Omega}\mathcal B[H_H].}
\]
This was checked directly by matching both asymptotic prefactors.

The generic construction allows arbitrary complex beta,gamma,delta,eta on compact alpha≠0 charts; it is therefore broader than the physical scalar family and Konno's factorization locus. It does not provide finite Gamma-only values, all Stokes matrices, an alpha=0 uniform limit, or automatically a regularized local Heun basis at resonance. Those remain distinct tasks. The exponential amplitude is invariant under the dominant-solution ambiguity, but it is still path/sector dependent through the selected canonical y_alg.

The other recessive infinity branch can also be constructed by the same theorem: substitute y=e^(−alpha*z)\widetilde y. The canonical parameters become (−alpha,beta,gamma,delta,eta), and the algebraic power for \widetilde y is A/alpha−beta−gamma−2. Its recessive rays satisfy Re(alpha*c)>0. Comparing these constructions after specified analytic continuations is the next step toward Stokes data; the present coefficient formula does not identify those transitions automatically.

## 8. Forward-only coefficient evaluation: Lagrange identity and tail enclosure

There is no need to propagate the infinity solution backward if only its exponential connection coefficient is required. Let h denote the selected source solution of the original CHE, H=z^(-rho)h, and
\[
q=JH=Jz^{-\rho}h,\qquad p=J'/J=\alpha+u/z+v/(z-1).
\]
A direct gauge transformation gives the adjoint equation
\[
q''-p q'+Wq=0,\qquad W=V-p'
=\frac{C+u}{z^2}+\frac D{z(z-1)}+\frac v{(z-1)^2}.
\]
This has a stable free forward derivative factor J(z)/J(t) on the same recessive ray. Define
\[
L(z,t)=J(t)^{-1}\int_t^zJ(s)ds.
\]
The forward integral equation is
\[
q(z)=q_R+L(z,z_R)q_R'
-\int_{z_R}^zL(z,t)W(t)q(t)dt.
\]
Note the minus sign in front of the W integral. Under the J-ratio bound above, |L|≤1/sigma. Consequently
\[
\sup_{\text{tail}}|q|\le Q_R\exp(\mathcal W_R/\sigma),
\qquad Q_R=|q_R|+|q_R'|/\sigma,
\quad \mathcal W_R\ge\int_{z_R}^{\infty_c}|W(t)||dt|.
\]
Dominated convergence in the derivative equation gives q'→0. Since L(z,t) has a finite limit as z→infinity, dominated convergence in the forward integral equation proves that q→q_infinity. The exact Wronskian expression in §6 then gives mathcal B[h]=q_infinity: p→alpha, m→1,m'→0.

For the finite polynomial P_N define the computable scalar
\[
\mathcal B_N(z_R)=\frac{[pP_N+P_N']q-P_Nq'}{\alpha}\bigg|_{z_R}.
\]
The Lagrange identity is exact:
\[
\frac d{dz}\left([pP_N+P_N']q-P_Nq'\right)=q\mathcal E P_N=qR_N.
\]
Because P_N→1 and P_N'→0,
\[
\boxed{\mathcal B[h]-\mathcal B_N(z_R)=
\frac1\alpha\int_{z_R}^{\infty_c}q(t)R_N(t)dt,}
\]
\[
\boxed{|\mathcal B[h]-\mathcal B_N(z_R)|
\le\frac{\mathcal F_R}{|\alpha|}
\left(|q_R|+\frac{|q_R'|}{\sigma}\right)
\exp(\mathcal W_R/\sigma).}
\]
Thus one forward source propagation and a finite polynomial suffice; no outgoing-state propagation is required. The bound is dependent on the actual coefficient state and can improve near a zero, but its numerical use requires errors in that state to be included.

Explicit forward-potential bounds are
\[
\mathcal W_R\le\frac{|C+u|}{R}
+|D|\log\frac R{R-1}+\frac{|v|}{R-1}
\]
on a radial tail cr, r≥R>1, and
\[
\mathcal W(s)\le\frac{|C+u|}{a}\left(\frac\pi2-\arctan\frac sa\right)
+\frac{|D|+|v|}{a-1}\left(\frac\pi2-\arctan\frac s{a-1}\right)
\]
on the vertical tail a+it, t≥s.

If terminal approximations have errors E_q,E_q', the finite expression has additional error
\[
\frac{|pP_N+P_N'|E_q+|P_N|E_{q'}}{|\alpha|}.
\]
Use |q_R|≤|qtilde_R|+E_q and |q_R'|≤|qtilde_R'|+E_q' in the residual tail bound. Errors in P_N and its coefficients also need inclusion if evaluated in finite arithmetic.

### Scalar Schwarzschild specialization

For the physically scaled horizon state h=e^(iOmega*x)u_H, h_xx−2iOmega*h_x−U h=0. Its forward kernel is
\[
K_f(s,t)=\frac{e^{2i\Omega[x(s)-x(t)]}-1}{2i\Omega},\qquad
h(s)=h_R+K_f(s,R)h_{x,R}+\int_R^sK_f(s,t)v(t)h(t)dt.
\]
With |K_f|≤k and the potential tail V_R,
\[
\sup|h|\le(|h_R|+k|h_{x,R}|)e^{kV_R}.
\]
The same dominated-convergence argument proves h_x→0 and h→A_in. If P_N and its residual are the scalar polynomials from the Schwarzschild note, define
\[
A_N(z_R)=\frac{(2i\Omega h-h_x)P_N+h fP_N'}{2i\Omega}.
\]
Then
\[
\frac d{dz}\left[(2i\Omega h-h_x)P_N+h fP_N'\right]=hR_N,
\]
\[
\boxed{A_{in}-A_N(z_R)=\frac1{2i\Omega}\int_{z_R}^{\infty}hR_Ndt,
\qquad
|A_{in}-A_N|\le\frac{F_R}{2|\Omega|}
(|h_R|+k|h_{x,R}|)e^{kV_R}.}
\]
Terminal-state error adds
[|2iOmega*P_N+fP_N'|E_h+|P_N|E_hx]/(2|Omega|), and the state magnitudes in the tail bound must also be enlarged by their errors. In particular, a tiny computed h near a root is not itself a certified coefficient bound. The method's stability or cost advantage requires a controlled numerical comparison, although its free forward exponential decays rather than grows on the chosen ray.
