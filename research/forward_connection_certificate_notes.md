# Forward connection-coefficient certificate: targeted comparison

19 November 2026. This supplements `generic_che_literature_notes.md` and the independent branch audit `adjoint_coefficient_prior_notes.md`.

## Candidate under investigation

For L=D²+pD+V, let q solve the formal adjoint equation

\[
q''-pq'+(V-p')q=0.
\]

If P_N is a finite formal infinity expansion and R_N=LP_N, the bilinear concomitant

\[
B_N=(pP_N+P_N')q-P_Nq'
\]

satisfies B_N′=qR_N. Under the required endpoint limits, the coefficient b is therefore

\[
b=\frac{B_N(R)}{\alpha}+\frac1\alpha\int_R^\infty q(t)R_N(t)\,dt.
\]

The specific proposed improvement is a computable bound on this coefficient tail from the already propagated regular/horizon state, without first enclosing an exact outgoing solution. The identity is an elementary instance of Green/Lagrange duality. The mathematical distinction, if any, is in explicit CHE constants, conditioning, cost, and a useful uniform parameter domain.

## Direct forward Jost methods are old

S. A. Rakityansky, S. A. Sofianos, K. Amos, **A method for calculating the Jost function for analytic potentials**, *Nuovo Cimento B* 111 (1996), 363–378, DOI [10.1007/BF02724658](https://doi.org/10.1007/BF02724658), [arXiv:nucl-th/9505022](https://arxiv.org/pdf/nucl-th/9505022).

Equations (12)–(19) use variable constants to evolve finite-radius incoming/outgoing coefficients of the regular solution on a rotated complex ray. Equation (19) expresses the incoming coefficient as an integral of an outgoing reference wave times the potential times the regular solution. Equation (24) takes its infinity limit to obtain the Jost function; the relevant sector is Im(kr)>0. Thus direct forward computation of one amplitude, with no outgoing solution integrated backward, is an established technique. The checked article demonstrates empirical accuracy rather than the proposed explicit high-order residual tail certificate. Section IV treats Coulomb tails using Coulomb reference waves.

S. A. Sofianos, S. A. Rakityansky, S. E. Massen, **Jost function for singular potentials**, *Phys. Rev. A* 60 (1999), 337, [arXiv:nucl-th/9901023](https://arxiv.org/pdf/nucl-th/9901023), DOI [10.1103/PhysRevA.60.337](https://doi.org/10.1103/PhysRevA.60.337). Equations (5)–(10) state the same forward coefficient equations, limiting Jost function, and complex rotation explicitly. It extends origin initialization to repulsive singular potentials. These papers are close method-level precedents, though their physical endpoint differs from a black-hole horizon.

## Existing black-hole forward integrations

K. Glampedakis, N. Andersson, **Quick and dirty methods for studying black-hole resonances**, [arXiv:gr-qc/0304030](https://arxiv.org/pdf/gr-qc/0304030), §IIIA. A horizon-initialized logarithmic derivative is evolved forward, then changed to a Prüfer phase and continued to infinity. Equations (24)–(30) extract the scattering ratio from the limiting phase; complex contours handle QNMs. This already avoids an outgoing backward solve in one of its two schemes. It uses a nonlinear phase and balanced-wave contour, not the candidate's scaled linear state or rigorous coefficient-tail bound. Numerical validation is empirical.

Correia et al., arXiv:2511.11794v2, equations (2.27) and (2.40), already give the incoming amplitude as a regular-solution-weighted potential integral, including rotated contours; see `qnm_connection_literature_notes.md`. Replacing a leading reference wave by a higher-order asymptotic test is the natural residual generalization of that identity.

## Important comparison with Chen/Motygin matching

C. Chen, Z. Cao, J. Jing, **Efficient and stable computation of gravitational-wave fluxes from generic Kerr orbits via a unified Heun-function framework**, [arXiv:2605.09250v4](https://arxiv.org/pdf/2605.09250), §IIA, pp.4–5.

Equation (18) chooses the far-field radius by a least-term heuristic. The algorithm then transports horizon and asymptotic solutions toward a matching point and solves a 2×2 system. However, step (i) moves the matching point all the way to the far-field radius for demanding parameters. In that configuration the backward chain has zero length: a finite asymptotic series is matched directly against the forward horizon state. Therefore “no outgoing ODE solve” alone is not a robust distinction even from this implementation. The candidate's stronger distinction is the direct coefficient error bound based on the measured horizon state, instead of an outgoing-state enclosure followed by projection. Chen et al.'s reported errors here are numerical comparisons, not interval certificates.

## Output residual identities and bounds

Kievsky et al., [arXiv:1002.3756](https://arxiv.org/pdf/1002.3756), equations (4)–(7), explicitly derive scattering amplitudes as integrals of the physical solution against differential residuals of reference functions. Kato's 1951 **Upper and Lower Bounds of Scattering Phases**, DOI [10.1143/ptp/6.3.394](https://doi.org/10.1143/ptp/6.3.394), provides rigorous output bounds from weighted residual norms in real self-adjoint scattering. Neither checked result supplies the proposed complex-ray CHE certificate. See `adjoint_coefficient_prior_notes.md` for exact scopes and additional sources.

Olver (1997) supplies the broad rank-one residual/derivative-error machinery. The extremal-Kerr Chen et al. paper arXiv:2610.39095v3 Appendix B already implements explicit infinity residual certificates and finite continuation. A fair numerical comparison must propagate its outgoing-state error into the same coefficient before claiming the new bound is sharper.

## New rigorous Stokes paper: relevant but different scope

M. Loday-Richaud, M. Mezzarobba, P. Remy, **Rigorous numerical computation of the Stokes multipliers for linear differential equations with single level one**, [arXiv:2601.04901v1](https://arxiv.org/pdf/2601.04901), January 2026.

This gives implemented arbitrary-precision local Stokes matrices with rigorous errors, through Borel-plane regular-singular connection problems. The introduction p.4 explicitly distinguishes transporting these local automorphisms to a common global basis, which remains a connection problem requiring effective summation. It does not directly subsume the candidate's regular-to-irregular one-entry certificate; it does rule out broad claims that rigorous irregular-singular computation itself is new.

## Defensible assessment

No exact published equivalent of the particular formula

\[
\frac{F_R}{|\alpha|}\left(|q(R)|+\frac{|q'(R)|}{\sigma}\right)e^{W_R/\sigma}
\]

was located in the inspected CHE/QNM sources. This bounded search result is not proof of novelty. The formula should be derived and audited independently, including the definition of F_R, the ray parametrization, the strict damping margin σ, and interval errors in q(R),q′(R).

A defensible possible contribution is an explicit high-order CHE coefficient-tail theorem whose evaluated bounds are materially sharper or cheaper near QNMs than projected outgoing-state bounds, with a stable forward implementation and independent certified validation. The lack of an outgoing backsolve, the coefficient limit, the residual identity, and Grönwall estimation are individually established ideas. Whether their CHE specialization is publishable depends on the quantitative result rather than on rebranding those principles.
