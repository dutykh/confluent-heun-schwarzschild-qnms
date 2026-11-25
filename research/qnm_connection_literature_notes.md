# Targeted literature audit: a computable QNM connection coefficient

Search date: 19 November 2026. This supplements `literature_notes.md`. Primary papers were opened wherever noted; reported theorems and numerical certificates have been inspected for scope, not independently reproduced. This is a targeted search, not a proof of exhaustive absence or priority.

## Main finding

The proposed core construction is well supported but substantially precedented: a horizon/infinity Wronskian, Volterra or Born integral equations, complex radial contours, analyticity on frequency patches, explicit endpoint error estimates, and even certified bounded QNM counts all occur in the literature. Two particularly close sources are Correia et al. (2025 preprint, 2026 publication) for Schwarzschild rotated Volterra integrals, and Chen et al. (November 2026 preprint) for endpoint Volterra estimates plus interval propagation and a certified extremal-Kerr count. They must be discussed explicitly before any novelty claim.

## Closest primary sources

### 1. Correia–Gopalka–Isabella–Wolz: Schwarzschild rotated Volterra connection

Miguel Correia, Tushar Gopalka, Giulia Isabella, Anna M. Wolz, *Analyticity of the black hole S-matrix*, JHEP **10 (2026) 150**, published 14 October 2026. Preprint v2 dated 26 November 2025.

- DOI: https://doi.org/10.1007/JHEP10(2026)150
- Primary PDF: https://arxiv.org/pdf/2511.11794
- Primary HTML: https://arxiv.org/html/2511.11794v2

PDF §2.3, (2.16)–(2.23): Volterra iteration, factorial bounds, uniform convergence and holomorphic Jost solutions. Equations (2.25), (2.27) give

\[
W[\phi_L,\phi_R]=2i\omega A_{\rm in},\qquad
A_{\rm in}=1-\frac1{2i\omega}\int_{-\infty}^{\infty}V\chi_L\,dx.
\]

Section 2.4, (2.38)–(2.40), gives rotated-contour factorial bounds and the continued integral for \(\operatorname{Im}(e^{i\theta}\omega)>0\). Section 4.1 applies this to scalar, electromagnetic and gravitational Schwarzschild potentials; singularities at \(x=i(2n+1)\pi R_s\) constrain the contour, with \(|\theta|<\pi/2\). Distinguishes the incoming coefficient/transmission from reflection and its different Stokes analyticity. This is direct prior art for the conceptual Schwarzschild construction. It does not present a certified numerical QNM root solver. **PDF numbering differs from HTML numbering.**

### 2. Chen–Zhang–Cao–Jing: rotated endpoint integrals and certified count

Changkai Chen, Xiaohua Zhang, Zhoujian Cao, Jiliang Jing, *Complete quasinormal modes of Type-D black holes. II. Extremal Kerr black holes*, arXiv **2610.39095v3**, 3 November 2026. Preprint; no publication claimed here.

- https://arxiv.org/abs/2610.39095
- https://arxiv.org/pdf/2610.39095
- https://arxiv.org/html/2610.39095v3

Appendix B is exceptionally close to the proposed numerical proof framework: endpoint ODE (B2); ray condition (B5); double Volterra integral and contraction (B6); noncritical endpoint determinants (B7); weighted residual contraction (B8); amplitude/derivative remainder bounds (B9); Taylor propagation bounds (B12). Parameters avoid \(\omega=0\) and \(2M\omega=m\), and the angular branch is separately isolated. Section IV.A reports an interval-certified count of exactly one simple damped mode in the conventional-sheet \((\ell,m)=(3,3)\) rectangle \(0.05<\Re\omega<1.48\), \(-3<\Im\omega<-0.03\), with \(M=1\). Appendix B controls continuous boundary arcs, not only samples. The larger mode survey uses convergence comparisons explicitly distinguished from certified frequency errors. Appendix C supplies local critical-frequency matching and exclusion bounds. Thus generic claims that “complex rays + Volterra errors + certified QNM count” are new would be untenable. **The underlying equation is double-confluent, and the certified domain is specific.**

### 3. Chen–Jing–Cao–Wang: direct incoming-amplitude QNM computation

Changkai Chen, Jiliang Jing, Zhoujian Cao, Mengjie Wang, *Complete quasinormal modes of type-D black holes*, Phys. Rev. D **112, 103036 (2025)**, published 19 November 2025.

- https://doi.org/10.1103/f8m8-vr4l
- https://arxiv.org/abs/2506.14635
- https://arxiv.org/html/2506.14635v5
- Author code/data: https://github.com/IronChen1/Complete-QNMs

Section II.4 uses analytic continuation and the vanishing incoming amplitude. Equations (49)–(55) give the asymptotic ratio condition, including the radial phase \(\arg r=2k\pi+\pi/2-\arg\omega\). Section II.5 supplies asymptotic amplitudes and renormalized angular momentum. The stated errors in §II.5.3, (69)–(71), are incoming/reflected amplitude ratios, normalized incoming amplitudes, and comparison with reference frequencies. These are useful residual/convergence diagnostics, not by themselves rigorous frequency enclosures. Their special-frequency discussion recognizes gamma/Pochhammer singularities in amplitude verification. Claims of complete spectra and efficiency must be read within their explored sectors and numerical checks. This is a direct computational competitor, including branch-cut continuation; do not treat it as merely a formal connection theorem.

### 4. Chen–Zhang–Cao–Jing: newly submitted RN continuation

*Complete quasinormal modes of Type-D black holes. III. Reissner–Nordstrom spectra from zero to extremal charge*, arXiv **2611.24647v1**, submitted 18 November 2026.

- https://arxiv.org/abs/2611.24647
- https://arxiv.org/html/2611.24647v1

Sections III.A–B factor the horizon phase, impose \(h(0)=1\), continue through ordinary Taylor disks, and extract the incoming infinity amplitude on specified radial/frequency sheets. A monodromy projector gives alternative horizon data away from integer index differences. Appendix A explicitly separates convergent horizon initialization and Taylor transport from generally divergent infinity asymptotic series: minimum-term truncation is described as a diagnostic, not a rigorous bound. Resonant horizon indices require a separate frequency-continued construction. It treats gravitational RN and the Schwarzschild limit, rather than arbitrary five-parameter CHE. This very recent preprint reinforces that a single incoming amplitude plus branch-aware continuation is already an active established approach.

### 5. Donninger–Schlag–Soffer: rigorous Schwarzschild Jost Volterra analysis

Roland Donninger, Wilhelm Schlag, Avy Soffer, *A proof of Price's Law on Schwarzschild black hole manifolds for all angular momenta*, Adv. Math. **226 (2011), 484–540**.

- https://doi.org/10.1016/j.aim.2010.06.026
- https://arxiv.org/abs/0908.4292
- https://arxiv.org/pdf/0908.4292

Section 3, Lemma 3.1 and (8), constructs Jost solutions using \(K=(e^{2iz(y-x)}-1)V(y)/(2iz)\). Equation (10) identifies the Wronskian with the incoming coefficient. Appendix B, Lemmas B.1–B.3, supplies existence, exponential bounds, differentiation and parameter dependence for Volterra equations. This is precise prior art for the upper-half-plane construction and differentiability. The purpose is scattering/decay, not lower-half-plane validated QNM computation. Section 6 discusses zero-frequency exceptions; Remark 6.2 records the scalar zero-frequency horizon solution through a Legendre polynomial. This source prevents attributing the preliminary UHP existence argument to the present project.

### 6. Leaver: convergent general spheroidal expansions and amplitudes

E. W. Leaver, *Solutions to a generalized spheroidal wave equation: Teukolsky's equations in general relativity, and the two-center problem in molecular quantum mechanics*, J. Math. Phys. **27 (1986), 1238–1265**.

- https://doi.org/10.1063/1.527130
- Author manuscript: https://edleaver.com/Misc/EdLeaver/Publications/SolutionsToGeneralizedSpheroidalWaveEquation.pdf

Section VI supplies bilateral Coulomb expansions, recurrence/minimality condition (112)–(113), convergence discussion (114)–(116), and explicit asymptotic phases/amplitudes (134)–(136). The auxiliary exponent must select the correct branch; half-integer roots are often spurious. Section VII gives confluent-hypergeometric expansions, with slower practical convergence. The introduction carefully distinguishes numerically tested Coulomb/Jaffé representations from untested second-solution formulas. An asymptotic term-ratio estimate is not automatically an interval-certified finite-tail bound. His 1985 QNM continued fraction and 1986 Schwarzschild Green-function paper remain essential baselines:

- https://doi.org/10.1098/rspa.1985.0119
- https://doi.org/10.1103/PhysRevD.34.384 (with 1988 erratum)

### 7. MST and exact conformal-block/Painlevé amplitudes

Mano–Suzuki–Takasugi, PTP **95 (1996), 1079–1096**, https://doi.org/10.1143/PTP.95.1079 and https://arxiv.org/abs/gr-qc/9603020; Mano–Takasugi, PTP **97 (1997), 213–232**, https://doi.org/10.1143/PTP.97.213 and https://arxiv.org/abs/gr-qc/9611014. Hypergeometric/Coulomb series and matching amplitudes are exact constructive representations with an auxiliary renormalized angular momentum. They already compute incident/reflected/transmitted amplitudes. Auxiliary-branch and exceptional-parameter handling matter.

Bonelli–Iossa–Panea Lichtig–Tanzini, PRD **105, 044047 (2022)**, https://arxiv.org/html/2105.04483v2: §4 connection matrices, §5.2 outgoing QNM condition (79)–(87); gamma factors combine with irregular conformal blocks/NS free energy, and accessory inversion remains necessary. The paper explicitly presents an instanton expansion; §4.3 numerical checks concern regular–regular matching. It is not a finite gamma-only formula in the original CHE coefficients or a certified truncation-error algorithm.

Carneiro da Cunha–Cavalcante, *Confluent conformal blocks and the Teukolsky master equation*, PRD **102, 105013 (2020)**, https://arxiv.org/abs/1906.10638: connection matrices through a Painlevé-V Riemann–Hilbert problem. This is an exact implicit/special-function framework, not an absence of a connection solution. The initial audit contains further Fredholm/recurrence references.

### 8. Complex scaling and older scattering theory

Ogawa–Hirose–Morikawa, *Complex scaling approach to quasinormal modes of Schwarzschild and Reissner–Nordström black holes*, PTEP **2026(9), 093E03**, https://doi.org/10.1093/ptep/ptag160 and https://arxiv.org/html/2604.20442v2. Section 3 rotates tortoise coordinates and discretizes in Gaussian bases. It notes branch restrictions and the practical Gaussian bound \(\theta<\pi/4\). Equations (71)–(73) quantify spread under parameter variation, explicitly a stability diagnostic rather than a rigorous finite-basis error. Thus different from a certified scalar connection evaluator, while defeating novelty of contour rotation itself.

Rakityansky–Sofianos–Amos, *A method of calculating the Jost function for analytic potentials*, Nuovo Cimento B **111 (1996), 363–380**, https://doi.org/10.1007/BF02724658 and https://arxiv.org/abs/nucl-th/9505022. Primary abstract verifies variable constants plus complex coordinate rotation for direct Jost functions, including resonance zeros and Coulomb tails. Older, general scattering precedent; not a Schwarzschild-specific certified solver.

Bachelot–Motet-Bachelot, *Les résonances d'un trou noir de Schwarzschild*, Ann. IHP Phys. Théor. **59 (1993), 3–68**, https://www.numdam.org/item/AIHPA_1993__59_1_3_0/ and https://www.numdam.org/item/AIHPA_1993__59_1_3_0.pdf. Primary abstract verifies rigorous Schwarzschild scattering, meromorphic continuation and numerical resonance work. Full theorem-level continuation assumptions were not audited in this pass.

## Novelty assessment and useful differentiation

These are judgments from comparison, not claims proved by a search.

| Proposed ingredient | Assessment |
|---|---|
| One horizon–infinity connection coefficient; zero equals QNM | Standard scattering formulation; Leaver/MST/CFT/Heun work already implements it |
| Wronskian at finite matching point | Standard; fixes normalization and avoids an ill-conditioned dominant-wave subtraction |
| UHP Volterra equation and factorial convergence | Established, including precise Schwarzschild analysis |
| Rotated contour to continue into damped-frequency domain | Established in general resonance theory and explicitly Schwarzschild |
| Computable endpoint residual, derivative and transport bounds | Explicit in November 2026 extremal-Kerr work; new specialization alone is insufficient |
| Rouché/argument-principle certified root counts | Standard complex analysis, with a direct recent black-hole implementation |
| Generic CHE coefficient algorithm with explicit uniform bounds and reproducible code | Potentially useful contribution, but must compare exact formulas, scope, sharpness, cost and exception handling with existing work |
| Uniform passage through resonant horizon indices / cut-bank limits | Possible substantive target if a precise missing theorem is identified and proved |
| Lower cost or stronger precision/domain certificate than existing solvers | Defensible only after measured comparisons with equivalent normalizations and domains |

A useful immediate deliverable is an independent, fully specified scalar-Schwarzschild benchmark with normalized endpoint solutions and transparent error accounting. Its scientific value does not require claiming a first connection solution. A publishable general-CHE advance would need a precise additional theorem or demonstrated numerical advantage beyond the above sources.

## Error language for the current project

Keep four levels separate: (1) an analytic identity; (2) a convergent approximation with mathematical remainder bounds; (3) floating-point implementation with precision/contour/refinement agreement; (4) a validated enclosure accounting for tail, transport, quadrature, rounding and parameter dependence. Small Wronskian residual alone neither encloses a root nor excludes nearby extra roots. A frequency error assertion needs derivative/nondegeneracy control or a validated contour argument. At resonant or critical parameters, multiply/divide normalization factors only after their zeros and poles are accounted for.

## Search limits and remaining audit items

The recent Ferreira–Castro paper https://arxiv.org/abs/2609.17727 was located and its primary abstract retrieved: it claims a canonical-CHE recurrence/regular-to-irregular limit defining one spectral Wronskian for weak-field Ernst scalar QNMs and quasibound states. Full primary text could not be opened in this pass, so detailed theorem/error claims are not relied upon. This is another close comparator for a future full paper.

The 2026 extremal-Kerr certificate deserves independent reproduction before using its numerical constants as trusted benchmarks. This audit confirms the content and stated scope of its primary source, not the correctness of every proof or computation. No broad “unsolved” or “first certified method” claim is justified by this search.
