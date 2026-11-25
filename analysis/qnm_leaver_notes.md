# Independent scalar Schwarzschild continued-fraction benchmark

The executable `qnm_leaver.py` derives and checks its recurrence from the radial ODE; it neither reads nor uses the proposed contour determinant. Output is `qnm_leaver_results.json`.

Conventions: time dependence exp(−iωt), Ω=2Mω, z=r/(2M), x=z+log(z−1), scalar ℓ=0,1,2. Roots have positive real part and negative imaginary part. Overtone labels n=0,1,2 use the corresponding continued-fraction inversion and ordering by damping; this calculation does not prove completeness of the spectrum.

With t=(z−1)/z, u=exp(iΩz)t^(−iΩ)(1−t)^(−iΩ)Σaₙtⁿ, the exact coefficients are αₙ=(n+1)(n+1−2iΩ), βₙ=−2n²+(8iΩ−2)n+8Ω²+4iΩ−ℓ(ℓ+1)−1, γₙ=(n−2iΩ)². Exact gauge substitution and all recurrence coefficients are checked by SymPy. The analytic derivative of the inverted fraction is checked independently with mpmath differentiation.

The finite continued fraction uses a_(N+1)=0, with no asymptotic-tail acceleration. The depth sequence is 50,100,200,400,800,..., up to the first depth at least 800 whose root differs from the previous depth by less than 1e−19. The deepest required solve was N=25600. Roots were solved at 60 decimal digits and repeated at 85 digits at the same final depth. The final root was also evaluated, without re-solving, at twice that depth.

| ℓ | n | Ω (rounded to15 significant digits) | Final depth | Last depth change | Estimated frequency correction at double depth |
|---:|---:|---|---:|---:|---:|
| 0 | 0 | 0.220909878160839 − 0.209791434173762 i | 3200 | 6.72969412099e-23 | 9.21455706062e-32 |
| 0 | 1 | 0.172233836672798 − 0.696104893612921 i | 12800 | 1.48486711026e-21 | 6.7841580664e-30 |
| 0 | 2 | 0.151483871070352 − 1.20215718007161 i | 25600 | 1.31252143824e-20 | 1.51468272435e-28 |
| 1 | 0 | 0.585872266534565 − 0.195319977827156 i | 800 | 6.53274858178e-24 | 2.8564315968e-33 |
| 1 | 1 | 0.528897301209665 − 0.612514783118094 i | 1600 | 9.28312186629e-23 | 8.61718852222e-32 |
| 1 | 2 | 0.459078669862603 − 1.08026685003821 i | 3200 | 4.03189447981e-22 | 6.32740512158e-31 |
| 2 | 0 | 0.967287744421426 − 0.193517551956576 i | 800 | 4.21282394851e-32 | 6.03920105408e-45 |
| 2 | 1 | 0.927701158039531 − 0.591207873975925 i | 800 | 2.16906805249e-25 | 1.1087486664e-35 |
| 2 | 2 | 0.861088108753152 − 1.01711680430855 i | 800 | 7.34132660164e-20 | 5.75329877023e-28 |

The double-depth correction is |F_(2N)(Ω_N)/F′_(2N)(Ω_N)|, a local Newton diagnostic, **not a certified truncation bound**. All last depth changes are below 1e−19; the largest double-depth frequency-correction diagnostic is below 5.8e−28. The largest 60-versus85-digit root difference is below 1.4e−61 after the final extension. High-precision finite-fraction residuals and individual arithmetic comparisons are retained in the JSON. Residuals within a fixed truncated equation alone do not measure infinite-fraction truncation error.

These numbers furnish an independent empirical baseline for the proposed connection-coefficient solver. Agreement would support implementation correctness; it would not by itself prove the latter solver’s continuation prescription or provide rigorous error bounds.

Reproduce with `python3 analysis/qnm_leaver.py --output analysis/qnm_leaver_results.json` (mpmath and sympy required).
