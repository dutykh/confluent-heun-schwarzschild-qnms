# Reproducing the coefficient enclosures and the nine radial zero counts

Use the unpacked project root as the working directory. The verifier uses Python integer and Fraction arithmetic, the directed algebraic operations in mpmath 1.3.0, and MPFI/MPFR for elementary functions. Its ctypes interface targets the Linux LP64 ABI with libmpfi.so.0, libmpfr.so.6 and their GMP dependency. The calculation was run with Python 3.12.14, MPFI 1.5.4 and MPFR 4.2.2. No numerical root finder is required to test a proposed disk. The canonical verifier runs full recomputations in a temporary directory. Individual low-level generation commands below write directly to their result paths and should be run in a working copy.

## Commands

```bash
# Verify the backend and run regression checks
python -m cert.elementary
python -m cert.test_transport
python -m cert.test_quadrature
python -m cert.test_isolate

# Verify strict zero counts and error ratios using the stored sample enclosures
python -m cert.reproduce --check-saved

# Regenerate the horizon states, finite evolutions, corrections, tails and counts
python -m cert.reproduce --workers 4

# Independently regenerate the two generic CHE examples
python cert/generic.py
```

The first stage of the full scalar run writes 27 point amplitudes and three same-frequency comparisons. The second stage bounds the amplitude on each outer frequency disk with rational inequalities and checks the strict Rouché comparison from three samples. Any failed comparison stops the run. At a fixed proposed center the entire calculation can be rerun separately, for example

```bash
python -m cert.campaign 00 center
python -m cert.campaign 00 plus
python -m cert.campaign 00 minus
python -m cert.check 00
```

## Interpreting outputs

A mode key is the pair of digits ℓ and n. All frequencies are Ω = 2Mω and all nine target radii are 10⁻²⁸. The independent continued-fraction and spectral calculations assign the conventional overtone labels. The mathematical statement counts multiplicity within each specified open disk. The boundary contains no zero.

- `campaign.py` fixes the exact rational centers, derivative-sample shift, contour choices and expansion orders.
- `results/XX_{center,plus,minus}.json` gives the complex enclosure of the exact amplitude at each frequency. Endpoints are rational strings. The diagnostics retain the finite estimate and tail contribution.
- `results/XX_disk.json` gives the outer analytic bound, error in the slope approximation and both sides of the strict root-isolation inequality. The `proved` flag follows from exact rational comparisons. Decimal hints aid reading.
- `results/00_fixed_endpoint_errors.json` bounds both actual absolute errors and their quotient at the same specified center frequency for three pairs (R,N). Their reference is an enclosure of the exact amplitude, which can differ from zero at the center.
- `results_generic.json` contains the normalized generic coefficient rectangles and the proof diagnostics for cases A and B. Exact endpoints define the enclosures. The elapsed time field depends on the host.
- `transport_notes.md`, `generic_notes.md` and `proof_root_isolation.md` supply detailed derivations and implementation notes. The main report also derives the analytic quadrature bounds.

The scalar certificates cover the physical continuation and nine disks specified in the report. Thresholds, branch-cut limits and zeros elsewhere require separate analysis.

## Fixed frequency neighborhoods and exceptional limits

To reproduce the full manuscript in one run, use the following command from the root of a fresh copy. The worker count controls the two independent collections of amplitude evaluations.

```bash
python -B -m cert.verify --full --workers 3
```

To recompute the rational zero counts, interpolation enclosures and error ratios and check the published exceptional-case inequalities from the saved endpoints:

```bash
python -B -m cert.verify --check-saved
```

The latter command uses the archived input intervals. The full run regenerates them from the equations. Additional modules can also be run individually.

| Command after `python -m` | Quantity reproduced |
| --- | --- |
| `cert.neighborhoods --workers 4` | The 72 additional point amplitudes, nine complex derivative disks, strict boundary lower bounds and annulus exclusions. |
| `cert.lateral` | Both continuations at −i/4 and the regularized continuations at −i/2, at the centers and on the two stated complex squares. |
| `cert.threshold` | Monopole Wronskian on the punctured lower half-disk of radius 1/1024, including both limiting values on the cut. |
| `cert.resonant_horizon` | Enclosed physical horizon data for ℓ=0,1,2 and the first three positive integer resonances in λ=2iΩ. |
| `cert.resonant` | Four accessory-parameter sign brackets and a common box for the resolved origin data. |

The proof of uniqueness and simplicity for the real accessory brackets uses the analytic Sturm identity in the manuscript. The static Legendre Wronskian for all multipoles is an analytic evaluation. Output files are named `results_neighborhoods.json`, `neighborhood_samples.json`, `results_lateral.json`, `results_threshold.json`, `results_resonant_horizon.json`, and `results_resonant.json`. The module `cert.test_neighborhoods` checks the interpolation inequalities, and `cert.test_lateral` repeats the two cut evaluations with a different spatial contour and truncation order.

The canonical saved-data check compares the recomputed disk JSON, fixed-endpoint
error quotients, and the nine derivative neighborhoods with the released
records. Differences raise an exception. It also checks the eight release
tables, the historical 68-grid spectral record, the accessory comparison
identity and the Schultze recurrence comparison. In full mode it reruns all
five transport, quadrature, isolation, neighborhood, and lateral regression
modules and compares the complete regenerated certificate JSON records after
discarding the three explicit wall-time fields. It recomputes the Schultze
comparison through its maximum recorded cutoff at both precisions.
