# Schwarzschild QNM connection coefficient

Research snapshot, 19 November 2026. Prepared for Dr. Davide Batic.

The accompanying standalone report derives the single scalar incoming amplitude
for time dependence `exp(-i omega t)` and dimensionless frequency
`Omega = 2 M omega`. The physical scalar equation is
`u_xx + (Omega^2 - U) u = 0`, with `x = z + log(z-1)`,
`z = r/(2M)`, and `U = (1-1/z)[ell(ell+1)/z^2 + 1/z^3]`.

The outgoing contour evaluator continues the real scattering solution into
`Re(Omega)>0`, using the upper complex-radius deformation specified in the report.
It does not provide negative-imaginary-axis lateral values or zero-frequency
limits. The regularized horizon module can evaluate its endpoint vector on the
axis, but that alone is not a complete axis QNM calculation.

## Main files

- `reports/Confluent_Heun_QNM_Coefficient.tex`: standalone mathematical report,
  including bibliography, proofs, exact normalization, truncation bounds,
  literature comparison, and numerical results.
- `reports/Confluent_Heun_QNM_Coefficient.pdf`: compiled report.
- `analysis/qnm_contour.py`: SciPy/DOP853 incoming-amplitude evaluator and
  nine-mode campaign; computes analytic infinity-tail estimates.
- `analysis/qnm_contour_mp.py`: independent mpmath Taylor contour evaluator.
- `analysis/qnm_contour_mp_overtone_run.py`: higher-overtone multiprecision check.
- `analysis/qnm_leaver.py`: independently derived inverted continued fraction.
- `analysis/qnm_spectral.py`: independent Chebyshev spectral campaign; imports
  `qnm_spectral_chebcoef.py` and `qnm_spectral_refine.py`.
- `analysis/qnm_horizon_regularized.py`: division-free reciprocal-Gamma horizon
  evaluator, including exact and near-resonant frequencies.
- `analysis/qnm_scattering_checks.py`: real-frequency flux normalization check.
- `analysis/qnm_compare.py`: cross-method comparison of saved numerical records.
- `analysis/*results.json`, `analysis/*notes.md`: full settings, numerical
  histories, method details, and limitations.
- `research/`: detailed primary-source literature and exceptional-frequency
  audit notes.
- `verification/scalar_resonance_verify.py`: symbolic and numerical identities.
- `reports/Historical_environment.json`, `requirements.txt`, `MANIFEST_SHA256.txt`: environment and
  file provenance.

## Reproduction

Use Python 3.12 or a compatible Python 3 environment with the listed packages.
From this directory:

```bash
python -m pip install -r requirements.txt
python analysis/qnm_contour.py --campaign
python analysis/qnm_leaver.py
python analysis/qnm_contour_mp.py --campaign
python analysis/qnm_contour_mp_overtone_run.py
che_replay=$(mktemp -d)
OPENBLAS_NUM_THREADS=1 python -B analysis/qnm_spectral_release.py --recompute "$che_replay/spectral.json"
python -B analysis/qnm_spectral_release.py --check-result "$che_replay/spectral.json"
python analysis/qnm_horizon_regularized.py --test
python analysis/qnm_scattering_checks.py
python -B analysis/qnm_compare.py --spectral "$che_replay/spectral.json" --output-dir "$che_replay"
```

The release driver calculates all 68 tracked finite-matrix roots, including degrees
400, 480 and 640, both precision repeats and nine residual diagnostics.
With no grid or precision overrides, `qnm_spectral.py --output PATH` invokes the same
full generator. Custom-resolution runs retain their own partial-grid output.
The release driver selects terminal roots using grid differences computed from
the decimal root strings in multiprecision. Run generation commands in a working copy to preserve
historical timings and snapshots. Each command above uses the working directory
at the root of this archive.

For a single incoming amplitude without rerunning a campaign:

```python
import sys
sys.path.insert(0, 'analysis')
from qnm_contour import incoming_amplitude
A, diagnostics = incoming_amplitude(
    0.22-0.21j, ell=0, anchor=1.5, radius=100,
    order=48, rtol=3e-14, diagnostics=True)
print(A)
print(diagnostics['bounds'])
```

The simple tail bounds returned by this interface require
`max(0,-Im(Omega)) <= 2 Re(Omega) sqrt(anchor*(anchor-1))`.
Choose a fixed anchor that satisfies this on an entire frequency patch if using
analytic root-counting arguments. The report gives a more general existence
bound throughout the right half-plane.

## What the results establish

Nine low scalar Schwarzschild modes were found by contour matching and compared
with independently implemented continued-fraction and spectral calculations.
Three contour roots were repeated in multiprecision, including the sensitive
`ell=0, n=2` mode. The code does not import published QNM tables as input.
The continued-fraction reference is a finite-depth numerical result, not assumed
exact truth. Numerical comparisons and refinement histories are kept separately.

The report proves analytic truncation estimates for the outgoing endpoint,
finite Volterra expansion, and horizon series. **The implementations do not
provide interval-certified propagation, rounding, or QNM frequency enclosures.**
A tiny tail estimate or amplitude residual is not a bound on the full error.
No spectral completeness or performance superiority is claimed.

The preliminary nodal spectral experiment is retained under `qnm_spectral_nodal*`
for transparency; its poorly converged entries are not used in the final table.
Exact coefficient assembly and higher precision were both changed in the final
spectral implementation. Observed stability is not a condition-number comparison.

## Prior art and research positioning

This is an independent derivation and reproducible computational specialization.
It is not claimed to be the first connection-coefficient QNM method, first
rotated-Volterra Schwarzschild construction, or first rigorous endpoint-error
architecture. The report discusses Correia et al., Chen et al., Leaver/MST,
Donninger–Schlag–Soffer, Casals–Ottewill, Motygin, and related work explicitly.
Specific novelty or an efficiency advantage remains to be demonstrated.
