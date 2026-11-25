#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""The complex frequency plane of the article.

The figure shows the nine isolated scalar quasinormal frequencies, read from the stored
records cert/results/<mode>_center.json, the lines that bound the regions where the outgoing
solution along the upward ray z = a + is is damped (p > 0 and q rho_a < p, with
Omega = p - iq and rho_a = 1/(2 sqrt(a(a-1)))), and the special points of the negative
imaginary axis.

Run from the root of the repository:

    python3 docs/figures/frequency_plane.py
    pdftocairo -svg docs/figures/frequency-plane.pdf docs/figures/frequency-plane.svg

The first command writes the PDF, the second the SVG shown in the README.
"""
import json
from fractions import Fraction
from math import sqrt
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

MODES = ("00", "01", "02", "10", "11", "12", "20", "21", "22")
ANCHORS = (1.5, 2.0, 3.0, 5.0)
RE_RANGE = (-0.74, 1.10)
IM_RANGE = (-1.36, 0.10)


def centre(mode):
    """The rational centre and the anchor of the ray used for one mode."""
    record = json.loads((ROOT / "cert" / "results" / f"{mode}_center.json").read_text())
    re, im = (float(Fraction(value)) for value in record["omega_exact"])
    return re, im, float(Fraction(str(record["anchor"])))


def slope(a):
    """Magnitude of the slope of the line q rho_a = p, that is 2 sqrt(a(a-1))."""
    return 2.0 * sqrt(a * (a - 1.0))


def main():
    hs.configure(family="latin-modern")
    fig, ax = plt.subplots(figsize=(hs.WIDTH, 5.35), layout="constrained")

    # Guides: the real axis and the upper half of the imaginary axis.
    ax.plot([-0.05, RE_RANGE[1]], [0.0, 0.0], color=hs.MIST, lw=0.6, zorder=0)
    ax.plot([0.0, 0.0], [0.0, IM_RANGE[1]], color=hs.MIST, lw=0.6, zorder=0)

    # The branch cut on the negative imaginary axis.
    ax.plot([0.0, 0.0], [0.0, IM_RANGE[0]], color=hs.OCHRE, lw=2.0,
            solid_capstyle="butt", zorder=1)
    ax.text(-0.025, -1.02, "branch cut", ha="right", va="center", color=hs.OCHRE)

    # Lines bounding the damped regions of the four anchors used in the proofs.
    for k, a in enumerate(ANCHORS):
        s = slope(a)
        y_end = IM_RANGE[0]
        ax.plot([0.0, -y_end / s], [0.0, y_end], color=hs.SLATE, lw=0.8,
                linestyle=hs.DASHES[(k + 1) % len(hs.DASHES)], zorder=1)
        label = f"$a={a:g}$"
        ax.text(1.29 / s - 0.014, -1.29, label, ha="right", va="center", color=hs.SLATE)

    # The nine frequencies, one marker and one colour for each angular momentum.
    points = [centre(mode) for mode in MODES]
    for ell in range(3):
        xs = [p[0] for mode, p in zip(MODES, points) if int(mode[0]) == ell]
        ys = [p[1] for mode, p in zip(MODES, points) if int(mode[0]) == ell]
        ax.plot(xs, ys, linestyle="none", marker=hs.MARKERS[ell], markersize=6.0,
                markerfacecolor="white", markeredgecolor=hs.SERIES[ell],
                markeredgewidth=1.2, label=rf"$\ell={ell}$", zorder=3)

    # Special points of the negative imaginary axis and the threshold half-disk.
    for y in (-0.25, -0.5):
        ax.plot([0.0], [y], linestyle="none", marker="o", markersize=4.5,
                markerfacecolor=hs.PAPER, markeredgecolor=hs.INK, markeredgewidth=1.0,
                zorder=4)
    notes = (
        ((0.0, 0.0), (-0.72, 0.035),
         "$0<|\\Omega|\\leqslant 2^{-10}$:\nno monopole resonance"),
        ((0.0, -0.25), (-0.72, -0.27), "$\\Omega=-i/4$:\nthe two banks differ"),
        ((0.0, -0.5), (-0.72, -0.60), "$\\Omega=-i/2$:\nthe horizon pole\ncancels ($\\ell=0$)"),
    )
    for target, where, text in notes:
        ax.annotate(text, xy=target, xytext=where, ha="left", va="center", color=hs.INK,
                    arrowprops=dict(arrowstyle="-", color=hs.SLATE, lw=0.6,
                                    shrinkA=2.0, shrinkB=3.0))

    ax.set_xlim(*RE_RANGE)
    ax.set_ylim(*IM_RANGE)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\mathrm{Re}\,\Omega$")
    ax.set_ylabel(r"$\mathrm{Im}\,\Omega$")
    hs.legend_outside(ax, where="above", ncol=3)
    hs.save(fig, HERE / "frequency-plane.pdf")


if __name__ == "__main__":
    main()
