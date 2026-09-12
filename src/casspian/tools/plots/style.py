"""The one style of every diagnostic figure. SPEC_02 Step 5.

Fonts, line widths and the color of each quantity are defined here and nowhere else. matplotlib
only: no seaborn, no external stylesheet. The style is applied through a context manager so a
caller's own matplotlib settings are left as they were.

**The footer.** Every figure carries the product file name, its `casspian_git_commit`, the first
twelve characters of its SHA-256, and the generation time. The generation time is the only line
that changes from one rendering of the same file to the next, so it is drawn alone in a band of
`FOOTER_TIME_BAND` of the figure height at the very bottom, and nothing else is drawn there. A
comparison of two renderings masks that band and compares everything else.
"""

from __future__ import annotations

from contextlib import contextmanager

import matplotlib

RC = {
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8.5,
    "legend.fontsize": 7,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "lines.linewidth": 1.4,
    "lines.markersize": 4,
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
}

#: The color of each quantity, the same in every figure that shows it.
COLOR = {
    "temperature": "#c0392b",
    "wind": "#1f77b4",
    "H2": "#2c7fb8",
    "He": "#7b3294",
    "NH3": "#1a9641",
    "molar_mass": "#8c510a",
    "mean_refractivity": "#01665e",
    "number_density": "#4d4d4d",
    "refractivity": "#01665e",
    "gravity_newton": "#555555",
    "gravity_effective": "#d95f02",
    "centrifugal": "#7570b3",
    "wind_term": "#1f77b4",
    "G_phi": "#e7298a",
    "psi": "#66a61e",
    "geoid_nowind": "#999999",
    "geoid_wind": "#1f77b4",
    "dynamical_height": "#d95f02",
    "recovered": "#c0392b",
    "geopotential": "#5e3c99",
    "height": "#b2abd2",
    "hydrostatic": "#e66101",
    "anchor": "#333333",
    "latitude": "#333333",
    "band": "#9ecae1",
    "parameterized": "#ffe9b8",
    "extrapolated": "#f6d6d6",
    "extended": "#e0ecf4",
}

#: A species not listed above takes its color from this cycle, in the order the file lists it.
SPECIES_CYCLE = ("#e6ab02", "#a6761d", "#666666", "#1b9e77")

#: Markers for the per level provenance of a trace species, by flag meaning.
PROVENANCE_MARKER = {"measured": "o", "interpolated": "s", "assumed": "x", "extrapolated": "^"}

FIGSIZE = (11.0, 8.5)
FIGSIZE_WIDE = (11.0, 5.0)
FOOTER_TIME_BAND = 0.025
FOOTER_TOP = 0.065
FOOTER_FONTSIZE = 6.5

ANCHOR_LINE = {"color": COLOR["anchor"], "linestyle": "--", "linewidth": 0.9}
LATITUDE_LINE = {"color": COLOR["latitude"], "linestyle": ":", "linewidth": 1.0}


@contextmanager
def styled():
    """Apply the diagnostics style for the duration of a figure's creation and saving."""
    with matplotlib.rc_context(RC):
        yield


def species_color(name: str, index: int) -> str:
    return COLOR.get(name, SPECIES_CYCLE[index % len(SPECIES_CYCLE)])


def pressure_axis(ax, pressure_Pa) -> None:
    """Pressure on the vertical axis: log scale, decreasing upward, labelled in mbar.

    The file carries Pa; the axis shows mbar, which is how every source prints pressure. The
    conversion is for display only.
    """
    import numpy as np

    p = np.asarray(pressure_Pa, dtype="float64")
    ax.set_yscale("log")
    ax.set_ylim(float(np.nanmax(p)) / 100.0 * 1.15, float(np.nanmin(p)) / 100.0 / 1.15)
    ax.set_ylabel("pressure (mbar)")


def anchor_line(ax, anchor_Pa, label: bool = True) -> None:
    """The anchor isobar, a horizontal line on every profile panel."""
    ax.axhline(float(anchor_Pa) / 100.0, **ANCHOR_LINE,
               label=f"anchor isobar {float(anchor_Pa) / 100.0:g} mbar" if label else None)


def latitude_line(ax, latitude_deg, label: bool = True) -> None:
    """The profile latitude, a vertical line on every latitude panel."""
    ax.axvline(float(latitude_deg), **LATITUDE_LINE,
               label=f"phi_c {float(latitude_deg):.4f} deg" if label else None)


def layout(fig, rows: int = 2) -> None:
    """Fixed margins, so the footer bands are the same on every figure."""
    # A one row figure is short, so its margins take a larger fraction of the height: room for
    # a panel title raised by a twin axis below the suptitle, and for the x label above the footer.
    fig.subplots_adjust(left=0.07, right=0.97, top=0.86 if rows > 1 else 0.78,
                        bottom=FOOTER_TOP + (0.06 if rows > 1 else 0.12), hspace=0.50,
                        wspace=0.30)


def footer(fig, file_name: str, commit: str, sha12: str, generated_at: str) -> str:
    """Draw the footer and return its traceable part, the text above the time band."""
    identity = f"{file_name}    casspian_git_commit {commit}    sha256 {sha12}"
    fig.text(0.01, FOOTER_TIME_BAND + 0.008, identity, fontsize=FOOTER_FONTSIZE, color="0.3",
             ha="left", va="bottom")
    fig.text(0.01, 0.004, f"generated {generated_at}", fontsize=FOOTER_FONTSIZE, color="0.3",
             ha="left", va="bottom")
    return identity
