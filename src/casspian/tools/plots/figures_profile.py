"""The kind `profile` figures, F5 and F6. SPEC_03 Step 3 deliverable 3.

F5 is the geopotential against the hydrostatic pressure, with the height above the anchor isobar
on a twin axis; it renders for every profile. F6 is the hydrostatic closure,
`pressure_Pa / pressure_tabulated_Pa - 1` against pressure, with the boundary and gauge levels
marked and the print-rounding envelope of the anchor drawn as a band; it renders only when the
tabulated pair is present (closure mode). The across-latitude panel of SPEC_02's F5 is dropped
until SPEC_04 builds the surface.

**The envelope, from the anchor's own file.** Two terms per level, added:

* half a unit in the last printed figure of the tabulated pressure, divided by the pressure, where
  the level is not on the grid its kind T declares in `pressure_grid_rule`, and zero where it is
  (after SPEC_03 Step 0 every level but the last is on the grid, and a grid value carries no
  print rounding);
* half a unit in the last printed figure of the tabulated height, divided by the local scale
  height `R T / (M g)`, with `T` the tabulated temperature, `M` the mean molar mass, and `g` the
  magnitude `|dPhi/dh|` of the profile's own geopotential against its height.

The printed precision is read from the numbers as the kind T file carries them in the units the
source prints (mbar, km): for each column, the finest decimal place any of its values shows in its
shortest representation. No precision is written into this module.
"""

from __future__ import annotations

import re

import numpy as np
from matplotlib.figure import Figure

from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.tools.plots import style


def _decimals(value: float) -> int:
    text = repr(float(value))
    if "e" in text or "E" in text:
        return 0
    return len(text.split(".")[1].rstrip("0")) if "." in text else 0


def _half_unit(values) -> float:
    """Half a unit in the last printed figure of a tabulated column.

    A column is printed to one precision, and a value such as 90.0 km shows none of its printed
    zeros once it is a float, so the precision is the finest one any value of the column shows.
    """
    values = np.asarray(values, dtype="float64")
    finest = max(_decimals(v) for v in values[np.isfinite(values)])
    return 0.5 * 10.0 ** (-finest)


class _Profile:
    """The profile, its first anchor, and the anchor's embedded kind T, read once."""

    def __init__(self, tree):
        self.root = tree.to_dataset(inherit=False)
        self.slug = str(self.root.attrs.get("profile_or_run", ""))
        anchors = tree["anchors"]
        self.anchor_slug = sorted(anchors.children)[0]
        anchor = anchors[self.anchor_slug]
        self.thermo = anchor["inputs/thermo"].to_dataset(inherit=False)
        self.p = np.asarray(self.root["pressure_Pa"].values, dtype="float64")
        self.phi = np.asarray(self.root["geopotential_m2s2"].values, dtype="float64")
        self.h = np.asarray(self.root["height_above_anchor_isobar_m"].values, dtype="float64")
        self.gauge = int(self.root["gauge_level_index"].values)
        self.boundary = int(self.root["boundary_level_index"].values)

    def envelope(self):
        """The two terms of the print-rounding envelope, per level, as fractions of pressure."""
        p_tab = np.asarray(self.root["pressure_tabulated_Pa"].values, dtype="float64")
        T = np.asarray(self.root["temperature_tabulated_K"].values, dtype="float64")
        M = np.asarray(self.root["mean_molar_mass_kg_mol"].values, dtype="float64")
        printed_mbar = np.asarray(self.thermo["pressure_printed_Pa"].values, dtype="float64") / 100.0
        height_km = np.asarray(self.thermo["height_m"].values, dtype="float64") / 1000.0

        rule = str(self.thermo.attrs.get("pressure_grid_rule", ""))
        match = re.search(r"10\^\(k/(\d+)\) mbar", rule)
        on_grid = np.zeros(p_tab.shape, dtype=bool)
        if match:
            denominator = int(match.group(1))
            k = np.round(denominator * np.log10(p_tab / 100.0))
            on_grid = np.abs(100.0 * 10.0 ** (k / denominator) / p_tab - 1.0) < 1e-12
        pressure_term = np.where(on_grid, 0.0, _half_unit(printed_mbar) * 100.0 / p_tab)

        g = np.abs(np.gradient(self.phi, self.h))
        scale_height = MOLAR_GAS_CONSTANT * T / (M * g)
        height_term = _half_unit(height_km) * 1000.0 / scale_height
        return pressure_term, height_term, on_grid, scale_height


def _mark_levels(ax, pr: _Profile, pressure):
    ax.axhline(pressure[pr.boundary] / 100.0, color="0.3", linestyle="-.", linewidth=0.9,
               label=f"boundary level {pr.boundary}, {pressure[pr.boundary] / 100.0:g} mbar")
    # The gauge is defined by the tabulated level, so it is labeled with the declared
    # gauge_isobar_Pa, not with the hydrostatic pressure found there (REVIEW_03_step3 ruling 6).
    gauge_Pa = float(pr.root["gauge_isobar_Pa"].values)
    ax.axhline(pressure[pr.gauge] / 100.0, **style.ANCHOR_LINE,
               label=f"gauge level {pr.gauge}, {gauge_Pa / 100.0:g} mbar (Phi = 0)")


def figure_5(pr: _Profile):
    fig = Figure(figsize=style.FIGSIZE_WIDE)
    ax = fig.subplots(1, 1)
    style.layout(fig, rows=1)
    ax.plot(pr.phi, pr.p / 100.0, color=style.COLOR["geopotential"], label="geopotential")
    ax.set_xlabel("geopotential (m2/s2)")
    twin = ax.twiny()
    twin.plot(pr.h / 1e3, pr.p / 100.0, color=style.COLOR["height"], linestyle="--",
              label="height above the anchor isobar")
    twin.set_xlabel("height above the anchor isobar, along the local vertical (km)",
                    color=style.COLOR["height"])
    twin.grid(False)
    style.pressure_axis(ax, pr.p)
    _mark_levels(ax, pr, pr.p)
    handles, labels = ax.get_legend_handles_labels()
    extra = twin.get_legend_handles_labels()
    ax.legend(handles + extra[0], labels + extra[1], loc="lower left")
    ax.set_title("geopotential against the hydrostatic pressure, with height", pad=22)
    fig.suptitle(f"F5. Geopotential, {pr.slug} (anchor {pr.anchor_slug})", fontsize=10)
    return fig, {"pressure_Pa": pr.p, "geopotential_m2s2": pr.phi,
                 "height_above_anchor_isobar_m": pr.h}


def figure_6(pr: _Profile):
    fig = Figure(figsize=style.FIGSIZE_WIDE)
    ax = fig.subplots(1, 1)
    style.layout(fig, rows=1)
    p_tab = np.asarray(pr.root["pressure_tabulated_Pa"].values, dtype="float64")
    residual = pr.p / p_tab - 1.0
    pressure_term, height_term, on_grid, scale_height = pr.envelope()
    band = pressure_term + height_term
    ax.fill_betweenx(p_tab / 100.0, -band, band, color=style.COLOR["band"], alpha=0.45,
                     linewidth=0, label="print-rounding envelope of the anchor")
    ax.axvline(0.0, color="0.6", linewidth=0.8)
    ax.plot(residual, p_tab / 100.0, marker="o", markersize=2.5,
            color=style.COLOR["hydrostatic"], label="p_hydro / p_tab - 1")
    style.pressure_axis(ax, p_tab)
    _mark_levels(ax, pr, p_tab)
    ax.set_xlabel("fractional difference")
    ax.set_title("hydrostatic closure against the tabulated pressure")
    ax.legend(loc="lower left")
    fig.suptitle(f"F6. Hydrostatic closure, {pr.slug} (anchor {pr.anchor_slug})", fontsize=10)
    inside = np.abs(residual) <= band
    return fig, {"pressure_tabulated_Pa": p_tab, "residual": residual, "envelope": band,
                 "envelope_pressure_term": pressure_term, "envelope_height_term": height_term,
                 "on_grid": on_grid, "scale_height_m": scale_height,
                 "inside_envelope": inside}


FIGURES = (
    ("F5", "geopotential", figure_5, None),
    ("F6", "hydrostatic", figure_6, "pressure_tabulated_Pa"),
)


def build_all(tree):
    """F5, and F6 when the tabulated pair is present. Returns `(figures, skipped, data)`."""
    pr = _Profile(tree)
    figures, skipped, data = [], {}, {}
    for key, name, maker, needs in FIGURES:
        if needs is not None and needs not in pr.root.variables:
            skipped[f"{key}_{name}"] = f"needs {needs}(level), which a profile carries in closure mode"
            continue
        fig, values = maker(pr)
        figures.append((key, name, fig))
        data[key] = values
    return figures, skipped, data
