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
        self.tree = tree
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
    # In transfer mode F5 gains the across-latitude panel of SPEC_04 Step 5 deliverable 4, so the
    # figure is two panels wide there and one in closure mode.
    transfer = _is_transfer(pr)
    fig = Figure(figsize=((style.FIGSIZE_WIDE[0] * 1.6, style.FIGSIZE_WIDE[1]) if transfer
                          else style.FIGSIZE_WIDE))
    ax = fig.add_subplot(1, 2, 1) if transfer else fig.subplots(1, 1)
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
    data = {"pressure_Pa": pr.p, "geopotential_m2s2": pr.phi,
            "height_above_anchor_isobar_m": pr.h}
    if transfer:
        data.update(figure_5_panel(pr, fig))
    return fig, data


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


def _is_transfer(pr: _Profile) -> bool:
    """Whether the product is a transfer run, by the mode its own attributes declare."""
    return str(pr.root.attrs.get("mode", "")) == "transfer"


def figure_5_panel(pr: _Profile, fig):
    """F5's across-latitude panel: the reference surface and the traced isobars in (phi, Phi).

    SPEC_04 Step 5 deliverable 4. Every anchor's isobars, the reference surface on the same
    latitude nodes, and the three latitudes that matter marked: each anchor, the gauge and the
    target. The isobars are drawn thinned, because 66 curves over 400 nodes is a black band and the
    shape is what the panel is for.
    """
    isobars = pr.tree["isobars"].to_dataset(inherit=False)
    surface = pr.tree["reference_surface"].to_dataset(inherit=False)
    ax = fig.add_subplot(1, 2, 2)
    latitude = np.asarray(isobars["latitude_planetocentric_deg"].values, dtype="float64")
    drawn = {}
    for name in sorted(str(n) for n in isobars.data_vars):
        if not name.startswith("geopotential_") or not name.endswith("_m2s2"):
            continue
        if name == "geopotential_gauge_to_target_m2s2":
            continue
        slug = name[len("geopotential_"):-len("_m2s2")]
        curves = np.asarray(isobars[name].values, dtype="float64")
        step = max(1, curves.shape[0] // 12)
        for row in range(0, curves.shape[0], step):
            ax.plot(latitude, curves[row] / 1e6, color=style.COLOR["geopotential"],
                    linewidth=0.6, alpha=0.55,
                    label=f"isobars of {slug}" if row == 0 else None)
        drawn[slug] = curves
    # The reference surface is the gauge isobar, so in the (phi, Phi) plane it is the line
    # Phi = 0 at every latitude; what varies with latitude is its radius, which the
    # `reference_surface` group carries and the right axis states at the two ends.
    ax.plot(latitude, np.zeros_like(latitude), color=style.ANCHOR_LINE["color"], linewidth=1.2,
            label="the reference surface, Phi = 0")
    radius = np.asarray(surface["reference_surface_radius_m"].values, dtype="float64")
    ax.set_xlabel(f"planetocentric latitude (degrees); r0 runs {radius[0] / 1e3:,.0f} to "
                  f"{radius[-1] / 1e3:,.0f} km across these nodes")

    # The marks go in the axes' own coordinates so that they sit inside the panel whatever the
    # data range is, and the target and the gauge are labeled apart because at M = 1 they and the
    # one anchor are the same latitude.
    marks = [(float(pr.root["latitude_planetocentric_deg"].values), "target",
              style.COLOR["height"], "--", 0.04)]
    for slug in drawn:
        marks.append((float(pr.tree[f"anchors/{slug}/transfer"]
                            .attrs["latitude_planetocentric_deg"]), slug, "0.35", ":", 0.16))
    marks.append((float(pr.root["gauge_latitude_planetocentric_deg"].values), "gauge",
                  style.ANCHOR_LINE["color"], "-.", 0.10))
    for where, text, color, dash, height in marks:
        ax.axvline(where, color=color, linestyle=dash, linewidth=1.0)
        ax.annotate(f"{text} {where:g}", (where, height), xycoords=("data", "axes fraction"),
                    textcoords="offset points", xytext=(3, 0), fontsize=7, color=color)
    ax.set_ylabel("geopotential (millions of m2/s2)")
    ax.set_title("the traced isobars across latitude")
    ax.legend(loc="best", fontsize=7)
    return {"latitude_planetocentric_deg": latitude,
            **{f"isobars_{k}": v for k, v in drawn.items()}}


def figure_6_transfer(pr: _Profile):
    """F6 in transfer mode: the pressure identity residual against the label pressure.

    Decision H's internal check of the whole chain, drawn where closure mode draws the hydrostatic
    closure against the tabulated pressure. There is no print-rounding envelope here: the label is
    a number the chain carries, not a printed measurement, so the residual is compared against
    nothing but zero and its own largest magnitude.
    """
    fig = Figure(figsize=style.FIGSIZE_WIDE)
    ax = fig.subplots(1, 1)
    style.layout(fig, rows=1)
    labels = np.asarray(pr.root["pressure_label_Pa"].values, dtype="float64")
    residual = np.asarray(pr.root["pressure_identity_residual"].values, dtype="float64")
    worst = int(np.argmax(np.abs(residual)))
    ax.axvline(0.0, color="0.6", linewidth=0.8)
    ax.plot(residual, labels / 100.0, marker="o", markersize=2.5,
            color=style.COLOR["hydrostatic"], label="p_produced / p_label - 1")
    ax.plot([residual[worst]], [labels[worst] / 100.0], marker="o", markersize=6,
            markerfacecolor="none", color=style.COLOR["hydrostatic"],
            label=f"largest {residual[worst]:+.2e} at level {worst}")
    style.pressure_axis(ax, labels)
    _mark_levels(ax, pr, labels)
    ax.set_xlabel("fractional difference")
    ax.set_title("the pressure identity against the isobar label (decision H)")
    ax.legend(loc="lower left")
    fig.suptitle(f"F6. Pressure identity, {pr.slug} (anchor {pr.anchor_slug})", fontsize=10)
    return fig, {"pressure_label_Pa": labels, "pressure_identity_residual": residual,
                 "largest_abs": float(np.abs(residual).max()), "largest_level": worst}


def figure_7(pr: _Profile):
    """F7, transfer mode only: what was delivered, beside each anchor's own.

    Three panels: the delivered `T(p)` and each anchor's on the same axes; `N` against the
    geopotential, delivered and at the gauge; and the isobar shift against the label pressure. At
    M >= 2 a fourth panel draws every `D_ij`, which is empty at M = 1 and is left out there rather
    than drawn blank.
    """
    estimate = pr.tree["estimate"].to_dataset(inherit=False)
    differences = [str(n) for n in estimate.data_vars if str(n).startswith("D_")]
    columns = 4 if differences else 3
    fig = Figure(figsize=(style.FIGSIZE_WIDE[0], style.FIGSIZE_WIDE[1]))
    axes = fig.subplots(1, columns)
    style.layout(fig, rows=1)
    labels = np.asarray(pr.root["pressure_label_Pa"].values, dtype="float64")
    T = np.asarray(pr.root["temperature_K"].values, dtype="float64")
    N = np.asarray(pr.root["refractivity"].values, dtype="float64")
    N_gauge = np.asarray(pr.root["refractivity_gauge"].values, dtype="float64")
    phi = np.asarray(pr.root["geopotential_m2s2"].values, dtype="float64")

    axes[0].plot(T, labels / 100.0, color=style.COLOR["hydrostatic"], label="delivered")
    for slug in sorted(pr.tree["anchors"].children):
        thermo = pr.tree[f"anchors/{slug}/inputs/thermo"].to_dataset(inherit=False)
        axes[0].plot(np.asarray(thermo["temperature_K"].values, dtype="float64"),
                     np.asarray(thermo["pressure_Pa"].values, dtype="float64") / 100.0,
                     linewidth=0.9, linestyle="--", label=f"{slug}, tabulated")
    style.pressure_axis(axes[0], labels)
    axes[0].set_xlabel("temperature (K)")
    axes[0].set_title("delivered T against each anchor's")
    axes[0].legend(loc="best", fontsize=7)

    axes[1].plot(N, phi / 1e6, color=style.COLOR["geopotential"], label="delivered at the target")
    axes[1].plot(N_gauge, np.asarray(pr.root["geopotential_gauge_m2s2"].values,
                                     dtype="float64") / 1e6,
                 linestyle="--", color=style.COLOR["height"], label="exp C at the gauge")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("refractivity")
    axes[1].set_ylabel("geopotential (millions of m2/s2)")
    axes[1].set_title("refractivity against geopotential")
    axes[1].legend(loc="best", fontsize=7)

    shift = phi - np.asarray(pr.root["geopotential_gauge_m2s2"].values, dtype="float64")
    axes[2].axvline(0.0, color="0.6", linewidth=0.8)
    axes[2].plot(shift, labels / 100.0, marker="o", markersize=2.0,
                 color=style.COLOR["geopotential"])
    style.pressure_axis(axes[2], labels)
    axes[2].set_xlabel("geopotential (m2/s2)")
    axes[2].set_title("the isobar shift, gauge to target")

    data = {"pressure_label_Pa": labels, "temperature_K": T, "refractivity": N,
            "refractivity_gauge": N_gauge, "isobar_shift_m2s2": shift}
    if differences:
        union = np.asarray(estimate["label_pressure_Pa"].values, dtype="float64")
        axes[3].axvline(0.0, color="0.6", linewidth=0.8)
        for name in sorted(differences):
            value = np.asarray(estimate[name].values, dtype="float64")
            axes[3].plot(value, union / 100.0, marker="o", markersize=2.0,
                         label=name.replace("_", " "))
            data[name] = value
        floor = float(pr.tree["estimate"].attrs.get("identity_floor", np.nan))
        if np.isfinite(floor):
            axes[3].axvspan(-floor, floor, color=style.COLOR["band"], alpha=0.45, linewidth=0,
                            label=f"identity floor {floor:.1e}")
        style.pressure_axis(axes[3], union)
        axes[3].set_xlabel("difference in ln N")
        axes[3].set_title("D_ij between the anchors (Eq. A33)")
        axes[3].legend(loc="best", fontsize=7)

    fig.suptitle(f"F7. Delivered profile, {pr.slug} at "
                 f"{float(pr.root['latitude_planetocentric_deg'].values):g} degrees", fontsize=10)
    return fig, data


FIGURES = (
    ("F5", "geopotential", figure_5, None),
    ("F6", "hydrostatic", figure_6, "pressure_tabulated_Pa"),
)

#: SPEC_04 Step 5 deliverable 4: what a transfer product draws. F5 keeps its name and gains a
#: panel, F6 draws the pressure identity where the closure draws the hydrostatic closure, and F7 is
#: the delivered profile beside each anchor's.
FIGURES_TRANSFER = (
    ("F5", "geopotential", figure_5, None),
    ("F6", "pressure_identity", figure_6_transfer, "pressure_label_Pa"),
    ("F7", "delivered", figure_7, "pressure_label_Pa"),
)


def build_all(tree):
    """F5 and F6 in the mode the product declares, and F7 in transfer mode.

    Returns `(figures, skipped, data)`.
    """
    pr = _Profile(tree)
    figures, skipped, data = [], {}, {}
    for key, name, maker, needs in (FIGURES_TRANSFER if _is_transfer(pr) else FIGURES):
        if needs is not None and needs not in pr.root.variables:
            skipped[f"{key}_{name}"] = f"needs {needs}(level), which a profile carries in closure mode"
            continue
        fig, values = maker(pr)
        figures.append((key, name, fig))
        data[key] = values
    return figures, skipped, data
