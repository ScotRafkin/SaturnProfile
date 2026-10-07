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
from casspian.lib.windfield import WindField
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
def _five_isobars(pr: _Profile, labels):
    """The five levels F8 draws along: the top, 10 mbar, the gauge, 1 bar and the bottom.

    Named by their label pressure, and duplicates dropped, so a profile whose gauge is 10 mbar
    draws four curves rather than the same one twice.
    """
    wanted = [("the top level", 0),
              ("10 mbar", int(np.argmin(np.abs(labels - 1.0e3)))),
              ("the gauge", pr.gauge),
              ("1 bar", int(np.argmin(np.abs(labels - 1.0e5)))),
              ("the bottom", labels.size - 1)]
    seen, out = set(), []
    for name, level in wanted:
        if level not in seen:
            seen.add(level)
            out.append((f"{name}, {labels[level] / 100.0:.3g} mbar", level))
    return out


def figure_8(pr: _Profile):
    """F8, transfer mode only: the transfer along the isobars, and T on altitude and on radius.

    SPEC_04 v0.18 deliverable 4, the author's request on viewing F7. Four panels: `S/g` against
    latitude on five isobars over the span the transfer crossed, from the `isobars` group's
    `shear_kernel_<slug>_per_rad`; the temperature change accumulated along the same isobars; and
    the delivered `T` against altitude and against radius, each beside the first anchor's own `T`
    from the kind N file. On the altitude panel both profiles are above their own datum level
    (v0.21), the anchor's located on its kind N heights by the log-linear rule the datum itself is
    placed with, so the panel compares the thickness of the two columns; on the radius panel the
    offset between them is the reference surface between the two latitudes, which is the point of
    that panel.

    **The staircase in the first panel is the wind file, not the mesh** (SPEC_04 section 18
    ruling 5). Under decision L the wind is linear between the file's half-degree nodes, so
    `S = 2 Omega cos(phi) du/dphi` is constant on each interval and steps at every node; a mesh at
    0.05 degrees resolves those steps exactly and a finer one draws the same staircase. The
    alternation over a few degrees of the flank is the source curve's own node-to-node roughness,
    and the kinks in the second panel are its integral, which is why the accumulated offset is
    smooth: the alternation cancels. Smoothing belongs to the wind tool, not to the model, which
    reads what it is given (SPEC_00 section 3.4).

    Nothing is recomputed. The accumulated change is read off `ln N_k(phi)`, which is the integral
    of the kernel the first panel draws: on an isobar the pressure is the label, so
    `T = p_label R_bar / (k_B N)` and the change from the anchor's latitude is
    `T(phi) - T(phi_a)` with `R_bar` the run's own composition on that label, from the product's
    `inputs/composition` group. The curve therefore starts at zero at the anchor and ends at the
    delivered temperature less the model's temperature at the anchor on the same isobar, which is
    the offset F7 draws against the anchor's tabulated column.
    """
    from casspian.forward.production import mean_properties_on_labels

    isobars = pr.tree["isobars"].to_dataset(inherit=False)
    labels = np.asarray(pr.root["pressure_label_Pa"].values, dtype="float64")
    latitude = np.asarray(isobars["latitude_planetocentric_deg"].values, dtype="float64")
    slug = pr.anchor_slug
    kernel = np.asarray(isobars[f"shear_kernel_{slug}_per_rad"].values, dtype="float64")
    ln_N = np.asarray(isobars[f"ln_refractivity_{slug}"].values, dtype="float64")
    anchor_deg = float(pr.tree[f"anchors/{slug}/transfer"]
                       .attrs["latitude_planetocentric_deg"])
    target_deg = float(pr.root["latitude_planetocentric_deg"].values)
    gauge_deg = float(pr.root["gauge_latitude_planetocentric_deg"].values)
    span = (latitude >= min(gauge_deg, target_deg)) & (latitude <= max(gauge_deg, target_deg))
    at_anchor = int(np.argmin(np.abs(latitude - anchor_deg)))

    # R_bar on each isobar's label, at every latitude of the span, from the run's composition.
    composition = pr.tree["inputs/composition"]
    R_bar = np.stack([mean_properties_on_labels(composition, float(deg), labels)[0]
                      for deg in latitude[span]], axis=1)
    R_bar_anchor = mean_properties_on_labels(composition, anchor_deg, labels)[0]
    T_anchor_column = np.asarray(pr.root["temperature_K"].values, dtype="float64")

    fig = Figure(figsize=(style.FIGSIZE_WIDE[0] * 1.6, style.FIGSIZE_WIDE[1]))
    axes = fig.subplots(1, 4)
    style.layout(fig, rows=1)
    chosen = _five_isobars(pr, labels)
    data = {"latitude_planetocentric_deg": latitude[span]}
    # Under a wind that does not vary along a column the kernel is the same on every isobar and
    # the curves lie on one another, which is the file's own statement and worth seeing; the widths
    # step down so that a coincident set still shows every member.
    for position, (name, level) in enumerate(chosen):
        color = f"C{position}"
        width = 2.6 - 0.45 * position
        axes[0].plot(latitude[span], kernel[level][span], color=color, linewidth=width,
                     alpha=0.85, label=name)
        # T on the isobar, up to the constant p_label / k_B, which cancels in the difference.
        ratio = (np.exp(ln_N[level][at_anchor] - ln_N[level][span])
                 * R_bar[level] / R_bar_anchor[level])
        delta_T = T_anchor_column[level] * (ratio - 1.0)
        axes[1].plot(latitude[span], delta_T, color=color, linewidth=width, alpha=0.85,
                     label=name)
        data[f"shear_kernel_level_{level}"] = kernel[level][span]
        data[f"delta_temperature_level_{level}"] = delta_T
    for ax, title, ylabel in ((axes[0], "the shear term S/g along the isobars (Eq. A27)",
                               "S/g (per radian)"),
                              (axes[1], "the temperature change accumulated along them",
                               "T(phi) - T(anchor) (K)")):
        ax.axvline(anchor_deg, **style.ANCHOR_LINE)
        ax.axvline(target_deg, color=style.COLOR["height"], linestyle="--", linewidth=1.0)
        ax.axhline(0.0, color="0.6", linewidth=0.8)
        ax.set_xlabel("planetocentric latitude (degrees)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc="best", fontsize=6.5)

    altitude = np.asarray(pr.root["altitude_m"].values, dtype="float64")
    radius = np.asarray(pr.root["radius_m"].values, dtype="float64")
    datum = float(pr.root["datum_isobar_Pa"].values)
    anchor = pr.tree[f"anchors/{slug}"].to_dataset(inherit=False)
    T_anchor = np.asarray(pr.thermo["temperature_K"].values, dtype="float64")
    h_anchor = np.asarray(anchor["height_above_anchor_isobar_m"].values, dtype="float64")
    r_anchor = np.asarray(anchor["radius_m"].values, dtype="float64")
    # Both profiles above their own datum level (v0.21). The anchor's is located on its own kind N
    # heights by the rule the datum itself is placed with, log-linear in pressure, so the panel
    # compares the thickness of the two columns and not the two files' reference levels.
    p_anchor = np.asarray(pr.thermo["pressure_Pa"].values, dtype="float64")
    rising = np.argsort(np.log(p_anchor))
    h_anchor_datum = float(np.interp(np.log(datum), np.log(p_anchor)[rising], h_anchor[rising]))
    axes[2].plot(T_anchor_column, altitude / 1000.0, color=style.COLOR["temperature"],
                 label=f"delivered at {target_deg:g} degrees")
    axes[2].plot(T_anchor, (h_anchor - h_anchor_datum) / 1000.0, linestyle="--", linewidth=0.9,
                 color="0.35", label=f"{slug}, tabulated")
    axes[2].axhline(0.0, color="0.6", linewidth=0.8)
    axes[2].set_ylabel(f"height above the {datum / 100:g} mbar level (km)")
    axes[2].set_title("the delivered T on altitude, both on their own datum")
    axes[3].plot(T_anchor_column, radius / 1000.0, color=style.COLOR["temperature"],
                 label=f"delivered at {target_deg:g} degrees")
    # The offset named is the one at the gauge isobar, where the delivered radius is r0 at the
    # target and the anchor's is its own anchor isobar radius, so the difference is the reference
    # surface between the two latitudes and nothing else. At the bottom level it is 1,602 km, which
    # is the same quantity plus the two columns' thicknesses and is not what the words say.
    surface_km = abs(radius[pr.gauge] - r_anchor[pr.gauge]) / 1000.0
    axes[3].plot(T_anchor, r_anchor / 1000.0, linestyle="--", linewidth=0.9, color="0.35",
                 label=f"{slug}, tabulated, {surface_km:.0f} km below at the gauge")
    axes[3].text(0.03, 0.03, "the gap at the gauge isobar is the reference surface\n"
                             "between the two latitudes", transform=axes[3].transAxes,
                 fontsize=6.5, color="0.25", va="bottom",
                 bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1.5))
    axes[3].set_ylabel("radius (km)")
    axes[3].set_title("and on radius, from the center")
    for ax in (axes[2], axes[3]):
        ax.set_xlabel("temperature (K)")
        ax.legend(loc="best", fontsize=6.5)
    data.update({"altitude_m": altitude, "radius_m": radius, "temperature_K": T_anchor_column,
                 "anchor_temperature_K": T_anchor, "anchor_height_m": h_anchor,
                 "anchor_radius_m": r_anchor,
                 "levels_drawn": np.asarray([level for _, level in chosen])})
    fig.suptitle(f"F8. The transfer along the isobars, {pr.slug} at {target_deg:g} degrees, "
                 f"anchor {slug} at {anchor_deg:.3f}", fontsize=10)
    return fig, data


def _bracketing(nodes, low, high):
    """The nodes from the last at or below `low` to the first at or above `high`, as a mask."""
    nodes = np.asarray(nodes, dtype="float64")
    below, above = nodes[nodes <= low], nodes[nodes >= high]
    start = below.max() if below.size else nodes.min()
    stop = above.min() if above.size else nodes.max()
    return (nodes >= start) & (nodes <= stop)


def figure_9(pr: _Profile):
    """F9, transfer mode only: the wind the run assumed, from the product's own inputs.

    SPEC_04 v0.19 deliverable 4, the author's request. `u_reference(phi)` over the file's whole
    latitude range, with the anchors, the gauge and the target marked and the span the transfer
    crossed shaded, and `u_total(phi, p)` as filled contours over the mesh's latitude range and the
    delivered profile's pressure range. The file's `source` and `vertical_structure` are in the
    second panel's title, so that a run under a sheared hypothesis shows the shear it was given.
    Nothing is recomputed: both panels are the `inputs/wind` group as the run read it.

    The author's direction of 28 September 2026 is that the posterior wind of the combination
    specification is drawn here as a second `u(phi, p)` panel once that specification produces it.
    Nothing here anticipates it.
    """
    wind = pr.tree["inputs/wind"].to_dataset(inherit=False)
    isobars = pr.tree["isobars"].to_dataset(inherit=False)
    latitude = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    pressure = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    # SPEC_08 section 6 ruling 10: the source line is drawn only where the file carries it.
    u_reference = (np.asarray(wind["u_reference_ms"].values, dtype="float64")
                   if "u_reference_ms" in wind.variables else None)
    u_total = np.asarray(wind["u_total_ms"].values, dtype="float64")
    reference_Pa = float(wind["reference_level_pressure_Pa"])
    # SPEC_05 Step 3 deliverable 5: u_total at the source's reference pressure, read by the model's
    # own rule, beside the source wind; where they agree the two lines lie on top of each other.
    u_total_at_reference = np.asarray(
        WindField(wind).wind_at(np.radians(latitude), np.full(latitude.shape, reference_Pa)),
        dtype="float64")
    mesh_latitude = np.asarray(isobars["latitude_planetocentric_deg"].values, dtype="float64")
    labels = np.asarray(pr.root["pressure_label_Pa"].values, dtype="float64")
    target_deg = float(pr.root["latitude_planetocentric_deg"].values)
    gauge_deg = float(pr.root["gauge_latitude_planetocentric_deg"].values)

    fig = Figure(figsize=style.FIGSIZE_WIDE)
    axes = fig.subplots(1, 2)
    style.layout(fig, rows=1)
    if u_reference is not None:
        axes[0].plot(latitude, u_reference, color=style.COLOR["wind"], linewidth=1.0,
                     label=f"source wind, assigned to {reference_Pa / 100:g} mbar")
    axes[0].plot(latitude, u_total_at_reference, color=style.COLOR["gravity_effective"],
                 linewidth=1.0, linestyle="--", label=f"u_total at {reference_Pa / 100:g} mbar")
    axes[0].axvspan(mesh_latitude.min(), mesh_latitude.max(), color=style.COLOR["band"],
                    alpha=0.35, linewidth=0,
                    label=f"the mesh, {mesh_latitude.min():g} to {mesh_latitude.max():g} degrees")
    for slug in sorted(pr.tree["anchors"].children):
        anchor_deg = float(pr.tree[f"anchors/{slug}/transfer"]
                           .attrs["latitude_planetocentric_deg"])
        axes[0].axvline(anchor_deg, color="0.35", linestyle=":", linewidth=0.9,
                        label=f"anchor {slug}, {anchor_deg:.3f}")
    axes[0].axvline(gauge_deg, **style.ANCHOR_LINE, label=f"gauge {gauge_deg:.3f}")
    axes[0].axvline(target_deg, color=style.COLOR["height"], linestyle="--", linewidth=1.0,
                    label=f"target {target_deg:g}")
    axes[0].set_xlabel("planetocentric latitude (degrees)")
    axes[0].set_ylabel("u (m/s)")
    axes[0].set_title("the wind at the reference level")
    axes[0].legend(loc="best", fontsize=6.5)

    # The wind grid's nodes that bracket the transfer's range, so that a range narrower than one
    # cell of the wind grid (a target at the anchor's own latitude) still has a cell to draw
    # (SPEC_05 Step 4, run 3b).
    inside = _bracketing(latitude, mesh_latitude.min(), mesh_latitude.max())
    within = _bracketing(pressure, labels.min(), labels.max())
    field = u_total[np.ix_(inside, within)]
    low, high = float(np.nanmin(field)), float(np.nanmax(field))
    if low == high:
        # A constant field has no range to contour, and the automatic levels would invent one at
        # round-off, a scale that looks like structure (SPEC_05 Step 3, the author's request).
        filled = axes[1].contourf(latitude[inside], pressure[within] / 100.0, field.T,
                                  levels=[low - 1.0, low + 1.0])
        bar = fig.colorbar(filled, ax=axes[1], label="u_total (m/s)")
        bar.set_ticks([low])
        axes[1].text(0.5, 0.5, f"u_total = {low:g} m/s everywhere", transform=axes[1].transAxes,
                     ha="center", va="center", fontsize=9, color="white")
    else:
        filled = axes[1].contourf(latitude[inside], pressure[within] / 100.0, field.T, levels=18)
        fig.colorbar(filled, ax=axes[1], label="u_total (m/s)")
    axes[1].set_yscale("log")
    axes[1].invert_yaxis()
    axes[1].axvline(target_deg, color=style.COLOR["height"], linestyle="--", linewidth=1.0)
    axes[1].axvline(gauge_deg, **style.ANCHOR_LINE)
    axes[1].set_xlabel("planetocentric latitude (degrees)")
    axes[1].set_ylabel("pressure (mbar)")
    axes[1].set_title(f"u_total, {str(wind.attrs.get('vertical_structure', 'unstated'))}")
    source = str(wind.attrs.get("source", "source unstated"))
    if len(source) > 78:
        source = source[:78].rsplit(" ", 1)[0] + " ..."
    fig.suptitle(f"F9. The wind assumed, {pr.slug}. {source}", fontsize=9)
    data = {"latitude_planetocentric_deg": latitude,
            "u_total_at_reference_ms": u_total_at_reference,
            "reference_level_pressure_Pa": reference_Pa,
            "pressure_Pa": pressure[within], "u_total_ms": field,
            "vertical_structure": str(wind.attrs.get("vertical_structure", "")),
            "source": str(wind.attrs.get("source", ""))}
    if u_reference is not None:
        data["u_reference_ms"] = u_reference
    return fig, data


FIGURES_TRANSFER = (
    ("F5", "geopotential", figure_5, None),
    ("F6", "pressure_identity", figure_6_transfer, "pressure_label_Pa"),
    ("F7", "delivered", figure_7, "pressure_label_Pa"),
    ("F8", "along_isobars", figure_8, "pressure_label_Pa"),
    ("F9", "wind_assumed", figure_9, "pressure_label_Pa"),
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
