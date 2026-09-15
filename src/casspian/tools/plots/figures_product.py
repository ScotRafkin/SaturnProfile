"""The kind N figure set, F1 to F4. SPEC_02 Step 5.

Every figure is drawn from the product file and the inputs embedded in it. Where a panel needs a
quantity the file does not carry (gravity along the profile, the two geoids across latitude, the
recovered temperature), it calls `lib` on the embedded inputs and its title says so. Nothing
computed here is written to any product.

F5 and F6 draw forward model fields and belong to kind `profile` (`figures_profile`, SPEC_03 Step 3,
SPEC_02 v0.10 decision 10 superseded); kind N renders F1 to F4 and reports nothing skipped.
"""

from __future__ import annotations

import tomllib

import numpy as np
from matplotlib.figure import Figure

from casspian.lib import geoid as gd
from casspian.lib import reduction as red
from casspian.lib.gravity import G_phi_eff, g_eff_radial, g_eff_vector, g_newton
from casspian.tools.plots import figures_inputs as fi
from casspian.tools.plots import style

class _Product:
    """The product and its embedded inputs, read once."""

    def __init__(self, tree):
        self.root = tree.to_dataset(inherit=False)
        self.thermo = tree["inputs/thermo"].to_dataset(inherit=False)
        self.composition = tree["inputs/composition"]
        self.wind = tree["inputs/wind"].to_dataset(inherit=False)
        gravity = tree["inputs/gravity"].to_dataset(inherit=False)
        rotation = tree["inputs/rotation"].to_dataset(inherit=False)
        self.record = dict(tree["reduction_record"].attrs)
        self.manifest = tomllib.loads(str(tree["manifest"].attrs["text"]))
        self.slug = str(self.root.attrs.get("profile_or_run", ""))
        self.constants = (
            float(rotation["angular_rate_rad_s"]),
            float(gravity["GM_m3s2"]),
            np.asarray(gravity["J"].values, dtype="float64"),
            np.asarray(gravity["degree"].values),
            float(gravity["normalization_radius_m"]),
        )
        self.p = np.asarray(self.thermo["pressure_Pa"].values, dtype="float64")
        self.anchor_Pa = float(self.root["anchor_isobar_pressure_Pa"].values)
        self.phi_c_deg = float(self.root["latitude_planetocentric_deg"].values)
        self.phi_c = np.radians(self.phi_c_deg)
        self.r0 = float(self.root["anchor_isobar_radius_m"].values)
        self.radius = np.asarray(self.root["radius_m"].values, dtype="float64")

    def u_at_levels(self):
        """`u` at `phi_c` on the product's levels, linear in log pressure between wind levels."""
        p_w, u_w = fi.wind_profile_at(self.wind, self.phi_c_deg)
        return np.interp(np.log(self.p), np.log(p_w), u_w)

    def reference_wind(self):
        lat, u, _, _ = fi.wind_reference_column(self.wind)
        lat_rad = np.radians(lat)

        def u_of_phi(phi):
            phi = np.asarray(phi, dtype="float64")
            return np.reshape(np.interp(phi, lat_rad, u), phi.shape)

        return lat, u_of_phi


def _title(fig, text):
    fig.suptitle(text, fontsize=10)


def figure_1(pr: _Product):
    fig = Figure(figsize=style.FIGSIZE)
    (a, b), (c, d) = fig.subplots(2, 2)
    style.layout(fig)
    data = {}
    data.update(fi.panel_temperature(a, pr.thermo, pr.anchor_Pa))

    p_w, u_w = fi.wind_profile_at(pr.wind, pr.phi_c_deg)
    b.plot(u_w, p_w / 100.0, color=style.COLOR["wind"], label="u at phi_c, embedded kind W")
    style.pressure_axis(b, pr.p)
    style.anchor_line(b, pr.anchor_Pa)
    span = max(5.0, float(np.nanmax(np.abs(u_w))) * 1.2)
    b.set_xlim(-span, span)
    b.set_xlabel("zonal wind (m/s)")
    b.set_title(f"u(p) at phi_c ({pr.wind.attrs.get('vertical_structure', 'vertical structure')})")
    b.legend(loc="lower left")
    # Raised clear of the panel's own tick labels, so the inset's axis label sits inside.
    inset = b.inset_axes([0.58, 0.20, 0.39, 0.36])
    lat, u_ref, _, _ = fi.wind_reference_column(pr.wind)
    inset.plot(lat, u_ref, color=style.COLOR["wind"], linewidth=0.9)
    style.latitude_line(inset, pr.phi_c_deg, label=False)
    inset.set_xlim(-90, 90)
    inset.tick_params(labelsize=6)
    inset.set_title("u(phi), reference level", fontsize=6.5)
    inset.set_xlabel("planetocentric lat. (deg)", fontsize=6, labelpad=1)
    data["u_at_phi_c_ms"] = u_w

    data.update(fi.panel_mole_fractions(c, pr.composition, pr.anchor_Pa))
    data.update(fi.panel_mean_properties(
        d, pr.composition, pr.anchor_Pa,
        mean_refractivity=np.asarray(pr.root["mean_refractivity_m3"].values),
        mean_molar_mass=np.asarray(pr.root["mean_molar_mass_kg_mol"].values)))
    _title(fig, f"F1. Inputs of the {pr.slug} reduction, phi_c = {pr.phi_c_deg:.4f} deg")
    return fig, data


def figure_2(pr: _Product):
    fig = Figure(figsize=style.FIGSIZE)
    (a, b), (c, d) = fig.subplots(2, 2)
    style.layout(fig)
    u = pr.u_at_levels()
    r, phi = pr.radius, pr.phi_c
    g_n = g_newton(r, phi, *pr.constants[1:])
    g_rigid = g_eff_radial(0.0, r, phi, *pr.constants)
    g_eff = g_eff_radial(u, r, phi, *pr.constants)
    G_phi = G_phi_eff(u, r, phi, *pr.constants)
    psi = np.degrees(g_eff_vector(u, r, phi, *pr.constants)[3])
    p = pr.p / 100.0

    a.plot(g_n, p, color=style.COLOR["gravity_newton"], label="g_N")
    a.plot(g_eff, p, color=style.COLOR["gravity_effective"], linestyle="--", label="g_eff")
    a.set_xlabel("radial gravity, inward (m/s2)")
    a.set_title("g_N and g_eff (lib.gravity on the embedded G, R, W)")
    b.plot(g_n - g_rigid, p, color=style.COLOR["centrifugal"], label="rigid centrifugal part")
    b.set_xlabel("rigid centrifugal part of g_N - g_eff (m/s2)", color=style.COLOR["centrifugal"])
    b.set_title("the difference, rigid rotation and wind (lib.gravity)", pad=22)
    # The wind part is a few parts in 1e4 of the rigid part, so it gets its own scale; on a
    # shared axis it would lie on zero and say nothing.
    wind_axis = b.twiny()
    wind_axis.plot(g_rigid - g_eff, p, color=style.COLOR["wind_term"], linestyle="--",
                   label="wind part")
    wind_axis.set_xlabel("wind part (m/s2)", color=style.COLOR["wind_term"])
    wind_axis.grid(False)
    # From zero and without an offset: a nearly constant wind part would otherwise be drawn
    # as a steep curve across a range of parts in 1e10, which is round-off, not structure.
    wind_axis.ticklabel_format(axis="x", useOffset=False, style="sci", scilimits=(-3, 3))
    largest = float(np.nanmax(np.abs(g_rigid - g_eff)))
    if largest > 0.0:
        wind_axis.set_xlim(-0.05 * largest, 1.5 * largest)
    c.plot(G_phi, p, color=style.COLOR["G_phi"], label="G_phi")
    c.set_xlabel("latitudinal gravity, along increasing latitude (m/s2)")
    c.set_title("G_phi (lib.gravity)")
    d.plot(psi, p, color=style.COLOR["psi"], label="psi")
    d.set_xlabel("tilt of the local vertical, psi (deg)")
    d.set_title("psi (lib.gravity)")
    for ax, where in ((a, "center"), (b, "center"), (c, "upper left"), (d, "upper left")):
        style.pressure_axis(ax, pr.p)
        style.anchor_line(ax, pr.anchor_Pa, label=ax is a)
        handles, labels = ax.get_legend_handles_labels()
        if ax is b:
            extra = wind_axis.get_legend_handles_labels()
            handles, labels = handles + extra[0], labels + extra[1]
        ax.legend(handles, labels, loc=where)
    _title(fig, f"F2. Gravity along the {pr.slug} profile at phi_c = {pr.phi_c_deg:.4f} deg")
    return fig, {"pressure_Pa": pr.p, "radius_m": r, "u_ms": u, "g_newton_ms2": g_n,
                 "g_eff_ms2": g_eff, "g_rigid_ms2": g_rigid, "G_phi_ms2": G_phi, "psi_deg": psi}


def figure_3(pr: _Product):
    fig = Figure(figsize=style.FIGSIZE)
    (a, b), (c, d) = fig.subplots(2, 2)
    style.layout(fig)
    lat, u_of_phi = pr.reference_wind()
    phi = np.radians(lat)
    r_anchor = float(pr.record["anchor_radius_m"])
    rule = str(pr.record["anchor_rule"])
    tol_m = float(pr.manifest.get("geoid", {}).get("convergence_m", 1.0))
    # The no wind surface passes through the anchor radius where the rule anchors it.
    reference_latitude = 0.0 if rule == "equatorial_radius" else np.pi / 2
    r_ref, _, _ = gd.reference_geoid(phi, r_anchor, *pr.constants, tol_m=tol_m,
                                     anchor_latitude=reference_latitude)
    marched = gd.wind_geoid(phi, r_anchor, rule, u_of_phi, *pr.constants, tol_m=tol_m)
    r_wind = marched.radius
    u = u_of_phi(phi)
    g_n = g_newton(r_wind, phi, *pr.constants[1:])
    g_eff = g_eff_radial(u, r_wind, phi, *pr.constants)
    psi = np.degrees(g_eff_vector(u, r_wind, phi, *pr.constants)[3])

    a.plot(lat, g_n, color=style.COLOR["gravity_newton"], label="g_N")
    a.plot(lat, g_eff, color=style.COLOR["gravity_effective"], linestyle="--", label="g_eff")
    a.set_ylabel("radial gravity (m/s2)")
    a.set_title("g_N and g_eff on the wind geoid (lib.gravity, lib.geoid)")
    b.plot(lat, r_ref / 1e3, color=style.COLOR["geoid_nowind"], label="no wind reference geoid")
    b.plot(lat, r_wind / 1e3, color=style.COLOR["geoid_wind"], linestyle="--",
           label=f"wind geoid ({rule})")
    b.plot([pr.phi_c_deg], [pr.r0 / 1e3], marker="o", color=style.COLOR["anchor"],
           linestyle="none", label=f"r0 {pr.r0 / 1e3:.3f} km (from the file)")
    north = float(pr.record["polar_radius_north_m"]) / 1e3
    south = float(pr.record["polar_radius_south_m"]) / 1e3
    # Both polar radii annotated inside the axes, below the surface and clear of the legend.
    low = min(north, south)
    b.annotate(f"north pole {north:.2f} km", xy=(90, north), xytext=(22, low + 180),
               fontsize=6.5, arrowprops={"arrowstyle": "-", "color": "0.5"})
    b.annotate(f"south pole {south:.2f} km", xy=(-90, south), xytext=(-62, low + 180),
               fontsize=6.5, arrowprops={"arrowstyle": "-", "color": "0.5"})
    b.set_ylabel("radius of the anchor isobar (km)")
    b.set_title("the two geoids (lib.geoid)")
    c.plot(lat, (r_wind - r_ref) / 1e3, color=style.COLOR["dynamical_height"])
    c.set_ylabel("dynamical height, wind minus no wind (km)")
    c.set_title("dynamical height of the wind geoid")
    d.plot(lat, psi, color=style.COLOR["psi"], label="psi")
    d.set_ylabel("psi (deg)")
    d.set_title("tilt of the local vertical (lib.gravity)")
    for ax, where in ((a, "upper center"), (b, "upper left"), (c, None), (d, "lower right")):
        style.latitude_line(ax, pr.phi_c_deg, label=ax is a)
        ax.set_xlim(-90, 90)
        ax.set_xlabel("planetocentric latitude (deg)")
        if where is not None:
            ax.legend(loc=where, fontsize=6.5)
    _title(fig, f"F3. Gravity and shape across latitude on the anchor isobar "
                f"({pr.anchor_Pa / 100.0:g} mbar), {pr.slug}")
    return fig, {"latitude_deg": lat, "r_nowind_m": r_ref, "r_wind_m": r_wind,
                 "g_newton_ms2": g_n, "g_eff_ms2": g_eff, "psi_deg": psi,
                 "polar_north_m": marched.polar_north_m, "polar_south_m": marched.polar_south_m}


def _band(ax, value, uncertainty, p, color):
    finite = np.isfinite(uncertainty)
    if finite.any():
        ax.fill_betweenx(p[finite] / 100.0, (value - uncertainty)[finite],
                         (value + uncertainty)[finite], color=color, alpha=0.25, linewidth=0,
                         label="+- 1 sigma")
        return ""
    return "; uncertainty not stated (NaN)"


def figure_4(pr: _Product):
    fig = Figure(figsize=style.FIGSIZE)
    (a, b), (c, d) = fig.subplots(2, 2)
    style.layout(fig)
    root = pr.root
    n = np.asarray(root["number_density_m3"].values, dtype="float64")
    dn = np.asarray(root["number_density_uncertainty_m3"].values, dtype="float64")
    N = np.asarray(root["refractivity"].values, dtype="float64")
    dN = np.asarray(root["refractivity_uncertainty"].values, dtype="float64")
    R_bar = np.asarray(root["mean_refractivity_m3"].values, dtype="float64")
    T = np.asarray(pr.thermo["temperature_K"].values, dtype="float64")
    p = pr.p

    note = _band(a, n, dn, p, style.COLOR["number_density"])
    a.plot(n, p / 100.0, color=style.COLOR["number_density"], label="n")
    a.set_xscale("log")
    a.set_xlabel("number density (m-3)")
    a.set_title(f"number density{note}")
    note = _band(b, N, dN, p, style.COLOR["refractivity"])
    b.plot(N, p / 100.0, color=style.COLOR["refractivity"], label="N")
    b.set_xscale("log")
    b.set_xlabel("refractivity N (unscaled)")
    b.set_title(f"refractivity{note}", pad=22)
    # REVIEW_02_step5: a band of a few percent vanishes on a log axis, so the fractional
    # uncertainty is drawn as its own curve on a linear top axis, in percent, beside the band.
    fractional = 100.0 * dN / N
    percent_axis = b.twiny()
    percent_axis.plot(fractional, p / 100.0, color=style.COLOR["refractivity"], linestyle=":",
                      label="uncertainty / N (top axis)")
    percent_axis.set_xlabel("fractional uncertainty of N (percent)",
                            color=style.COLOR["refractivity"])
    percent_axis.grid(False)
    if np.any(np.isfinite(fractional)):
        percent_axis.set_xlim(0.0, 1.5 * float(np.nanmax(fractional)))
    for ax in (a, b):
        style.pressure_axis(ax, p)
        style.anchor_line(ax, pr.anchor_Pa, label=ax is a)
        handles, labels = ax.get_legend_handles_labels()
        if ax is b:
            extra = percent_axis.get_legend_handles_labels()
            handles, labels = handles + extra[0], labels + extra[1]
        ax.legend(handles, labels, loc="lower left")

    c.plot(N, pr.radius / 1e3, color=style.COLOR["refractivity"], label="N")
    c.axhline(pr.r0 / 1e3, **style.ANCHOR_LINE, label=f"anchor isobar, r0 {pr.r0 / 1e3:.3f} km")
    c.set_xscale("log")
    c.set_xlabel("refractivity N (unscaled)")
    c.set_ylabel("radius (km)")
    c.set_title("refractivity against absolute radius")
    c.legend(loc="upper right")

    recovered = red.temperature_from_refractivity(p, R_bar, N)
    fraction = recovered / T - 1.0
    worst = float(np.nanmax(np.abs(fraction)))
    limit = max(2.0e-12, 1.2 * worst)
    d.axvspan(-1.0e-12, 1.0e-12, color=style.COLOR["band"], alpha=0.4, linewidth=0,
              label="+- 1e-12, the acceptance bound")
    d.axvline(0.0, color="0.6", linewidth=0.8)
    d.plot(fraction, p / 100.0, marker="o", markersize=2.5, color=style.COLOR["recovered"],
           label="T recovered / T tabulated - 1")
    d.set_xlim(-limit, limit)
    d.set_xlabel("fractional difference")
    d.set_title(f"T = p R_bar / (k_B N) against the table (lib.reduction); max {worst:.1e}")
    style.pressure_axis(d, p)
    style.anchor_line(d, pr.anchor_Pa, label=False)
    d.legend(loc="lower left")
    _title(fig, f"F4. The {pr.slug} product")
    return fig, {"pressure_Pa": p, "recovered_temperature_fraction": fraction,
                 "recovered_temperature_max_abs": worst}


FIGURES = (
    ("F1", "inputs", figure_1),
    ("F2", "gravity_profile", figure_2),
    ("F3", "gravity_latitude", figure_3),
    ("F4", "product", figure_4),
)


def build_all(tree):
    """F1 to F4. Returns `(figures, skipped, data)`; nothing is skipped for kind N."""
    pr = _Product(tree)
    figures, data = [], {}
    for key, name, maker in FIGURES:
        fig, values = maker(pr)
        figures.append((key, name, fig))
        data[key] = values
    return figures, {}, data
