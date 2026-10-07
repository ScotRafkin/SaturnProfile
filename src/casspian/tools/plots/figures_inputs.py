"""Panels of the input kinds, and the single figure views of a kind W, T or C file.

SPEC_02 Step 5. The panels here are shared: F1 of the kind N set draws the same temperature,
composition and mean property panels from the embedded inputs, and a new input file offered on
its own gets the same panels as its quick look before it enters a manifest.
"""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from casspian.lib import composition as comp
from casspian.lib import reduction as red
from casspian.tools.plots import style


# ---------------------------------------------------------------------------
# Shared panels
# ---------------------------------------------------------------------------


def panel_temperature(ax, thermo, anchor_Pa=None):
    """`T(p)`, with its uncertainty band where the file states one."""
    p = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
    T = np.asarray(thermo["temperature_K"].values, dtype="float64")
    dT = (np.asarray(thermo["temperature_uncertainty_K"].values, dtype="float64")
          if "temperature_uncertainty_K" in thermo.variables else np.full(p.shape, np.nan))
    finite = np.isfinite(dT)
    if finite.any():
        ax.fill_betweenx(p[finite] / 100.0, (T - dT)[finite], (T + dT)[finite],
                         color=style.COLOR["band"], alpha=0.6, linewidth=0,
                         label="+- uncertainty")
    ax.plot(T, p / 100.0, color=style.COLOR["temperature"], label="T, tabulated")
    style.pressure_axis(ax, p)
    if anchor_Pa is not None:
        style.anchor_line(ax, anchor_Pa)
    ax.set_xlabel("temperature (K)")
    note = "" if finite.any() else "; uncertainty not stated (NaN)"
    ax.set_title(f"temperature{note}")
    ax.legend(loc="lower left")
    return {"pressure_Pa": p, "temperature_K": T}


def _species_names(composition):
    return [str(v) for v in composition["species"].dataset["species_name"].values]


def _column(composition, latitude_deg):
    """The composition's column at `latitude_deg`, and the latitude actually drawn.

    A panel draws `x(p)`, so it needs one column of the field kind C now is (SPEC_04 decision
    O). A caller that knows the latitude passes it; a kind C file offered on its own has no
    latitude in its context, and the panel draws the node nearest the equator and says so.
    """
    if not hasattr(composition, "to_dataset"):
        return composition, latitude_deg
    nodes = comp.latitudes_of(composition)
    if latitude_deg is None:
        latitude_deg = float(nodes[int(np.argmin(np.abs(nodes)))])
    return comp.column_at(composition, latitude_deg), float(latitude_deg)


def panel_mole_fractions(ax, composition, anchor_Pa=None, latitude_deg=None):
    """`x_i(p)` on a log axis, the per level provenance of a flagged species by marker.

    Kind C is a field on `(level, latitude)` (SPEC_04 decision O), so the panel draws one
    column and `latitude_deg` says which, by the interpolant of decision L. The title names it,
    since a single column of a field that varies with latitude is not the whole file.
    """
    composition, latitude_deg = _column(composition, latitude_deg)
    root = composition.to_dataset(inherit=False) if hasattr(composition, "to_dataset") else composition
    p = np.asarray(root["pressure_Pa"].values, dtype="float64")
    names = _species_names(composition)
    data = {"pressure_Pa": p}
    for i, name in enumerate(names):
        x = np.asarray(root[f"x_{name}"].values, dtype="float64")
        positive = x > 0.0
        color = style.species_color(name, i)
        flag = f"{name.lower()}_provenance"
        ax.plot(np.where(positive, x, np.nan), p / 100.0, color=color, label=f"x_{name}")
        if flag in root.variables and positive.any():
            meanings = str(root[flag].attrs.get("flag_meanings", "")).split()
            codes = np.asarray(root[flag].values)
            for code, meaning in enumerate(meanings):
                chosen = (codes == code) & positive
                if chosen.any():
                    ax.plot(x[chosen], p[chosen] / 100.0, linestyle="none",
                            marker=style.PROVENANCE_MARKER.get(meaning, "."), color=color,
                            label=f"{name} {meaning}")
        data[f"x_{name}"] = x
    ax.set_xscale("log")
    style.pressure_axis(ax, p)
    if anchor_Pa is not None:
        style.anchor_line(ax, anchor_Pa, label=False)
    ax.set_xlabel("mole fraction (zero not drawn on the log axis)")
    ax.set_title("composition"
                 + ("" if latitude_deg is None else f" at {latitude_deg:.4f} deg"))
    # The trace species sit at the bottom left, so the legend goes where no species is drawn.
    ax.legend(loc="upper left", ncol=2)
    return data


def panel_mean_properties(ax, composition, anchor_Pa=None, mean_refractivity=None,
                          mean_molar_mass=None, latitude_deg=None):
    """`m_bar(p)` and `R_bar(p)` on twin axes.

    Taken from the product when it carries them; for a kind C file alone they are formed with
    `lib.reduction.mean_over_species`, which the panel title says. `latitude_deg` names the
    column, as in `panel_mole_fractions`.
    """
    composition, latitude_deg = _column(composition, latitude_deg)
    root = composition.to_dataset(inherit=False)
    species = composition["species"].to_dataset(inherit=False)
    p = np.asarray(root["pressure_Pa"].values, dtype="float64")
    source = "from the product"
    if mean_refractivity is None or mean_molar_mass is None:
        names = _species_names(composition)
        x = np.stack([np.asarray(root[f"x_{n}"].values, dtype="float64") for n in names])
        mean_refractivity = red.mean_over_species(
            x, np.asarray(species["refractivity_per_molecule_m3"].values, dtype="float64"))
        mean_molar_mass = red.mean_over_species(
            x, np.asarray(species["molar_mass_kg_mol"].values, dtype="float64"))
        source = "by lib.reduction from the composition"
    ax.plot(np.asarray(mean_molar_mass) * 1e3, p / 100.0, color=style.COLOR["molar_mass"],
            label="mean molar mass")
    ax.set_xlabel("mean molar mass (g/mol)", color=style.COLOR["molar_mass"])
    twin = ax.twiny()
    twin.plot(mean_refractivity, p / 100.0, color=style.COLOR["mean_refractivity"],
              linestyle="--", label="mean refractivity")
    twin.set_xlabel("mean refractivity per molecule (m3)", color=style.COLOR["mean_refractivity"])
    twin.grid(False)
    style.pressure_axis(ax, p)
    if anchor_Pa is not None:
        style.anchor_line(ax, anchor_Pa, label=False)
    ax.set_title(f"mean properties, {source}", pad=22)
    return {"pressure_Pa": p, "mean_molar_mass_kg_mol": np.asarray(mean_molar_mass),
            "mean_refractivity_m3": np.asarray(mean_refractivity)}


def wind_reference_column(wind):
    """Latitude (deg, ascending) and `u` on the reference level, with its uncertainty and flags."""
    lat = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    p = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    column = int(np.flatnonzero(p == float(wind["reference_level_pressure_Pa"]))[0])
    order = np.argsort(lat)
    u = np.asarray(wind["u_total_ms"].values, dtype="float64")[:, column][order]
    du = np.asarray(wind["u_total_uncertainty_ms"].values, dtype="float64")[:, column][order]
    flags = np.asarray(wind["value_provenance"].values)[:, column][order]
    return lat[order], u, du, flags


def wind_profile_at(wind, latitude_deg, name="u_total_ms"):
    """`name(p)` at one planetocentric latitude, linear in latitude on the file's grid."""
    lat = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    p = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    u = np.asarray(wind[name].values, dtype="float64")
    order = np.argsort(lat)
    profile = np.array([np.interp(float(latitude_deg), lat[order], u[order, j])
                        for j in range(p.size)])
    p_order = np.argsort(p)
    return p[p_order], profile[p_order]


def _reference_from_parts(wind):
    """`wind_reference_column` for a file whose total does not give the reference level.

    SPEC_08 section 6 ruling 10: the file's own `u_reference_ms`, with no uncertainty band (kind W
    gives the reference part none), and the provenance of the reference level where it is a node
    of the file's pressure grid.
    """
    lat = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    p = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    order = np.argsort(lat)
    u = np.asarray(wind["u_reference_ms"].values, dtype="float64")[order]
    column = np.flatnonzero(p == float(wind["reference_level_pressure_Pa"]))
    flags = (np.asarray(wind["value_provenance"].values)[:, int(column[0])][order]
             if column.size else np.full(lat.size, -1))
    return lat[order], u, np.full(lat.size, np.nan), flags


def panel_wind_latitude(ax, wind, latitude_deg=None):
    """`u(phi)` on the reference level, with the uncertainty band and provenance shading.

    The total at the reference level where the file carries the total on that level, and the
    file's reference wind otherwise (SPEC_08 section 6 ruling 10).
    """
    p = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    on_total = ("u_total_ms" in wind.variables
                and bool(np.any(p == float(wind["reference_level_pressure_Pa"]))))
    lat, u, du, flags = (wind_reference_column(wind) if on_total
                         else _reference_from_parts(wind))
    meanings = str(wind["value_provenance"].attrs.get("flag_meanings", "")).split()
    shade = {"parameterized": style.COLOR["parameterized"],
             "extrapolated": style.COLOR["extrapolated"],
             "extended_by_source_assumption": style.COLOR["extended"],
             "interpolated": style.COLOR["extended"]}
    for code, meaning in enumerate(meanings):
        if meaning not in shade:
            continue
        inside = flags == code
        if not inside.any():
            continue
        edges = np.flatnonzero(np.diff(np.concatenate([[0], inside.astype(np.int8), [0]])))
        label = meaning.replace("_", " ")
        for a, b in zip(edges[::2], edges[1::2]):
            ax.axvspan(lat[a], lat[min(b, lat.size) - 1], color=shade[meaning], zorder=0,
                       label=label)
            label = None
    finite = np.isfinite(du)
    if finite.any():
        ax.fill_between(lat[finite], (u - du)[finite], (u + du)[finite],
                        color=style.COLOR["band"], alpha=0.6, linewidth=0, label="+- 1 sigma")
    ax.plot(lat, u, color=style.COLOR["wind"],
            label=f"{'u' if on_total else 'u_reference'} at "
                  f"{float(wind['reference_level_pressure_Pa']) / 100.0:g} mbar")
    ax.axhline(0.0, color="0.75", linewidth=0.8, zorder=0)
    if latitude_deg is not None:
        style.latitude_line(ax, latitude_deg)
    ax.set_xlim(-90, 90)
    ax.set_xlabel("planetocentric latitude (deg)")
    ax.set_ylabel("zonal wind (m/s)")
    return {"latitude_deg": lat, "u_reference_ms": u, "u_reference_uncertainty_ms": du}


# ---------------------------------------------------------------------------
# Single figure views
# ---------------------------------------------------------------------------


def _figure_wind(wind):
    fig = Figure(figsize=style.FIGSIZE_WIDE)
    left, right = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1.6, 1.0]})
    style.layout(fig, rows=1)
    data = panel_wind_latitude(left, wind)
    left.set_title("u(phi) on the reference level, uncertainty band and provenance shading")
    left.legend(loc="upper right", fontsize=6.5)
    # SPEC_08 section 6 ruling 10: the right panel draws the total where the file carries it, the
    # shear where it carries only that, and is left blank where it carries neither.
    part = next((name for name in ("u_total_ms", "u_shear_ms") if name in wind.variables), None)
    p = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    if part is None:
        right.text(0.5, 0.5, "the file carries no u_total or u_shear", transform=right.transAxes,
                   ha="center", va="center", fontsize=9)
    else:
        marker = "o" if p.size == 1 else None
        for latitude, dash in ((0.0, "-"), (30.0, "--"), (-30.0, "-."), (60.0, ":"), (-60.0, ":")):
            p, u = wind_profile_at(wind, latitude, part)
            right.plot(u, p / 100.0, linestyle=dash, marker=marker, label=f"{latitude:+.0f} deg")
            key = "u_shear_profile" if part == "u_shear_ms" else "u_profile"
            data[f"{key}_{latitude:+.0f}_ms"] = u
        right.legend(loc="lower right")
    style.pressure_axis(right, p)
    right.set_xlabel("u_shear (m/s)" if part == "u_shear_ms" else "zonal wind (m/s)")
    what = "u_shear(p)" if part == "u_shear_ms" else "u(p)"
    right.set_title(f"{what} at five latitudes ({wind.attrs.get('vertical_structure', '')})")
    fig.suptitle(f"Wind: {wind.attrs.get('title', '')}")
    return fig, data


def _figure_thermo(thermo):
    fig = Figure(figsize=style.FIGSIZE_WIDE)
    ax = fig.subplots(1, 1)
    style.layout(fig, rows=1)
    data = panel_temperature(ax, thermo)
    label = thermo.attrs.get("latitude_planetographic_deg")
    where = f" at {float(label):g} deg planetographic" if label is not None else ""
    fig.suptitle(f"Thermodynamic profile: {thermo.attrs.get('title', '')}{where}")
    return fig, data


def _figure_composition(composition):
    fig = Figure(figsize=style.FIGSIZE_WIDE)
    left, right = fig.subplots(1, 2)
    style.layout(fig, rows=1)
    data = panel_mole_fractions(left, composition)
    data.update(panel_mean_properties(right, composition))
    fig.suptitle(f"Composition: {composition.attrs.get('title', '')}")
    return fig, data


def build(kind, handle):
    """The single figure of an input kind. Returns `(figures, skipped, data)`."""
    maker = {"wind": _figure_wind, "thermo": _figure_thermo,
             "composition": _figure_composition}[kind]
    fig, data = maker(handle)
    return [(kind, kind, fig)], {}, {kind: data}
