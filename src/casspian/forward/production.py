"""The forward production (manuscript Appendix B7) and `casspian-forward`. SPEC_03 Step 4.

`produce` is the production as a function that serves any latitude: given a refractivity profile on
its levels, the run's composition, gravity, rotation and wind, a gauge level and a boundary pressure,
it returns the geopotential, pressure and temperature. It is pure orchestration of `lib.geopotential`
(Step 1) and `lib.hydrostatic` (Step 2) on arrays and the run's loaded inputs, with no file access,
so that SPEC_04 calls it unchanged at a transferred latitude.

`run` is the closure: read the namelist, the anchor and the run's inputs, and check the inputs
against the anchor's embedded copies (`lib.control`, Step 3); take the boundary pressure by the
namelist's rule; produce on the anchor's own levels at its own latitude; write kind `profile` with
the anchor and the inputs embedded and the closure statistics in `production_record`; render F5 and
F6 when the namelist asks. The closure is a test of the production, not the model (SPEC_03 decision
6). `casspian-forward` takes the namelist path, dispatches on its mode, and calls `run`.

Manuscript equations implemented, through `lib`: B2 along the local vertical (the geopotential from
the effective gravity magnitude), B4 (density from refractivity), B5 (pressure by the exact
log-linear layer integral from the top), B6 (temperature).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

import numpy as np
import xarray as xr

from casspian.lib import composition as comp
from casspian.lib import control as ctl
from casspian.lib import geopotential as gp
from casspian.lib import hydrostatic as hs
from casspian.lib import io as cio
from casspian.lib import reduction as red
from casspian.lib.constants import CODATA_RELEASE
from casspian.lib.schema import PROFILE_EPOCH_WITHOUT_DATE, uncertainty_companion
from casspian.refrac.anchor import wind_of_latitude

__all__ = ["TOOL", "Profile", "Production", "closure_statistics", "main", "mean_properties",
           "produce", "profile_from_anchor", "run"]

TOOL = "casspian-forward"

#: The layer rules of Steps 1 and 2, named in `production_record`.
GEOPOTENTIAL_RULE = ("the trapezoid on the effective gravity magnitude |g_eff| along the local "
                     "vertical, dPhi = (|g_eff|_k + |g_eff|_(k+1)) / 2 (h_(k+1) - h_k), the sums "
                     "taken outward from the gauge level (lib.geopotential, SPEC_03 Step 1)")
HYDROSTATIC_RULE = ("the exact integral of the log-linear density across each layer, "
                    "rho_k dPhi expm1(x) / x with x = ln(rho_(k+1) / rho_k), pressure summed from "
                    "the boundary at the top (lib.hydrostatic, SPEC_03 Step 2)")

#: The pressure bins of the closure statistics (SPEC_03 v0.13 Step 4), in Pa. A level at a bin
#: edge belongs to the deeper bin. The 2 mbar edge is the grid level 10**2.3 Pa (199.526 Pa),
#: the row Table I prints as 2.00 mbar, so that row sits on the edge and the bin above 2 mbar
#: holds the ten top levels, as the specification counts them.
_P_2_MBAR, _P_10_MBAR, _P_100_MBAR = 10.0 ** 2.3, 1000.0, 10000.0


@dataclass(frozen=True)
class Profile:
    """A refractivity profile on its levels at one latitude. Angles in radians."""

    radius_m: np.ndarray
    height_above_anchor_isobar_m: np.ndarray
    refractivity: np.ndarray
    latitude_planetocentric_rad: float


@dataclass(frozen=True)
class Production:
    """What `produce` forms on the levels and layers."""

    u_ms: float
    geopotential: gp.GeopotentialProfile
    mean_refractivity_m3: np.ndarray
    mean_molar_mass_kg_mol: np.ndarray
    number_density_m3: np.ndarray
    density_kg_m3: np.ndarray
    layer_mass_Pa: np.ndarray
    pressure_Pa: np.ndarray
    temperature_K: np.ndarray
    gauge_level_index: int
    boundary_pressure_Pa: float


def profile_from_anchor(anchor) -> Profile:
    """The profile a kind N anchor carries: radius, field-line height, refractivity, `phi_c`."""
    root = anchor.to_dataset(inherit=False)
    return Profile(
        radius_m=np.asarray(root["radius_m"].values, dtype="float64"),
        height_above_anchor_isobar_m=np.asarray(root["height_above_anchor_isobar_m"].values,
                                                dtype="float64"),
        refractivity=np.asarray(root["refractivity"].values, dtype="float64"),
        latitude_planetocentric_rad=float(np.radians(float(root["latitude_planetocentric_deg"].values))),
    )


def mean_properties(composition, latitude_deg=None):
    """`(R_bar, m_bar)` on the composition's levels, `sum x_i R_i` and `sum x_i M_i`.

    From the kind C root's mole fractions and its `species` group, in the order the group lists the
    species; `m_bar` in kg/mol as kinds C and N carry it.

    Kind C is a field on `(level, latitude)` for every use (SPEC_04 decision O), so
    `latitude_deg` names the column to take, by the interpolant of decision L. Closure mode
    passes the anchor's own latitude; for a field uniform in latitude, which the Lindal
    hypothesis is, the interpolation returns the source's column unchanged.
    """
    if latitude_deg is not None:
        composition = comp.column_at(composition, latitude_deg)
    root = composition.to_dataset(inherit=False)
    species = composition["species"].to_dataset(inherit=False)
    names = [str(name) for name in species["species_name"].values]
    x = np.stack([np.asarray(root[f"x_{name}"].values, dtype="float64") for name in names])
    R_i = np.asarray(species["refractivity_per_molecule_m3"].values, dtype="float64")
    M_i = np.asarray(species["molar_mass_kg_mol"].values, dtype="float64")
    return red.mean_over_species(x, R_i), red.mean_over_species(x, M_i)


def produce(profile: Profile, inputs, gauge: int, p_b: float) -> Production:
    """The production on one column. Pure: arrays in, arrays out, no file access.

    `inputs` carries the run's loaded `composition` (a DataTree on the profile's levels),
    `gravity`, `rotation` and `wind`; `gauge` is the level index where `Phi = 0`; `p_b` the
    pressure at the top level. In order: `u(phi_c)` from the wind's reference level; `|g_eff|`,
    `psi` and `Phi` on the levels (Step 1); `R_bar` and `m_bar` from the composition; `n = N /
    R_bar` and `rho = n m_bar` (B4); the layer masses and `p` from the top (B5); `T` (B6).
    """
    gravity, rotation = inputs.gravity, inputs.rotation
    constants = (
        float(rotation["angular_rate_rad_s"]),
        float(gravity["GM_m3s2"]),
        np.asarray(gravity["J"].values, dtype="float64"),
        np.asarray(gravity["degree"].values),
        float(gravity["normalization_radius_m"]),
    )
    phi_c = profile.latitude_planetocentric_rad
    u = float(wind_of_latitude(inputs.wind)(np.array([phi_c]))[0])
    geopotential = gp.geopotential_along_profile(
        u, profile.radius_m, profile.height_above_anchor_isobar_m, phi_c, gauge, *constants)
    R_bar, m_bar = mean_properties(inputs.composition, np.degrees(phi_c))
    N = profile.refractivity
    n = N / R_bar
    rho = hs.density(N, R_bar, m_bar)
    layer_mass = hs.layer_mass(rho, geopotential.geopotential_m2s2)
    pressure = hs.pressure_from_top(p_b, layer_mass)
    temperature = hs.temperature(pressure, N, R_bar)
    return Production(
        u_ms=u,
        geopotential=geopotential,
        mean_refractivity_m3=R_bar,
        mean_molar_mass_kg_mol=m_bar,
        number_density_m3=n,
        density_kg_m3=rho,
        layer_mass_Pa=layer_mass,
        pressure_Pa=pressure,
        temperature_K=temperature,
        gauge_level_index=int(gauge),
        boundary_pressure_Pa=float(p_b),
    )


def closure_statistics(pressure, pressure_tabulated, temperature, temperature_tabulated) -> dict:
    """The closure residual and its statistics (SPEC_03 v0.13 Step 4), by pressure bin.

    The residual is `p / p_tab - 1` at every level. The bins, in tabulated pressure: above 2 mbar
    (the ten top levels), 2 to 10 mbar, 10 to 100 mbar, below 100 mbar excluding the bottom row;
    the bottom row (the deepest level) alone; and the mean below 10 mbar. A level on a bin edge
    belongs to the deeper bin. The largest residual from 2 mbar down excluding the bottom row is
    recorded too, since the acceptance bounds it.
    """
    p = np.asarray(pressure, dtype="float64")
    p_tab = np.asarray(pressure_tabulated, dtype="float64")
    residual = p / p_tab - 1.0
    residual_T = np.asarray(temperature, dtype="float64") / np.asarray(temperature_tabulated,
                                                                     dtype="float64") - 1.0
    bottom = int(np.argmax(p_tab))
    not_bottom = np.arange(p_tab.size) != bottom
    above_2 = p_tab < _P_2_MBAR
    from_2 = (p_tab >= _P_2_MBAR) & not_bottom
    b_2_10 = (p_tab >= _P_2_MBAR) & (p_tab < _P_10_MBAR)
    b_10_100 = (p_tab >= _P_10_MBAR) & (p_tab < _P_100_MBAR)
    below_100 = (p_tab >= _P_100_MBAR) & not_bottom
    below_10 = p_tab >= _P_10_MBAR

    def worst(mask):
        return float(np.max(np.abs(residual[mask]))) if mask.any() else float("nan")

    return {
        "residual_pressure": residual,
        "residual_temperature": residual_T,
        "residual_max_abs_above_2mbar": worst(above_2),
        "residual_max_abs_2_to_10mbar": worst(b_2_10),
        "residual_max_abs_10_to_100mbar": worst(b_10_100),
        "residual_max_abs_below_100mbar_excluding_bottom_row": worst(below_100),
        "residual_max_abs_from_2mbar_down_excluding_bottom_row": worst(from_2),
        "residual_bottom_row": float(residual[bottom]),
        "residual_mean_below_10mbar": float(np.mean(residual[below_10])),
        "residual_temperature_minus_pressure_max_abs": float(np.max(np.abs(residual_T - residual))),
        "levels_above_2mbar": int(above_2.sum()),
        "bottom_row_index": bottom,
        "bin_rule": ("residual p / p_tab - 1 by tabulated pressure (SPEC_03 v0.13): above 2 mbar "
                     "p < 10**2.3 Pa (the grid level Table I prints as 2.00 mbar, on the edge); "
                     "2 to 10 mbar 10**2.3 <= p < 1000 Pa; 10 to 100 mbar "
                     "1000 <= p < 10000 Pa; below 100 mbar p >= 10000 Pa, the bottom row excluded; "
                     "the bottom row is the deepest level; the mean below 10 mbar over "
                     "p >= 1000 Pa; a level on an edge belongs to the deeper bin"),
    }


# ---------------------------------------------------------------------------
# The closure run and the product
# ---------------------------------------------------------------------------


def _attrs(units, long_name, provenance, **extra):
    out = {"units": units, "long_name": long_name, "provenance": provenance}
    out.update({k: v for k, v in extra.items() if v is not None})
    return out


def _declared_terms(anchor_root) -> list:
    """The uncertainty terms the anchor declares, from its companions' term lists."""
    terms = []
    for name in ("refractivity_uncertainty", "radius_uncertainty_m"):
        if name in anchor_root.variables:
            text = str(anchor_root[name].attrs.get("uncertainty_terms_included", ""))
            for term in text.replace(";", ",").replace(" ", ",").split(","):
                if term and term not in terms:
                    terms.append(term)
    return terms


def _resolved_namelist(namelist, inputs, product: Path) -> str:
    """The namelist as the run resolved it, every path relative to the product and hashed."""
    lines = [
        f"# The run namelist as {TOOL} resolved it; paths relative to {product.name}, with hashes",
        "[run]",
        f'name = "{namelist.name}"',
        f'mode = "{namelist.mode}"',
        f"solar_longitude_deg = {namelist.solar_longitude_deg!r}",
        f"date = {('\"' + namelist.date + '\"') if namelist.date else '\"none\"'}",
        f'namelist = "{cio.input_hash_entry(namelist.path, product)}"',
        "[anchor]",
        f'slug = "{namelist.anchors[0].slug}"',
        f'path = "{cio.input_hash_entry(namelist.anchors[0].path, product)}"',
        "[inputs]",
    ]
    lines += [f'{key} = "{cio.input_hash_entry(path, product)}"' for key, path in namelist.inputs.items()]
    lines += [
        "[hydrostatic_boundary]",
        f'p_b_rule = "{namelist.p_b_rule}"',
        f'p_b_location = "{namelist.p_b_location}"',
        "[isobars]",
        f"gauge_isobar_Pa = {namelist.gauge_isobar_Pa!r}",
        f"gauge_level_index = {inputs.gauge_level_index}",
        "[diagnostics]",
    ]
    lines += [f"{key} = {value!r}" for key, value in dict(namelist.diagnostics).items()]
    return "\n".join(lines) + "\n"


def _nodes(tree, prefix: str) -> dict:
    out = {}
    base = tree.path.rstrip("/")
    for node in tree.subtree:
        relative = node.path[len(base):] if base else ("" if node.path == "/" else node.path)
        out[prefix + relative] = node.to_dataset(inherit=False)
    return out


@dataclass(frozen=True)
class RunResult:
    product: Path
    production: Production
    statistics: MappingProxyType
    figures: object


def run(namelist_path) -> RunResult:
    """The closure run of SPEC_03 Step 4. Returns the product path, the production, the statistics
    and the render result (None when the namelist asks for no figures)."""
    namelist = ctl.read_run_namelist(namelist_path)
    if namelist.mode != "closure":
        raise ctl.ControlFileError(f"{namelist.path}: mode {namelist.mode!r} is not implemented.")
    inputs = ctl.load_run_inputs(namelist)
    anchor = inputs.anchor
    anchor_root = anchor.to_dataset(inherit=False)
    thermo = anchor["inputs/thermo"].to_dataset(inherit=False)
    p_tab = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
    T_tab = np.asarray(thermo["temperature_K"].values, dtype="float64")

    if namelist.p_b_rule != "anchor_profile_top":
        raise ctl.ControlFileError(f"{namelist.path}: p_b_rule {namelist.p_b_rule!r} is not implemented.")
    boundary = int(np.argmin(p_tab))
    if boundary != 0:
        raise ctl.ControlFileError(
            f"the anchor's thermo levels are not top down (the least pressure is level {boundary}); "
            "the production integrates from the top level 0."
        )
    p_b = float(p_tab[boundary])

    profile = profile_from_anchor(anchor)
    production = produce(profile, inputs, inputs.gauge_level_index, p_b)
    statistics = closure_statistics(production.pressure_Pa, p_tab, production.temperature_K, T_tab)

    product = namelist.product
    record = _production_record(namelist, inputs, anchor, production, statistics)
    dataset = _product_dataset(namelist, inputs, anchor_root, production, p_tab, T_tab, product)
    groups = _nodes(anchor, f"anchors/{namelist.anchors[0].slug}")
    for key in ctl.RUN_INPUT_KINDS:
        loaded = getattr(inputs, key)
        if isinstance(loaded, xr.DataTree):
            groups.update(_nodes(loaded, f"inputs/{key}"))
        else:
            groups[f"inputs/{key}"] = loaded
    groups["namelist"] = xr.Dataset(attrs={
        "text": namelist.text,
        "sha256": namelist.sha256,
        "resolved": _resolved_namelist(namelist, inputs, product),
    })
    groups["production_record"] = xr.Dataset(attrs=record)
    written = cio.write(product, dataset, "profile", groups=groups, created_by=TOOL)

    figures = None
    diagnostics = dict(namelist.diagnostics)
    if diagnostics.get("figures"):
        # Drawn from the file just written, never from this run's memory (SPEC_02 Step 5).
        from casspian.tools.plots import render

        figures = render(written, namelist.output_directory / "figures", diagnostics["format"],
                         int(diagnostics["dpi"]))
    return RunResult(product=written, production=production,
                     statistics=MappingProxyType(statistics), figures=figures)


def _product_dataset(namelist, inputs, anchor_root, production, p_tab, T_tab, product):
    terms = _declared_terms(anchor_root)
    unstated = ", ".join(terms)
    reason = ("not propagated: a first order companion carrying only the terms that pass through "
              "the integral would be read as a total; the Monte Carlo wrapper supplies the product's "
              "uncertainty (SPEC_03 decision 5)")

    def companion(units, of):
        return (("level",), np.full(p_tab.size, np.nan),
                _attrs(units, f"uncertainty on {of}", "modeled", uncertainty_kind="1sigma",
                       uncertainty_method="not propagated", uncertainty_terms_unstated=unstated,
                       uncertainty_note=reason))

    level = ("level",)
    copied = "copied from the anchor"
    variables = {
        "radius_m": (level, np.asarray(anchor_root["radius_m"].values, dtype="float64"),
                     _attrs("m", f"absolute planetocentric radius of the level, {copied}", "derived",
                            positive="up")),
        "height_above_anchor_isobar_m": (
            level, np.asarray(anchor_root["height_above_anchor_isobar_m"].values, dtype="float64"),
            _attrs("m", f"height above the anchor isobar along the local vertical, {copied}",
                   "measured", positive="up")),
        "refractivity": (level, np.asarray(anchor_root["refractivity"].values, dtype="float64"),
                         _attrs("1", f"refractivity, unscaled, {copied}", "derived")),
        "mean_refractivity_m3": (level, production.mean_refractivity_m3,
                                 _attrs("m3", "mean refractivity per molecule from the run's kind C",
                                        "modeled")),
        "mean_refractivity_uncertainty_m3": companion("m3", "mean_refractivity_m3"),
        "mean_molar_mass_kg_mol": (level, production.mean_molar_mass_kg_mol,
                                   _attrs("kg mol-1", "mean molar mass from the run's kind C", "modeled")),
        "mean_molar_mass_uncertainty_kg_mol": companion("kg mol-1", "mean_molar_mass_kg_mol"),
        "number_density_m3": (level, production.number_density_m3,
                              _attrs("m-3", "number density, N / R_bar (Eq. B4)", "modeled")),
        "number_density_uncertainty_m3": companion("m-3", "number_density_m3"),
        "pressure_Pa": (level, production.pressure_Pa,
                        _attrs("Pa", "hydrostatic pressure, summed from the boundary at the top "
                               "(Eq. B5)", "modeled", boundary_pressure_Pa=production.boundary_pressure_Pa,
                               positive="down", direction="increasing")),
        "pressure_uncertainty_Pa": companion("Pa", "pressure_Pa"),
        "temperature_K": (level, production.temperature_K,
                          _attrs("K", "temperature, p R_bar / (k_B N) (Eq. B6)", "modeled")),
        "temperature_uncertainty_K": companion("K", "temperature_K"),
        "geopotential_uncertainty_m2s2": companion("m2 s-2", "geopotential_m2s2"),
        "pressure_tabulated_Pa": (level, p_tab,
                                  _attrs("Pa", "the anchor's tabulated pressure, for the closure",
                                         "derived")),
        "temperature_tabulated_K": (level, T_tab,
                                    _attrs("K", "the anchor's tabulated temperature, for the closure",
                                           "derived")),
        "latitude_planetocentric_deg": ((), float(anchor_root["latitude_planetocentric_deg"].values),
                                        _attrs("degrees_north", f"planetocentric latitude, {copied}",
                                               "derived")),
        "psi_deg": ((), float(anchor_root["psi_deg"].values),
                    _attrs("degrees", f"tilt of the local vertical at the anchor, {copied}", "derived")),
        "gauge_isobar_Pa": ((), float(namelist.gauge_isobar_Pa),
                            _attrs("Pa", "the declared gauge isobar, where Phi = 0", "index")),
        "gauge_level_index": ((), np.int32(production.gauge_level_index),
                              _attrs("1", "the level of the gauge isobar", "index")),
        "boundary_pressure_Pa": ((), production.boundary_pressure_Pa,
                                 _attrs("Pa", "the boundary pressure p_b, the anchor's top tabulated "
                                        "pressure by p_b_rule", "index")),
        "boundary_level_index": ((), np.int32(0), _attrs("1", "the level of the boundary", "index")),
    }
    dataset = xr.Dataset(
        variables,
        coords={"geopotential_m2s2": (level, production.geopotential.geopotential_m2s2,
                                      _attrs("m2 s-2", "geopotential along the local vertical, zero on "
                                             "the gauge isobar", "modeled", positive="up",
                                             gauge_isobar_Pa=float(namelist.gauge_isobar_Pa)))},
    )
    anchor = inputs.anchor
    input_paths = [namelist.anchors[0].path] + [namelist.inputs[k] for k in ctl.RUN_INPUT_KINDS] + [namelist.path]
    dataset.attrs = {
        "title": f"{namelist.name} forward profile, kind profile",
        "profile_or_run": namelist.name,
        "role": "forward",
        "source": (f"{TOOL} closure of the {namelist.anchors[0].slug} anchor under the run "
                   f"{namelist.name}; {namelist.description}"),
        "epoch": namelist.date if namelist.date else PROFILE_EPOCH_WITHOUT_DATE,
        "solar_longitude_deg": float(namelist.solar_longitude_deg),
        "solar_longitude_source": "declared in the namelist",
        "anchor_solar_longitudes_deg": np.array([float(anchor.attrs["solar_longitude_deg"])]),
        "boundary_pressure_Pa": production.boundary_pressure_Pa,
        "latitude_planetocentric_absent_meaning": "point",
        "mode": namelist.mode,
        "input_hashes": cio.input_hashes(input_paths, product),
    }
    cio.history_append(
        dataset,
        f"{TOOL}: closure of {namelist.anchors[0].slug} on its {p_tab.size} levels at phi_c; "
        f"geopotential by {GEOPOTENTIAL_RULE.split(' (')[0]}; pressure by "
        f"{HYDROSTATIC_RULE.split(' (')[0]}; boundary p_b = {production.boundary_pressure_Pa!r} Pa "
        f"by {namelist.p_b_rule}; gauge at level {production.gauge_level_index}",
    )
    return dataset


def _production_record(namelist, inputs, anchor, production, statistics) -> dict:
    anchor_record = dict(anchor["reduction_record"].attrs)
    record = {
        "mode": namelist.mode,
        "run": namelist.name,
        "anchor_slug": namelist.anchors[0].slug,
        "gauge_isobar_Pa": float(namelist.gauge_isobar_Pa),
        "gauge_level_index": np.int32(production.gauge_level_index),
        "p_b_rule": str(namelist.p_b_rule),
        "p_b_location": namelist.p_b_location,
        "boundary_pressure_Pa": production.boundary_pressure_Pa,
        "boundary_level_index": np.int32(0),
        "u_at_phi_c_ms": production.u_ms,
        "geopotential_rule": GEOPOTENTIAL_RULE,
        "hydrostatic_rule": HYDROSTATIC_RULE,
        "closure_comparison_rule": (
            "each run input compared with the anchor's embedded copy, root and every group, after "
            f"dropping {sorted(ctl.CLOSURE_DROPPED_ATTRIBUTES)} and any attribute holding sha256: "
            "(SPEC_03 Step 3 deliverable 2)"),
        "codata_release": CODATA_RELEASE,
        "casspian_version": cio.package_version(),
        "casspian_git_commit": cio.git_commit(),
    }
    for comparison in inputs.closure:
        record[f"closure_comparison_{comparison.kind}"] = "identical" if comparison.identical else \
            "; ".join(comparison.differences)
        record[f"closure_dropped_{comparison.kind}"] = "; ".join(
            f"{group}: {', '.join(names)}" for group, names in comparison.dropped.items())
    for key in ("radius_projection_rule", "radius_projection_residual_m", "latitude_drift_neglected_deg",
                "radius_projection_note"):
        if key in anchor_record:
            record[f"anchor_{key}"] = anchor_record[key]
    for key, value in statistics.items():
        record[key] = value if isinstance(value, (str, np.ndarray)) else (
            np.int32(value) if isinstance(value, int) else float(value))
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Run a forward model namelist (SPEC_03: closure mode).")
    parser.add_argument("namelist", help="path to <run>.toml")
    args = parser.parse_args(argv)
    result = run(args.namelist)
    print(f"wrote {result.product}")
    if result.figures is not None:
        from casspian.tools.plots import describe

        for line in describe(result.figures):
            print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
