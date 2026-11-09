"""The three other Voyager profiles as anchors: the raw bundle (stage one) and kinds T and D (stage two).

SPEC_20 v0.3 section 2 deliverable 1, REPORT_20_preexecution section 1, REVIEW_20_step0. The source-
specific layer of SPEC_18's line for the profiles SPEC_15 traced from Lindal et al. (1985) Figs. 4, 5 and
6: Voyager 2 egress (31.2 S), Voyager 1 ingress (74.2 S) and Voyager 1 egress (3.0 S). Every tool from
the standard kinds on is the shared one, and Lindal's own files are read, never written.

**Stage one** (`casspian-voyager-raw`) reads a profile's scalars transcription and the traced CSV it
names, and writes the raw bundle:
* the rows in increasing pressure, each printed value converted once in exact decimal;
* one level at exactly 100 mbar, the anchor isobar, where the trace has none (the Voyager 2 egress),
  temperature and height linear in `ln p` between the neighbors, as Kliore's (SPEC_18 section 2);
* Lindal's ammonia (REVIEW_20_step0 ruling 2): his Table I values, linear in pressure onto this
  profile's levels inside the span he tabulates and blank outside it, so that kind C's declared rule
  (`tools/lindal/ammonia.py`) fills the rest exactly as it fills his.

**Stage two** (`casspian-voyager-inputs`) writes kind D under the profile's own name from the bundle's
transcription of Lindal's surfaces, then kind T and the manifest. Heights (ruling 1):
* where the trace carries Lindal's altitudes (the Voyager 2 egress), they are used as they are;
* otherwise, hydrostatic from `h = 0` at 100 mbar, `dh = -R T / (m_bar g) d ln p`, along the local
  vertical, with `m_bar` from the profile's own kind C and `g` the effective gravity of Lindal's G, R
  and W at `r = r0 + h cos psi`. `phi_c`, `r0` and `psi` are those `refrac` itself will find, from
  `refrac.anchor.freeze_anchor` on this profile's label and kind D, so the heights and the reduction
  agree. The integral is the trapezoid in `ln p`, repeated once with `g` at the first pass's heights.
`refrac` uses only `h - h_ref` (REPORT_20_preexecution section 1.2). sigma_T is NaN in kind T
(ruling 4); the traced bars stay in the raw bundle.

Manuscript equations implemented: none. No physical value appears in this module.
"""

from __future__ import annotations

import argparse
import csv
import math
import tomllib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.lib.gravity import g_eff_vector
from casspian.refrac.anchor import freeze_anchor
from casspian.tools import profile_inputs as pi
from casspian.tools.lindal.build_raw import SEASON_KEYS, _scaled, read_scalars
from casspian.tools.lindal.lindal_wind import mean_molar_mass

RAW_TOOL = "casspian-voyager-raw"
INPUTS_TOOL = "casspian-voyager-inputs"

#: The anchor isobar in Pa.
ANCHOR_ISOBAR_Pa = 1.0e4


# ---- stage one ------------------------------------------------------------------------------------

def read_trace(path: Path) -> dict:
    """The traced CSV in SI, increasing pressure. `height_m` only where the trace has altitudes."""
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.lstrip().startswith("#")))
    table = {
        "pressure_Pa": np.array([_scaled(r["pressure_mbar"], "100") for r in rows]),
        "temperature_K": np.array([_scaled(r["temperature_K"], "1") for r in rows]),
        "temperature_uncertainty_K": np.array([_scaled(r["temperature_uncertainty_K"] or "", "1") for r in rows]),
    }
    if "altitude_km" in rows[0]:
        table["height_m"] = np.array([_scaled(r["altitude_km"], "1000") for r in rows])
    order = np.argsort(table["pressure_Pa"])
    table = {k: v[order] for k, v in table.items()}
    if np.any(np.diff(table["pressure_Pa"]) <= 0):
        raise ValueError(f"{path}: pressures repeat")
    return table


def insert_level(table: dict, pressure_Pa: float) -> tuple[dict, bool]:
    """The table with a row at `pressure_Pa` if it has none: temperature and height linear in ln p,
    the traced bar NaN there. Returns the table and whether a row was inserted."""
    p = table["pressure_Pa"]
    if np.any(p == pressure_Pa):
        return table, False
    k = int(np.searchsorted(p, pressure_Pa))
    weight = math.log(pressure_Pa / p[k - 1]) / math.log(p[k] / p[k - 1])
    out = {}
    for name, column in table.items():
        if name == "pressure_Pa":
            value = pressure_Pa
        elif name in ("temperature_K", "height_m"):
            value = column[k - 1] + weight * (column[k] - column[k - 1])
        else:
            value = np.nan
        out[name] = np.insert(column, k, value)
    return out, True


def lindal_ammonia(bundle: Path, pressure_Pa) -> np.ndarray:
    """Lindal's Table I ammonia, linear in pressure onto `pressure_Pa` inside his tabulated span, NaN
    outside it, where kind C's declared rule fills (REVIEW_20_step0 ruling 2)."""
    raw = cio.read(bundle, "raw")
    try:
        table1 = raw["table1"].dataset
        p = np.asarray(table1["pressure_Pa"].values, dtype="float64")
        x = np.asarray(table1["nh3_mole_fraction"].values, dtype="float64")
    finally:
        raw.close()
    good = np.isfinite(x)
    p, x = p[good], x[good]
    order = np.argsort(p)
    p, x = p[order], x[order]
    at = np.asarray(pressure_Pa, dtype="float64")
    out = np.interp(at, p, x)
    out[(at < p[0]) | (at > p[-1])] = np.nan
    return out


def build_raw(scalars_toml) -> Path:
    """Build a profile's raw bundle from its scalars transcription. Returns its path."""
    scalars_toml = Path(scalars_toml).resolve()
    with open(scalars_toml, "rb") as handle:
        scalars = tomllib.load(handle)
    source, profile = scalars["source"], scalars["profile"]
    missing = [key for key in SEASON_KEYS if key not in source]
    if missing:
        raise ValueError(f"{scalars_toml}: [source] lacks {missing}")
    here = scalars_toml.parent
    trace_csv = (here / profile["trace"]).resolve()
    ammonia_bundle = (here / profile["ammonia_bundle"]).resolve()
    output = (here / profile["output"]).resolve()

    table, inserted = insert_level(read_trace(trace_csv), ANCHOR_ISOBAR_Pa)
    table["nh3_mole_fraction"] = lindal_ammonia(ammonia_bundle, table["pressure_Pa"])
    names = {
        "pressure_Pa": ("Pa", "pressure, as traced"),
        "temperature_K": ("K", "temperature, as traced"),
        "temperature_uncertainty_K": ("K", "half a drawn error bar, at the pressure where it is drawn; NaN elsewhere"),
        "height_m": ("m", "Lindal's altitude above 1 bar along the local vertical, as traced"),
        "nh3_mole_fraction": ("mol/mol", "Lindal's Table I ammonia, linear in pressure onto these levels "
                                         "inside his tabulated span; blank outside it"),
    }
    variables = {name: (("level",), table[name], {"units": units, "long_name": long_name,
                                                  "value_source": profile["value_source"]})
                 for name, (units, long_name) in names.items() if name in table}
    table1 = xr.Dataset(variables, attrs={"level_order": "increasing pressure",
                                         "inserted_level_Pa": ANCHOR_ISOBAR_Pa if inserted else np.nan})
    groups = {"table1": table1}
    groups.update(read_scalars(scalars_toml))
    data_inputs = [trace_csv, scalars_toml, ammonia_bundle]
    root = xr.Dataset(attrs={
        "title": f"{source['spacecraft']} {source['event']}, Lindal et al. (1985), raw bundle",
        "source": source["citation"],
        "epoch": str(source["observation_date"]),
        "solar_longitude_deg": float(source["solar_longitude_deg"]),
        "subsolar_latitude_deg": float(source["subsolar_latitude_deg"]),
        "solar_longitude_source": str(source["solar_longitude_source"]),
        "input_hashes": cio.input_hashes(data_inputs, output),
    })
    cio.history_append(root, f"{RAW_TOOL}: the trace in SI by exact decimal conversion, increasing pressure; "
                             + (f"one level inserted at {ANCHOR_ISOBAR_Pa:g} Pa, linear in ln p; " if inserted else "")
                             + "Lindal's Table I ammonia carried onto these levels")
    return cio.write(output, root, "raw", groups=groups, created_by=RAW_TOOL)


# ---- stage two ------------------------------------------------------------------------------------

def _read(path, kind):
    handle = cio.read(path, kind)
    try:
        return handle.load()
    finally:
        handle.close()


def hydrostatic_height(pressure, temperature, composition, anchor, constants) -> np.ndarray:
    """Height along the local vertical above the 100 mbar level, `dh = -R T / (m_bar g) d ln p`."""
    x = np.log(pressure)
    datum = int(np.flatnonzero(pressure == ANCHOR_ISOBAR_Pa)[0])
    phi_deg = float(np.degrees(anchor.phi_c_rad))
    m = mean_molar_mass(composition, pressure, np.array([phi_deg]))[:, 0]
    h = np.zeros_like(pressure)
    for _ in range(2):
        radius = anchor.r0_m + h * math.cos(anchor.psi_rad)
        g = np.asarray(g_eff_vector(np.full(h.shape, anchor.u_at_phi_c_ms), radius,
                                    np.full(h.shape, anchor.phi_c_rad), *constants)[2])
        rate = -MOLAR_GAS_CONSTANT * temperature / (m * g)
        h = np.concatenate([[0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(x))])
        h = h - h[datum]
    return h


def build_inputs(control_path, section: str = "stage_two"):
    """Write kinds D and T and the reduction manifest. Returns the three paths."""
    control = pi.read_control(control_path, section)
    prefix = control["prefix"]
    control_path = Path(control_path)
    thermo_output = Path(control["thermo_output"])
    raw_path = Path(control["raw_bundle"])
    if float(control["anchor_isobar_Pa"]) != ANCHOR_ISOBAR_Pa:
        raise ValueError(f"{control_path}: the anchor isobar is {ANCHOR_ISOBAR_Pa:g} Pa here (SPEC_20)")

    raw = cio.read(raw_path, "raw")
    try:
        table1 = raw["table1"].dataset
        pressure = np.asarray(table1["pressure_Pa"].values, dtype="float64")
        temperature = np.asarray(table1["temperature_K"].values, dtype="float64")
        traced_height = (np.asarray(table1["height_m"].values, dtype="float64")
                         if "height_m" in table1.variables else None)
        source = dict(raw["scalars/source"].attrs)
        latitude = dict(raw["scalars/latitude"].attrs)
        profile = dict(raw["scalars/profile"].attrs)
    finally:
        raw.close()

    pi.check_companions(control)
    geodesy_path = pi.write_geodesy(control, control_path, raw_path, INPUTS_TOOL)

    label = float(latitude["planetographic_deg"])
    if traced_height is not None:
        height = traced_height
        height_source = "Lindal's altitude above 1 bar along the local vertical, as traced (SPEC_15)"
        datum_note = "1 bar, Lindal's datum; refrac uses only the height above 100 mbar"
        method = "the traced altitudes"
    else:
        inputs = SimpleNamespace(gravity=_read(control["gravity_file"], "gravity"),
                                 rotation=_read(control["rotation_file"], "rotation"),
                                 wind=_read(control["wind_file"], "wind"),
                                 geodesy=_read(geodesy_path, "geodesy"),
                                 thermo=SimpleNamespace(attrs={"latitude_planetographic_deg": label}))
        manifest = SimpleNamespace(anchor_surface_Pa=float(control["geoid_anchor_surface_Pa"]),
                                   anchor_quantity=control["geoid_anchor_quantity"],
                                   anchor_rule=control["geoid_anchor_rule"],
                                   convergence_m=float(control["geoid_convergence_m"]),
                                   fixed_point_tolerance_deg=float(control["fixed_point_tolerance_deg"]),
                                   max_iterations=int(control["max_iterations"]))
        anchor = freeze_anchor(inputs, manifest)
        g, r = inputs.gravity, inputs.rotation
        constants = (float(r["angular_rate_rad_s"]), float(g["GM_m3s2"]), np.asarray(g["J"].values, dtype="float64"),
                     np.asarray(g["degree"].values), float(g["normalization_radius_m"]))
        composition = cio.read(control["composition_file"], "composition")
        try:
            height = hydrostatic_height(pressure, temperature, composition, anchor, constants)
        finally:
            composition.close()
        height_source = (f"hydrostatic from 100 mbar on the traced temperature, m_bar of this profile's kind C, "
                         f"g of Lindal's G, R and W at r0 {anchor.r0_m:.1f} m, psi {np.degrees(anchor.psi_rad):.5f} deg, "
                         f"phi_c {np.degrees(anchor.phi_c_rad):.5f} deg, as refrac finds them (REVIEW_20_step0 ruling 1)")
        datum_note = "the 100 mbar level"
        method = (f"heights hydrostatic from 100 mbar at phi_c {np.degrees(anchor.phi_c_rad):.5f} deg, "
                  f"r0 {anchor.r0_m:.1f} m")

    nan = np.full(pressure.shape, np.nan)
    none = "the trace's half bars are in the raw bundle; NaN here, as SPEC_18 (REVIEW_20_step0 ruling 4)"
    thermo = xr.Dataset({
        "pressure_Pa": (("level",), pressure,
                        pi.attrs("Pa", "pressure of the source profile", "derived", value_source=profile["value_source"],
                                 positive="down", direction="increasing")),
        "pressure_printed_Pa": (("level",), pressure.copy(),
                                pi.attrs("Pa", "pressure of the source profile as traced", "derived",
                                         value_source=profile["value_source"])),
        "pressure_uncertainty_Pa": (("level",), nan.copy(),
                                    pi.attrs("Pa", "uncertainty on pressure", "derived", uncertainty_kind="stated",
                                             uncertainty_method="the source states none")),
        "temperature_K": (("level",), temperature,
                          pi.attrs("K", "temperature of the source profile", "derived",
                                   value_source=profile["value_source"])),
        "temperature_uncertainty_K": (("level",), nan.copy(),
                                      pi.attrs("K", "uncertainty on temperature", "derived", uncertainty_kind="stated",
                                               uncertainty_method=none)),
        "height_m": (("level",), height,
                     pi.attrs("m", "altitude above the datum, along the local vertical", "derived",
                              value_source=height_source)),
        "height_uncertainty_m": (("level",), nan.copy(),
                                 pi.attrs("m", "uncertainty on height", "derived", uncertainty_kind="stated",
                                          uncertainty_method="the source states none")),
    })
    thermo.attrs.update({
        "title": f"{prefix} thermodynamic profile, kind T",
        **pi.season_of(source),
        "source": source["citation"],
        "vertical_coordinate": "pressure_Pa",
        "latitude_planetocentric_absent_meaning": "point",
        "latitude_planetographic_deg": label,
        "latitude_planetographic_deg_value_source": latitude["value_source"],
        "latitude_definition": latitude["convention_source"],
        "height_datum": datum_note,
        **pi.provenance(raw_path, control_path, thermo_output),
    })
    cio.history_append(thermo, f"{INPUTS_TOOL}: pressure and temperature from the raw bundle table1 group; "
                               f"{method}; uncertainties present and NaN; no composition and no geodesy in this file")
    thermo_path = pi.write_unless_unchanged(thermo_output, thermo, "thermo", INPUTS_TOOL)
    manifest_path = pi.write_manifest(control, thermo_path, geodesy_path, INPUTS_TOOL)
    return thermo_path, geodesy_path, manifest_path


def main_raw(argv=None) -> int:
    parser = argparse.ArgumentParser(prog=RAW_TOOL, description="Build a Voyager profile's raw bundle.")
    parser.add_argument("scalars", help="the profile's scalars transcription, occul_data/lindal/raw/<prefix>_scalars.toml")
    args = parser.parse_args(argv)
    print(f"wrote {build_raw(args.scalars)}")
    return 0


def main_inputs(argv=None) -> int:
    parser = argparse.ArgumentParser(prog=INPUTS_TOOL,
                                     description="Write a Voyager profile's kinds T and D and its reduction manifest.")
    parser.add_argument("control", help="the profile's build file")
    parser.add_argument("--section", default="stage_two", help="section (default stage_two)")
    args = parser.parse_args(argv)
    for written in build_inputs(args.control, args.section):
        print(f"ready {written}")
    return 0
