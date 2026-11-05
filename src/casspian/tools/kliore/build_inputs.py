"""Stage two of the Kliore tool: kinds T and D, and the reduction manifest.

SPEC_18 v0.3 section 5. Kind D and the manifest are the shared writers of `tools/profile_inputs.py`.
Kind T is Kliore's own, and it depends on the geoid, which is why this runs after kinds G, R and W:

* **The label latitude (ruling 1).** `refrac` reads a planetographic label and solves for `phi_c`;
  Kliore states 9.5 S planetocentric. The label is `phi_c + psi(phi_c)` on the geoid `refrac` marches
  (`refrac.anchor.geoid_setup`), so `refrac` returns 9.5 S.
* **Height (ruling 2).** `refrac` computes `r = r0 + (h - h_ref) cos psi`, reading `h` along the local
  vertical; Kliore's column is radial. `h = (r_local - r_local at 100 mbar) / cos psi`, so that kind
  N's radii keep his spacing. The datum is the 100 mbar level, which is the level `table1` inserted
  (SPEC_18 section 2), identified by its pressure.
* **Uncertainties (ruling 3)** are NaN here, as for Lindal. Kliore's ranges stay in the raw bundle.

Manuscript equations implemented: none. No physical value appears in this module.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.gravity import g_eff_vector
from casspian.refrac.anchor import geoid_setup
from casspian.tools import profile_inputs as pi

TOOL = "casspian-kliore-inputs"


def tilt(control, geodesy_path, latitude_deg: float) -> float:
    """`psi` in radians at the planetocentric `latitude_deg`, on the geoid `refrac` will march."""
    def read(path, kind):
        handle = cio.read(path, kind)
        try:
            return handle.load()
        finally:
            handle.close()

    inputs = SimpleNamespace(gravity=read(control["gravity_file"], "gravity"),
                             rotation=read(control["rotation_file"], "rotation"),
                             wind=read(control["wind_file"], "wind"),
                             geodesy=read(geodesy_path, "geodesy"))
    manifest = SimpleNamespace(anchor_surface_Pa=float(control["geoid_anchor_surface_Pa"]),
                               anchor_quantity=control["geoid_anchor_quantity"],
                               anchor_rule=control["geoid_anchor_rule"],
                               convergence_m=float(control["geoid_convergence_m"]))
    setup = geoid_setup(inputs, manifest)
    at = np.radians(np.array([latitude_deg]))
    radius = setup.march(at).radius
    return float(g_eff_vector(setup.u_of_phi(at), radius, at, *setup.constants)[3][0])


def build(control_path, section: str = "stage_two"):
    """Write kinds D and T and the reduction manifest. Returns the three paths."""
    control = pi.read_control(control_path, section)
    prefix = control["prefix"]
    control_path = Path(control_path)
    thermo_output = Path(control["thermo_output"])
    raw_path = Path(control["raw_bundle"])
    anchor_isobar = float(control["anchor_isobar_Pa"])

    raw = cio.read(raw_path, "raw")
    try:
        table1 = raw["table1"].dataset
        pressure = np.asarray(table1["pressure_Pa"].values, dtype="float64")
        temperature = np.asarray(table1["temperature_K"].values, dtype="float64")
        radius = np.asarray(table1["radius_local_m"].values, dtype="float64")
        source = dict(raw["scalars/source"].attrs)
        latitude = dict(raw["scalars/latitude"].attrs)
    finally:
        raw.close()

    pi.check_companions(control)
    geodesy_path = pi.write_geodesy(control, control_path, raw_path, TOOL)

    phi_c_deg = float(latitude["planetocentric_deg"])
    psi = tilt(control, geodesy_path, phi_c_deg)
    label_deg = phi_c_deg + float(np.degrees(psi))
    # The 100 mbar level, inserted at stage one (SPEC_18 section 2): the datum of the height.
    datum = int(np.flatnonzero(pressure == anchor_isobar)[0])
    height = (radius - radius[datum]) / np.cos(psi)

    nan = np.full(pressure.shape, np.nan)
    ranges = "not carried into kind T; Kliore's extreme ranges are in the raw bundle (SPEC_18 section 5)"
    thermo = xr.Dataset({
        "pressure_Pa": (("level",), pressure,
                        pi.attrs("Pa", "pressure of the source profile", "derived",
                                 value_source="Table 1, with the 100 mbar level inserted",
                                 positive="down", direction="increasing")),
        "pressure_printed_Pa": (("level",), pressure.copy(),
                                pi.attrs("Pa", "pressure of the source profile as printed",
                                         "derived",
                                         value_source="Table 1; the 100 mbar level is not printed")),
        "pressure_uncertainty_Pa": (("level",), nan.copy(),
                                    pi.attrs("Pa", "uncertainty on pressure", "derived",
                                             uncertainty_kind="stated", uncertainty_method=ranges)),
        "temperature_K": (("level",), temperature,
                          pi.attrs("K", "temperature of the source profile", "derived",
                                   value_source="Table 1, the 100 mbar level linear in ln p")),
        "temperature_uncertainty_K": (("level",), nan.copy(),
                                      pi.attrs("K", "uncertainty on temperature", "derived",
                                               uncertainty_kind="stated",
                                               uncertainty_method=ranges)),
        "height_m": (("level",), height,
                     pi.attrs("m", "altitude above the datum, along the local vertical", "measured",
                              value_source="the radial spacing of Table 1 over cos psi")),
        "height_uncertainty_m": (("level",), nan.copy(),
                                 pi.attrs("m", "uncertainty on height", "derived",
                                          uncertainty_kind="stated", uncertainty_method=ranges)),
    })
    thermo.attrs.update({
        "title": f"{prefix} thermodynamic profile, kind T",
        **pi.season_of(source),
        "source": source["citation"],
        "vertical_coordinate": "pressure_Pa",
        "latitude_planetocentric_absent_meaning": "point",
        "latitude_planetographic_deg": label_deg,
        "latitude_planetographic_deg_value_source": (
            f"computed from {phi_c_deg:g} deg planetocentric ({latitude['value_source']}) as "
            "phi_c + psi(phi_c) on the geoid refrac marches (SPEC_18 section 5 ruling 1)"),
        "latitude_definition": latitude["convention_source"],
        "height_datum": "the 100 mbar level, inserted at stage one (SPEC_18 section 2)",
        **pi.provenance(raw_path, control_path, thermo_output),
    })
    cio.history_append(
        thermo,
        f"{TOOL}: pressure and temperature from the raw bundle table1 group; height from its local "
        f"radius over cos psi, psi {np.degrees(psi):.5f} deg; the label latitude computed; "
        "uncertainties present and NaN; no composition and no geodesy in this file",
    )
    thermo_path = pi.write_unless_unchanged(thermo_output, thermo, "thermo", TOOL)
    manifest_path = pi.write_manifest(control, thermo_path, geodesy_path, TOOL)
    return thermo_path, geodesy_path, manifest_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Write the Kliore kind T and kind D files and the reduction manifest.")
    parser.add_argument("control", help="path to the TOML control file")
    parser.add_argument("--section", default="stage_two", help="section (default stage_two)")
    args = parser.parse_args(argv)
    for written in build(args.control, args.section):
        print(f"ready {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
