"""Stage two of the Lindal tool: kinds T and D, and the reduction manifest.

SPEC_01 Step 9. Completes `occul_data/lindal/` so that `refrac` (SPEC_02) can run. Reads only
the raw bundle and the control file, verifies that the four files the other tools wrote are
present and valid, and writes the two files nobody else writes plus the manifest that names
all six.

Manuscript equations implemented: none. This tool separates the raw bundle's transcription
into the two standard kinds and records what the reduction is to be run with.

**What separating means here.** The raw bundle deliberately holds the thermodynamic profile and
the fitted geoid together, because that is how the paper states them. SPEC_00 section 6.3 keeps
them apart in the model inputs, because the retrieval and the geoid fit are different
measurements with different error budgets. This tool is the boundary: kind T carries the profile
and no geodesy, kind D carries the surfaces and no profile, and each says which raw bundle it
came from.

**What is Lindal's here (SPEC_18 v0.3 section 5).** Kind T's reading of Table I and its pressure
grid. Kind D and the manifest are general and live in `tools/profile_inputs.py`.

No physical value appears in this module. Every number comes from the raw bundle; every choice
comes from the control file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.control import ControlFileError
from casspian.tools import profile_inputs as pi

TOOL = "casspian-lindal-inputs"


def build(control_path, section: str = "stage_two"):
    """Write kinds T and D and the reduction manifest. Returns the three paths."""
    control = pi.read_control(control_path, section)
    prefix = control["prefix"]
    control_path = Path(control_path)
    thermo_output = Path(control["thermo_output"])

    raw_path = Path(control["raw_bundle"])
    raw = cio.read(raw_path, "raw")
    try:
        table1 = raw["table1"].dataset
        try:
            pressure_printed = np.asarray(table1["pressure_printed_Pa"].values, dtype="float64")
            pressure_grid = dict(raw["scalars/pressure_grid"].attrs)
        except KeyError as exc:
            raise ControlFileError(
                f"{raw_path}: the raw bundle carries no pressure_printed_Pa or no "
                "scalars/pressure_grid group; SPEC_01 v0.23 Step 9 needs both. Rebuild the raw "
                "bundle with casspian-lindal-raw."
            ) from exc
        pressure = np.asarray(table1["pressure_Pa"].values, dtype="float64")
        temperature = np.asarray(table1["temperature_K"].values, dtype="float64")
        height = np.asarray(table1["height_m"].values, dtype="float64")
        source = dict(raw["scalars/source"].attrs)
        latitude = dict(raw["scalars/latitude"].attrs)
        geodesy = dict(raw["scalars/geodesy"].attrs)
        wind_scalars = dict(raw["scalars/wind_used_by_source"].attrs)
        gravity_scalars = dict(raw["scalars/gravity_used_by_source"].attrs)
    finally:
        raw.close()

    pi.check_companions(control)

    # ---- kind T ----------------------------------------------------------------------
    nan = np.full(pressure.shape, np.nan)
    longitudes = np.asarray(latitude["longitude_system_iii_deg"], dtype="float64")
    thermo = xr.Dataset(
        {
            "pressure_Pa": (("level",), pressure,
                            pi.attrs("Pa", "pressure of the source profile, on the declared grid",
                                     "derived",
                                     value_source="Table I, placed on the grid in "
                                                  "pressure_grid_rule",
                                     positive="down", direction="increasing")),
            "pressure_printed_Pa": (("level",), pressure_printed,
                                    pi.attrs("Pa", "pressure of the source profile as printed",
                                             "derived", value_source="Table I")),
            "pressure_uncertainty_Pa": (("level",), nan.copy(),
                                        pi.attrs("Pa", "uncertainty on pressure", "derived",
                                                 uncertainty_kind="stated",
                                                 uncertainty_method="the source states none")),
            "temperature_K": (("level",), temperature,
                              pi.attrs("K", "temperature of the source profile", "derived",
                                       value_source="Table I")),
            "temperature_uncertainty_K": (("level",), nan.copy(),
                                          pi.attrs("K", "uncertainty on temperature", "derived",
                                                   uncertainty_kind="stated",
                                                   uncertainty_method="the source states none")),
            "height_m": (("level",), height,
                         pi.attrs("m", "altitude above the datum, the measured quantity",
                                  "measured", value_source="Table I")),
            "height_uncertainty_m": (("level",), nan.copy(),
                                     pi.attrs("m", "uncertainty on height", "derived",
                                              uncertainty_kind="stated",
                                              uncertainty_method="the source states none")),
        }
    )
    thermo.attrs.update({
        "title": f"{prefix} thermodynamic profile, kind T",
        **pi.season_of(source),
        "source": source["citation"],
        "vertical_coordinate": "pressure_Pa",
        "latitude_planetocentric_absent_meaning": "point",
        "latitude_planetographic_deg": float(latitude["planetographic_deg"]),
        "latitude_planetographic_deg_value_source": latitude["value_source"],
        "latitude_planetographic_deg_swath": np.asarray(
            latitude["swath_planetographic_deg"], dtype="float64"),
        "latitude_planetographic_deg_swath_source": latitude["swath_source"],
        "latitude_planetographic_deg_uncertainty": float(latitude["uncertainty_deg"]),
        "latitude_planetographic_deg_uncertainty_kind": latitude["uncertainty_kind"],
        "latitude_definition": latitude["convention_source"],
        "longitude_swath_deg": longitudes,
        "longitude_source": latitude["longitude_source"],
        # SPEC_01 v0.23 Step 9: the grid declaration of the transcription, quoted.
        "pressure_grid_rule": (
            f"{pressure_grid['rule']}; rows printed as "
            f"{np.asarray(pressure_grid['excluded_printed_values_mbar']).tolist()} mbar are "
            f"excluded and kept as printed; basis: {pressure_grid['basis']}; value_source: "
            f"{pressure_grid['value_source']}"),
        # From the dedicated table the Step 9 review added to the transcription. It used to be
        # cut from the prose of `fit_inputs`, which is still carried below as the source's own
        # statement of everything that went into the geoid fit.
        "source_gravity_citation_value_source": gravity_scalars["value_source"],
        "source_pole_vector": gravity_scalars["pole_vector"],
        "source_fit_inputs": geodesy["fit_inputs"],
        "source_wind_vertical_structure": wind_scalars["vertical_structure"],
        **pi.provenance(raw_path, control_path, thermo_output),
    })
    cio.history_append(
        thermo,
        f"{TOOL}: pressure (on the declared grid, with the printed value beside it), temperature "
        "and height taken from the raw bundle table1 group; uncertainties present and NaN, the "
        "source states none; no composition and no geodesy in this file (SPEC_00 section 6.1)",
    )
    thermo_path = pi.write_unless_unchanged(thermo_output, thermo, "thermo", TOOL)

    geodesy_path = pi.write_geodesy(control, control_path, raw_path, TOOL)
    manifest_path = pi.write_manifest(control, thermo_path, geodesy_path, TOOL)
    return thermo_path, geodesy_path, manifest_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL,
        description="Write the Lindal kind T and kind D files and the reduction manifest.",
    )
    parser.add_argument("control", help="path to the TOML control file")
    parser.add_argument("--section", default="stage_two", help="section (default stage_two)")
    args = parser.parse_args(argv)
    for written in build(args.control, args.section):
        print(f"ready {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
