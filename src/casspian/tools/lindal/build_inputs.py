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

No physical value appears in this module. Every number comes from the raw bundle; every choice
comes from the control file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.tools.gravity.control import ControlFileError, load_section

TOOL = "casspian-lindal-inputs"

SECTION_KEYS = {
    "raw_bundle": True,
    "prefix": True,
    "description": True,
    "thermo_output": True,
    "geodesy_output": True,
    "manifest": True,
    "product": True,
    "gravity_file": True,
    "rotation_file": True,
    "wind_file": True,
    "composition_file": True,
    "anchor_isobar_Pa": True,
    "fixed_point_tolerance_deg": True,
    "max_iterations": True,
    "geoid_anchor_surface_Pa": True,
    "geoid_anchor_quantity": True,
    "geoid_anchor_rule": True,
    "geoid_convergence_m": True,
}

PATH_KEYS = (
    "raw_bundle", "thermo_output", "geodesy_output", "manifest",
    "gravity_file", "rotation_file", "wind_file", "composition_file",
)

#: The four files the other tools wrote, with the kind each must validate as.
COMPANIONS = (
    ("gravity_file", "gravity"),
    ("rotation_file", "rotation"),
    ("wind_file", "wind"),
    ("composition_file", "composition"),
)


def _attrs(units, long_name, provenance, **extra):
    out = {"units": units, "long_name": long_name, "provenance": provenance}
    for key, value in extra.items():
        if value is None or (isinstance(value, str) and not value):
            continue
        out[key] = value
    return out


def _repository_root(start: Path) -> Path:
    start = Path(start).resolve()
    for candidate in [start, *start.parents]:
        if (candidate / ".git").exists():
            return candidate
    return start.parent


def build(control_path, section: str = "stage_two"):
    """Write kinds T and D and the reduction manifest. Returns the three paths."""
    control = load_section(control_path, section, SECTION_KEYS, path_keys=PATH_KEYS)
    prefix = control["prefix"]
    root_dir = _repository_root(Path(control["raw_bundle"]))
    control_path = Path(control_path)

    raw_path = Path(control["raw_bundle"])
    raw = cio.read(raw_path, "raw")
    try:
        table1 = raw["table1"].dataset
        pressure = np.asarray(table1["pressure_Pa"].values, dtype="float64")
        temperature = np.asarray(table1["temperature_K"].values, dtype="float64")
        height = np.asarray(table1["height_m"].values, dtype="float64")
        source = dict(raw["scalars/source"].attrs)
        latitude = dict(raw["scalars/latitude"].attrs)
        geodesy = dict(raw["scalars/geodesy"].attrs)
        surfaces = {
            "surface_100mbar": dict(raw["scalars/geodesy/surface_100mbar"].attrs),
            "surface_1bar": dict(raw["scalars/geodesy/surface_1bar"].attrs),
        }
        top_boundary = dict(raw["scalars/top_boundary"].attrs)
        rotation_scalars = dict(raw["scalars/rotation"].attrs)
        wind_scalars = dict(raw["scalars/wind_used_by_source"].attrs)
        gravity_scalars = dict(raw["scalars/gravity_used_by_source"].attrs)
    finally:
        raw.close()

    raw_hash = cio.input_hash_entry(raw_path, relative_to=root_dir)
    control_hash = cio.input_hash_entry(control_path, relative_to=root_dir)
    provenance_hashes = cio.input_hashes([raw_path, control_path], relative_to=root_dir)

    # ---- the four companions must already exist and validate -------------------------
    companion_hashes = {}
    for key, kind in COMPANIONS:
        path = Path(control[key])
        if not path.exists():
            raise ControlFileError(
                f"{path} does not exist. Step 9 completes the directory; it does not build "
                f"the kind {kind} file. Run its tool first."
            )
        handle = cio.read(path, kind)
        try:
            companion_hashes[key] = cio.input_hash_entry(path, relative_to=root_dir)
        finally:
            handle.close()

    # ---- kind T ----------------------------------------------------------------------
    nan = np.full(pressure.shape, np.nan)
    longitudes = np.asarray(latitude["longitude_system_iii_deg"], dtype="float64")
    thermo = xr.Dataset(
        {
            "pressure_Pa": (("level",), pressure,
                            _attrs("Pa", "pressure of the source profile", "derived",
                                   value_source="Table I", positive="down",
                                   direction="increasing")),
            "pressure_uncertainty_Pa": (("level",), nan.copy(),
                                        _attrs("Pa", "uncertainty on pressure", "derived",
                                               uncertainty_kind="stated",
                                               uncertainty_method="the source states none")),
            "temperature_K": (("level",), temperature,
                              _attrs("K", "temperature of the source profile", "derived",
                                     value_source="Table I")),
            "temperature_uncertainty_K": (("level",), nan.copy(),
                                          _attrs("K", "uncertainty on temperature", "derived",
                                                 uncertainty_kind="stated",
                                                 uncertainty_method="the source states none")),
            "height_m": (("level",), height,
                         _attrs("m", "altitude above the datum, the measured quantity",
                                "measured", value_source="Table I")),
            "height_uncertainty_m": (("level",), nan.copy(),
                                     _attrs("m", "uncertainty on height", "derived",
                                            uncertainty_kind="stated",
                                            uncertainty_method="the source states none")),
        }
    )
    thermo.attrs.update({
        "title": f"{prefix} thermodynamic profile, kind T",
        "profile_or_run": prefix,
        "role": "reduction",
        "source": source["citation"],
        "vertical_coordinate": "pressure_Pa",
        "thermo_instance": "source_profile",
        "latitude_planetocentric_absent_meaning": "point",
        "latitude_planetographic_deg": float(latitude["planetographic_deg"]),
        "latitude_planetographic_deg_value_source": latitude["value_source"],
        "latitude_planetographic_deg_swath": np.asarray(
            latitude["swath_planetographic_deg"], dtype="float64"),
        "latitude_planetographic_deg_swath_source": latitude["swath_source"],
        "latitude_planetographic_deg_uncertainty": float(latitude["uncertainty_deg"]),
        "latitude_planetographic_deg_uncertainty_kind": latitude["uncertainty_kind"],
        "latitude_definition": latitude["convention_source"],
        "spacecraft": source["spacecraft"],
        "event": source["event"],
        "observation_date": source["observation_date"],
        "frequency_bands": source["frequency_bands"],
        # The source states the longitude as a swath; the scalar is its first end, where the
        # ingress data begin, and the pair is kept beside it.
        "longitude_deg": float(longitudes[0]),
        "longitude_swath_deg": longitudes,
        "longitude_system": rotation_scalars["system"],
        "longitude_source": latitude["longitude_source"],
        "height_datum": (
            "the 1 bar pressure level as defined by the source; Table I gives the local "
            "altitude of the measurements relative to it"
        ),
        "source_top_boundary": top_boundary["statement"],
        # From the dedicated table the Step 9 review added to the transcription. It used to be
        # cut from the prose of `fit_inputs`, which is still carried below as the source's own
        # statement of everything that went into the geoid fit.
        "source_gravity_citation": gravity_scalars["citations"],
        "source_gravity_citation_value_source": gravity_scalars["value_source"],
        "source_pole_vector": gravity_scalars["pole_vector"],
        "source_fit_inputs": geodesy["fit_inputs"],
        "source_rotation_system": rotation_scalars["system"],
        "source_wind_citation": wind_scalars["citations"],
        "source_wind_vertical_structure": wind_scalars["vertical_structure"],
        "raw_bundle": raw_hash,
        "input_hashes": provenance_hashes,
        "control_file": control_hash,
    })
    cio.history_append(
        thermo,
        f"{TOOL}: pressure, temperature and height taken from the raw bundle table1 group; "
        "uncertainties present and NaN, the source states none; no composition and no geodesy "
        "in this file (SPEC_00 section 6.1)",
    )
    thermo_path = cio.write(Path(control["thermo_output"]), thermo, "thermo", created_by=TOOL)

    # ---- kind D ----------------------------------------------------------------------
    order = ["surface_100mbar", "surface_1bar"]
    def column(key, scale=1.0):
        return np.array([float(surfaces[s][key]) * scale for s in order], dtype="float64")

    residual_range = np.asarray(geodesy["fit_residual_km"], dtype="float64") * 1e3
    geodesy_ds = xr.Dataset(
        {
            "radius_equatorial_m": (("surface",), column("radius_equatorial_km", 1e3),
                                    _attrs("m", "fitted equatorial radius", "derived")),
            "radius_equatorial_uncertainty_m": (
                ("surface",), column("radius_equatorial_unc_km", 1e3),
                _attrs("m", "uncertainty on the equatorial radius", "derived",
                       uncertainty_kind="1sigma")),
            "radius_polar_m": (("surface",), column("radius_polar_mean_km", 1e3),
                               _attrs("m", "fitted mean polar radius, the anchor of the "
                                      "construction", "derived")),
            "radius_polar_uncertainty_m": (
                ("surface",), column("radius_polar_mean_unc_km", 1e3),
                _attrs("m", "uncertainty on the mean polar radius", "derived",
                       uncertainty_kind="1sigma")),
            "oblateness": (("surface",), column("oblateness"),
                           _attrs("1", "fitted oblateness", "derived")),
            "oblateness_uncertainty": (("surface",), column("oblateness_unc"),
                                       _attrs("1", "uncertainty on the oblateness", "derived",
                                              uncertainty_kind="1sigma")),
            # The source gives a range across its fits, not a value per surface, so the
            # per surface variable is NaN and the range is an attribute. SPEC_00 section 5:
            # NaN is the honest statement of no information at that location.
            "fit_residual_m": (("surface",), np.full(len(order), np.nan),
                               _attrs("m", "standard deviation of the fit", "derived",
                                      value_source="stated as a range, see "
                                                   "fit_residual_range_m")),
        },
        coords={
            "surface_pressure_Pa": (("surface",), column("pressure_Pa"),
                                    _attrs("Pa", "pressure of the fitted surface", "index",
                                           positive="down", direction="increasing")),
        },
    )
    geodesy_ds.attrs.update({
        "title": f"{prefix} fitted reference surfaces, kind D",
        "profile_or_run": prefix,
        "role": "reduction",
        "source": source["citation"],
        "citation": source["citation"],
        "fit_inputs": geodesy["fit_inputs"],
        "fit_latitude_convention": geodesy["fit_latitude_convention"],
        "fit_method": geodesy["fit_method"],
        "fit_residual_range_m": residual_range,
        "fit_residual_source": geodesy["fit_residual_source"],
        "reference_geoid_polar_radius_m": float(geodesy["reference_geoid_polar_radius_km"]) * 1e3,
        "reference_geoid_source": geodesy["reference_geoid_source"],
        "polar_asymmetry_note": surfaces["surface_100mbar"].get("polar_asymmetry_note", ""),
        "surface_value_sources": "\n".join(
            f"{s}: {surfaces[s]['value_source']}" for s in order),
        "raw_bundle": raw_hash,
        "input_hashes": provenance_hashes,
        "control_file": control_hash,
    })
    cio.history_append(
        geodesy_ds,
        f"{TOOL}: the two fitted surfaces taken from the raw bundle scalars/geodesy groups; "
        "the fit residual is stated by the source as a range across its fits, so the per "
        "surface variable is NaN and the range is in fit_residual_range_m",
    )
    geodesy_path = cio.write(Path(control["geodesy_output"]), geodesy_ds, "geodesy",
                             created_by=TOOL)

    # ---- the manifest, SPEC_00 section 7.1 -------------------------------------------
    manifest_path = Path(control["manifest"])
    directory = manifest_path.parent
    names = {
        "thermo": thermo_path.name,
        "composition": Path(control["composition_file"]).name,
        "geodesy": geodesy_path.name,
        "gravity": Path(control["gravity_file"]).name,
        "rotation": Path(control["rotation_file"]).name,
        "wind": Path(control["wind_file"]).name,
    }
    for kind_name, file_name in names.items():
        if not (directory / file_name).exists():
            raise ControlFileError(f"the manifest would name {file_name}, which does not exist")
        if not file_name.startswith(f"{prefix}_"):
            raise ControlFileError(
                f"{file_name} does not carry the directory prefix {prefix!r} (SPEC_00 "
                "section 2.2)"
            )

    lines = [
        f"# Reduction manifest for the {prefix} profile (SPEC_00 section 7.1).",
        f"# Written by {TOOL}. Pointers and declared choices only; every number the reduction",
        "# consumes is in the files named below.",
        "",
        "[profile]",
        f'slug        = "{prefix}"',
        f'description = "{control["description"]}"',
        "",
        "[inputs]",
    ]
    width = max(len(k) for k in names)
    for kind_name, file_name in names.items():
        lines.append(f'{kind_name:<{width}} = "{file_name}"')
    lines += [
        "",
        "[anchor_isobar]",
        f"pressure_Pa = {control['anchor_isobar_Pa']:.6g}",
        "",
        "[anchor_latitude]",
        f"fixed_point_tolerance_deg = {control['fixed_point_tolerance_deg']:.6g}",
        f"max_iterations            = {int(control['max_iterations'])}",
        "",
        "[geoid]",
        f"anchor_surface_Pa = {control['geoid_anchor_surface_Pa']:.6g}",
        f'anchor_quantity   = "{control["geoid_anchor_quantity"]}"',
        f'anchor_rule       = "{control["geoid_anchor_rule"]}"'
        "   # SPEC_01 v0.16: how the one constant of Eq. B3 is fixed",
        f"convergence_m     = {control['geoid_convergence_m']:.6g}",
        "",
        "[output]",
        f'product = "{Path(control["product"]).name}"',
        "",
    ]
    manifest_path.write_bytes("\n".join(lines).encode("utf-8"))
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
        print(f"wrote {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
