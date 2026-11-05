"""What a profile's stage two shares: kind D from a raw bundle's geodesy, and the reduction manifest.

SPEC_18 v0.3 section 5, the general part of what was `tools/lindal/build_inputs.py`. The source's own
stage two reads its table and writes kind T; everything here operates on the standard form and is
called by every source's tool. A raw bundle's `scalars/geodesy` group and its two surface subgroups
are that form.

Manuscript equations implemented: none. No physical value appears in this module.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.control import (
    ANCHOR_QUANTITY_FOR_RULE,
    ControlFileError,
    load_section,
)
from casspian.lib.schema import WRITER_FILLED_GLOBALS, CasspianSchemaError

SECTION_KEYS = {
    # SPEC_01 v0.26: required, no default; written into kinds T and D.
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
    # SPEC_01 v0.20 Step 9 amendment: carried into the manifest's [diagnostics] section.
    "diagnostics_figures": False,
    "diagnostics_format": False,
    "diagnostics_dpi": False,
}

#: Attributes that record the writing of a file rather than its content. A kind T or D file
#: that differs from what this run would write only in these is the same file, verified and
#: kept (SPEC_01 v0.20 Step 9): the control file hash changes whenever a manifest choice does.
_WRITING_ATTRIBUTES = frozenset(WRITER_FILLED_GLOBALS) | {"history", "control_file",
                                                        "input_hashes"}

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


def _content(dataset, written):
    """`dataset` with the writing attributes dropped, and with only the globals in `written`.

    SPEC_R1: a global this tool no longer writes (`role`, `thermo_instance` and the source profile
    set) is not content, so a file written before SPEC_R1 that still carries one is kept.
    """
    kept = dataset.copy()
    kept.attrs = {k: v for k, v in dataset.attrs.items()
                  if k not in _WRITING_ATTRIBUTES and k in written}
    return kept


def write_unless_unchanged(path: Path, dataset, kind: str, tool: str) -> Path:
    """Write kind T or D, unless the file on disk already holds exactly this content.

    SPEC_01 v0.20 Step 9: a manifest choice changes nothing in kinds T and D, so they are
    verified and kept rather than rewritten, and their hashes and commits stand.
    """
    if path.exists():
        try:
            handle = cio.read(path, kind)
            try:
                existing = handle.load()
            finally:
                handle.close()
        except CasspianSchemaError:
            existing = None
        written = set(dataset.attrs)
        if existing is not None and _content(existing, written).identical(
                _content(dataset, written)):
            print(f"kept {path}: content unchanged, verified against what this run would write")
            return path
    return cio.write(path, dataset, kind, created_by=tool)


def attrs(units, long_name, provenance, **extra):
    out = {"units": units, "long_name": long_name, "provenance": provenance}
    for key, value in extra.items():
        if value is None or (isinstance(value, str) and not value):
            continue
        out[key] = value
    return out


def read_control(control_path, section: str) -> dict:
    """The stage two section of a build file, with its anchor and diagnostics choices checked."""
    control = load_section(control_path, section, SECTION_KEYS, path_keys=PATH_KEYS)
    expected_quantity = ANCHOR_QUANTITY_FOR_RULE.get(control["geoid_anchor_rule"])
    if expected_quantity is None or control["geoid_anchor_quantity"] != expected_quantity:
        raise ControlFileError(
            f"{control_path}: [{section}] geoid_anchor_rule = {control['geoid_anchor_rule']!r} "
            f"and geoid_anchor_quantity = {control['geoid_anchor_quantity']!r} do not belong "
            f"together; the manifest rules are {dict(ANCHOR_QUANTITY_FOR_RULE)} (SPEC_00 "
            "section 7.1 v0.15)."
        )
    if control.get("diagnostics_figures") and not all(
            k in control for k in ("diagnostics_format", "diagnostics_dpi")):
        raise ControlFileError(
            f"{control_path}: [{section}] diagnostics_figures = true needs diagnostics_format "
            "and diagnostics_dpi; no default is declared for them."
        )
    return control


def provenance(raw_path, control_path, product) -> dict:
    """SPEC_00 v0.18 section 5: every path recorded relative to the file that records it."""
    return {
        "raw_bundle": cio.input_hash_entry(raw_path, product),
        "input_hashes": cio.input_hashes([raw_path, control_path], product),
        "control_file": cio.input_hash_entry(control_path, product),
    }


def check_companions(control) -> None:
    """The four companions must already exist and validate."""
    for key, kind in COMPANIONS:
        path = Path(control[key])
        if not path.exists():
            raise ControlFileError(
                f"{path} does not exist. Step 9 completes the directory; it does not build "
                f"the kind {kind} file. Run its tool first."
            )
        handle = cio.read(path, kind)
        try:
            cio.input_hash_entry(path, Path(control["manifest"]))
        finally:
            handle.close()


def season_of(source) -> dict:
    """SPEC_01 v0.25 Step 9: kinds T and D carry the observation's date and season."""
    return {
        "epoch": str(source["observation_date"]),
        "solar_longitude_deg": float(source["solar_longitude_deg"]),
        "solar_longitude_source": str(source["solar_longitude_source"]),
    }


def write_geodesy(control, control_path, raw_path, tool: str) -> Path:
    """Kind D from the raw bundle's `scalars/geodesy` group.

    The group's own `citation`, where it has one, is the source of the surfaces: a profile may
    take its surfaces from another paper's fit (SPEC_18 section 1 item 7).
    """
    prefix = control["prefix"]
    geodesy_output = Path(control["geodesy_output"])
    raw = cio.read(raw_path, "raw")
    try:
        source = dict(raw["scalars/source"].attrs)
        geodesy = dict(raw["scalars/geodesy"].attrs)
        surfaces = {
            "surface_100mbar": dict(raw["scalars/geodesy/surface_100mbar"].attrs),
            "surface_1bar": dict(raw["scalars/geodesy/surface_1bar"].attrs),
        }
    finally:
        raw.close()

    order = ["surface_100mbar", "surface_1bar"]

    def column(key, scale=1.0):
        return np.array([float(surfaces[s][key]) * scale for s in order], dtype="float64")

    residual_range = np.asarray(geodesy["fit_residual_km"], dtype="float64") * 1e3
    citation = geodesy.get("citation", source["citation"])
    geodesy_ds = xr.Dataset(
        {
            "radius_equatorial_m": (("surface",), column("radius_equatorial_km", 1e3),
                                    attrs("m", "fitted equatorial radius", "derived")),
            "radius_equatorial_uncertainty_m": (
                ("surface",), column("radius_equatorial_unc_km", 1e3),
                attrs("m", "uncertainty on the equatorial radius", "derived",
                      uncertainty_kind="1sigma")),
            "radius_polar_m": (("surface",), column("radius_polar_mean_km", 1e3),
                               attrs("m", "fitted mean polar radius, the anchor of the "
                                     "construction", "derived")),
            "radius_polar_uncertainty_m": (
                ("surface",), column("radius_polar_mean_unc_km", 1e3),
                attrs("m", "uncertainty on the mean polar radius", "derived",
                      uncertainty_kind="1sigma")),
            "oblateness": (("surface",), column("oblateness"),
                           attrs("1", "fitted oblateness", "derived")),
            "oblateness_uncertainty": (("surface",), column("oblateness_unc"),
                                       attrs("1", "uncertainty on the oblateness", "derived",
                                             uncertainty_kind="1sigma")),
            # The source gives a range across its fits, not a value per surface, so the
            # per surface variable is NaN and the range is an attribute. SPEC_00 section 5:
            # NaN is the honest statement of no information at that location.
            "fit_residual_m": (("surface",), np.full(len(order), np.nan),
                               attrs("m", "standard deviation of the fit", "derived",
                                     value_source="stated as a range, see "
                                                  "fit_residual_range_m")),
        },
        coords={
            "surface_pressure_Pa": (("surface",), column("pressure_Pa"),
                                    attrs("Pa", "pressure of the fitted surface", "index",
                                          positive="down", direction="increasing")),
        },
    )
    geodesy_ds.attrs.update({
        "title": f"{prefix} fitted reference surfaces, kind D",
        **season_of(source),
        "source": citation,
        "citation": citation,
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
        **provenance(raw_path, control_path, geodesy_output),
    })
    cio.history_append(
        geodesy_ds,
        f"{tool}: the two fitted surfaces taken from the raw bundle scalars/geodesy groups; "
        "the fit residual is stated by the source as a range across its fits, so the per "
        "surface variable is NaN and the range is in fit_residual_range_m",
    )
    return write_unless_unchanged(geodesy_output, geodesy_ds, "geodesy", tool)


def write_manifest(control, thermo_path: Path, geodesy_path: Path, tool: str) -> Path:
    """The reduction manifest, SPEC_00 section 7.1."""
    prefix = control["prefix"]
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

    lines = [
        f"# Reduction manifest for the {prefix} profile (SPEC_00 section 7.1).",
        f"# Written by {tool}. Pointers and declared choices only; every number the reduction",
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
    if "diagnostics_figures" in control:
        lines += [
            "[diagnostics]",
            f"figures = {'true' if control['diagnostics_figures'] else 'false'}",
        ]
        if "diagnostics_format" in control:
            lines.append(f'format  = "{control["diagnostics_format"]}"')
        if "diagnostics_dpi" in control:
            lines.append(f"dpi     = {int(control['diagnostics_dpi'])}")
        lines.append("")
    manifest_path.write_bytes("\n".join(lines).encode("utf-8"))
    return manifest_path
