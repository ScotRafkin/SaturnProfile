"""Write a kind G harmonic set file from a `data_static` transcription.

SPEC_01 Step 3, SPEC_00 section 6.4. Driven by a TOML control file that names the static
source, the output path, and the prefix. No physical value appears in this module: every
number comes from the transcription the control file points at.

Manuscript equations implemented: none. The potential the coefficients belong to is
`lib.gravity` (Step 4); this tool records which convention they obey, as the code
`CASSPIAN-J1` defined in `lib.schema`.

Two things the source files carry that the netCDF keeps separate. The `fitted` flag per degree
becomes the `J_status` flag variable (0 fitted, 1 assumed, 2 constrained), which says how the
number was obtained. The `status` string per entry is the transcription status, whether a
second reader has checked the value against the paper, and it becomes a per-variable attribute.
They answer different questions and are not merged.
"""

from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.schema import HARMONIC_CONVENTION
from casspian.tools.gravity.control import (
    ControlFileError,
    load_section,
    reject_physical_values,
)

TOOL = "casspian-gravity-file"

SECTION_KEYS = {"source": True, "output": True, "prefix": True, "title": False}

#: SPEC_00 section 6.4. CF style flag pair for `J_status`.
J_STATUS_VALUES = np.array([0, 1, 2], dtype="int8")
J_STATUS_MEANINGS = "fitted assumed constrained"


def _collect(entries: dict, key: str) -> str:
    """Join a per degree string attribute, collapsing to one value when they all agree."""
    values = {degree: table.get(key) for degree, table in entries.items()}
    present = {d: v for d, v in values.items() if v}
    if not present:
        return ""
    unique = set(present.values())
    if len(unique) == 1:
        return unique.pop()
    return "\n".join(f"degree {d}: {present[d]}" for d in sorted(present, key=int))


def build(control_path, section: str = "gravity") -> Path:
    """Build a kind G file from the control file's named section."""
    control = load_section(control_path, section, SECTION_KEYS,
                           path_keys=("source", "output"))
    reject_physical_values(
        {k: v for k, v in control.items() if not isinstance(v, Path)}, control_path, section
    )
    source_path = Path(control["source"])
    output = Path(control["output"])
    prefix = control["prefix"]

    if not source_path.exists():
        raise ControlFileError(f"{source_path}: the static transcription does not exist")
    with open(source_path, "rb") as handle:
        static = tomllib.load(handle)

    meta = static.get("meta", {})
    entries = static.get("J", {})
    if not entries:
        raise ControlFileError(f"{source_path}: no [J.<degree>] tables")
    gm = static.get("GM", {})

    degrees = sorted(entries, key=int)
    degree_values = np.array([int(d) for d in degrees], dtype="int32")
    j_values = np.array([float(entries[d]["value"]) for d in degrees], dtype="float64")
    j_uncertainty = np.array(
        [float(entries[d].get("uncertainty", np.nan)) for d in degrees], dtype="float64"
    )
    j_status = np.array(
        [0 if entries[d].get("fitted", False) else 1 for d in degrees], dtype="int8"
    )

    def _absent(value) -> bool:
        # An array attribute such as flag_values is never "empty" in the scalar sense, and
        # testing it with `in (None, "")` raises rather than answering.
        if isinstance(value, np.ndarray):
            return value.size == 0
        return value is None or value == ""

    def attrs(units, long_name, provenance, **extra):
        out = {"units": units, "long_name": long_name, "provenance": provenance}
        out.update({k: v for k, v in extra.items() if not _absent(v)})
        return out

    dataset = xr.Dataset(
        {
            "J": (
                ("degree",),
                j_values,
                # SPEC_00 section 5: fitted to tracking data with no added physical
                # assumption, which is `derived` rather than `inferred`.
                attrs("1", "zonal harmonic coefficient, unscaled", "derived",
                      status=_collect(entries, "status"),
                      value_source=_collect(entries, "value_source")),
            ),
            "J_uncertainty": (
                ("degree",),
                j_uncertainty,
                attrs("1", "uncertainty on the zonal harmonic coefficient", "derived",
                      uncertainty_kind="1sigma",
                      uncertainty_method=meta.get("uncertainty_scaling",
                                                  "as published by the source")),
            ),
            "J_status": (
                ("degree",),
                j_status,
                attrs("1", "how the coefficient was obtained", "assumed",
                      flag_values=J_STATUS_VALUES, flag_meanings=J_STATUS_MEANINGS),
            ),
            "GM_m3s2": (
                (),
                float(gm["planet_value_m3s2"]),
                attrs("m3 s-2", "planet gravitational parameter", "derived",
                      status=gm.get("planet_status"), value_source=gm.get("planet_source")),
            ),
            "GM_uncertainty_m3s2": (
                (),
                float(gm.get("planet_uncertainty_m3s2", np.nan)),
                attrs("m3 s-2", "uncertainty on the planet gravitational parameter",
                      "derived", uncertainty_kind="1sigma"),
            ),
            "normalization_radius_m": (
                (),
                float(meta["normalization_radius_m"]),
                attrs("m", "radius the harmonic coefficients are normalized to", "assumed",
                      status=meta.get("normalization_status"),
                      value_source=meta.get("normalization_source")),
            ),
        },
        coords={
            "degree": (
                ("degree",),
                degree_values,
                # A label, not a measurement (SPEC_00 section 5, `index`).
                attrs("1", "degree of the zonal harmonic", "index"),
            )
        },
    )

    globals_ = {
        "title": control.get("title", f"{prefix} gravity harmonic set"),
        "profile_or_run": prefix,
        "role": "reduction",
        "source": meta.get("citation", str(source_path.name)),
        "GM_scope": gm.get("scope", "unstated"),
        "epoch": meta.get("epoch", "unstated"),
        "harmonic_convention": HARMONIC_CONVENTION,
        "input_hashes": cio.input_hashes([source_path, Path(control_path)]),
    }
    # The source file states its convention in prose. The netCDF asserts the code; the prose is
    # kept beside it so a reader can see what the transcription claimed.
    if meta.get("harmonic_convention"):
        globals_["harmonic_convention_source_statement"] = meta["harmonic_convention"]
    if meta.get("uncertainty_scaling"):
        globals_["uncertainty_scaling"] = meta["uncertainty_scaling"]
    # SPEC_01 Step 3: the planet GM goes into the file, the system GM into an attribute.
    if gm.get("system_value_m3s2") is not None:
        globals_["GM_system_m3s2"] = float(gm["system_value_m3s2"])
        globals_["GM_system_uncertainty_m3s2"] = float(gm.get("system_uncertainty_m3s2", np.nan))
        globals_["GM_system_source"] = gm.get("system_source", "")
    if gm.get("note"):
        globals_["GM_derivation_note"] = gm["note"]
    if meta.get("rotation_note"):
        globals_["rotation_note"] = meta["rotation_note"]
    dataset.attrs.update({k: v for k, v in globals_.items() if v not in (None, "")})

    cio.history_append(
        dataset,
        f"{TOOL}: transcribed from {source_path.name}; harmonic_convention asserted as "
        f"{HARMONIC_CONVENTION}; J values unscaled as transcribed",
    )
    return cio.write(output, dataset, "gravity", created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Write a kind G harmonic set file from a control file."
    )
    parser.add_argument("control", help="path to the TOML control file")
    parser.add_argument("--section", default="gravity", help="section to read (default gravity)")
    args = parser.parse_args(argv)
    print(f"wrote {build(args.control, args.section)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
