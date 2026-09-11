"""Write a kind R rotation system file from a `data_static` transcription.

SPEC_01 Step 3, SPEC_00 section 6.5. Tiny by design: its existence as a file lets a wind file
reference a rotation system by name and carry its value, and lets a reader compare the two.

Manuscript equations implemented: none, beyond the definition of angular rate from period,
`Omega = 2 pi / P`, which is recomputed here rather than transcribed. The static file carries
a rounded `angular_rate_rad_s` for a human reader; taking it at face value would lose the last
digits of the period, so the tool derives the rate and records that it did.
"""

from __future__ import annotations

import argparse
import math
import tomllib
from pathlib import Path

import xarray as xr

from casspian.lib import io as cio
from casspian.tools.gravity.control import (
    ControlFileError,
    load_section,
    reject_physical_values,
)

TOOL = "casspian-rotation-file"

SECTION_KEYS = {
    "source": True,
    "entry": True,
    "output": True,
    "prefix": True,
    "title": False,
}


def build(control_path, section: str = "rotation") -> Path:
    """Build a kind R file from the control file's named section."""
    control = load_section(control_path, section, SECTION_KEYS,
                           path_keys=("source", "output"))
    reject_physical_values(
        {k: v for k, v in control.items() if not isinstance(v, Path)}, control_path, section
    )
    source_path = Path(control["source"])
    output = Path(control["output"])
    prefix = control["prefix"]
    entry_name = control["entry"]

    if not source_path.exists():
        raise ControlFileError(f"{source_path}: the static transcription does not exist")
    with open(source_path, "rb") as handle:
        static = tomllib.load(handle)
    if entry_name not in static:
        raise ControlFileError(
            f"{source_path}: no [{entry_name}] entry; it holds {sorted(static)}"
        )
    entry = static[entry_name]

    period_s = float(entry["period_s"])
    angular_rate = 2.0 * math.pi / period_s
    transcribed_rate = entry.get("angular_rate_rad_s")

    def attrs(units, long_name, provenance, **extra):
        out = {"units": units, "long_name": long_name, "provenance": provenance}
        out.update({k: v for k, v in extra.items() if v not in (None, "")})
        return out

    dataset = xr.Dataset(
        {
            "period_s": (
                (),
                period_s,
                attrs("s", "rotation period of the system", "measured",
                      status=entry.get("status"), value_source=entry.get("value_source")),
            ),
            "angular_rate_rad_s": (
                (),
                angular_rate,
                attrs("rad s-1", "angular rate of the system", "derived",
                      value_source="computed as 2 pi / period_s by " + TOOL),
            ),
        },
        attrs={
            "title": control.get("title", f"{prefix} rotation system"),
            "profile_or_run": prefix,
            "role": "reduction",
            "source": entry.get("citation", source_path.name),
            "system_name": entry.get("name", entry_name),
            "epoch": entry.get("epoch", "unstated"),
            "citation": entry.get("citation", "unstated"),
            "input_hashes": cio.input_hashes([source_path, Path(control_path)]),
        },
    )
    if entry.get("note"):
        dataset.attrs["note"] = entry["note"]
    if transcribed_rate is not None:
        dataset.attrs["angular_rate_transcribed_rad_s"] = float(transcribed_rate)
        dataset.attrs["angular_rate_transcribed_note"] = (
            "the value written in the static file for a human reader; the variable holds "
            "2 pi / period_s at full precision and is the one to use"
        )

    cio.history_append(
        dataset,
        f"{TOOL}: period transcribed from {source_path.name} [{entry_name}]; "
        "angular_rate_rad_s recomputed as 2 pi / period_s",
    )
    return cio.write(output, dataset, "rotation", created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Write a kind R rotation system file from a control file."
    )
    parser.add_argument("control", help="path to the TOML control file")
    parser.add_argument("--section", default="rotation", help="section to read (default rotation)")
    args = parser.parse_args(argv)
    print(f"wrote {build(args.control, args.section)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
