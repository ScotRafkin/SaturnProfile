"""Stage one of the Lindal tool: the two ASCII transcriptions into the raw bundle.

SPEC_01 Step 2, SPEC_00 section 2.2.1. Reads `lindal_table1.csv` and `lindal_scalars.toml`
from a profile's `raw/` directory and writes `lindal_raw.nc` beside them. The bundle is the
transcription in netCDF form: it is the file a reviewer opens to ask what Lindal actually said.

Manuscript equations implemented: none. This tool transcribes and converts units. It
interprets nothing, fills nothing, and drops nothing.

Unit conversion. Table I is printed in mbar, ppm and km, and SPEC_00 section 5 requires SI, so
the three columns are converted once, here, and the conversions are named in `history`. The
conversion is done in exact decimal arithmetic on the text as written in the CSV, then rounded
once to float64. Multiplying the parsed double instead disagrees in 15 of the 207 transcribed
cells: 2.51 mbar becomes 250.99999999999997 Pa rather than 251.0, and 5.1 ppm becomes
5.0999999999999995e-06 rather than 5.1e-06. None of those cells is one the Step 2 acceptance
names, so the acceptance alone would not have caught it.

The scalars are not converted. Every table of the TOML becomes a sub-group of `scalars` with
that table's keys as attributes, numeric where the value is numeric and string otherwise, so
that each `value_source` sits beside the value it justifies.
"""

from __future__ import annotations

import argparse
import csv
import tomllib
from decimal import Decimal
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import io as cio

TOOL = "casspian-lindal-raw"

#: Column name in the CSV, variable name in the bundle, decimal scale factor, units,
#: long name, and provenance. SPEC_01 Step 2: `measured` for height, `derived` for the rest.
COLUMNS = (
    ("pressure_mbar", "pressure_Pa", "100", "Pa", "pressure", "derived"),
    ("temperature_K", "temperature_K", "1", "K", "temperature", "derived"),
    ("nh3_ppm", "nh3_mole_fraction", "1e-6", "1", "ammonia mole fraction", "derived"),
    ("altitude_km", "height_m", "1000", "m", "altitude above the source 1 bar level", "measured"),
)

CONVERSIONS = (
    "pressure converted from mbar to Pa by a factor of 100",
    "ammonia converted from ppm to mole fraction by a factor of 1e-6",
    "altitude converted from km to m by a factor of 1000",
)


def _scaled(text: str, scale: str) -> float:
    """Convert one printed value to SI, exactly in decimal, then once to float64.

    A blank cell is no information at that location and becomes NaN. SPEC_00 section 5: the
    reader never fills a NaN, and filling is the job of a tool under a declared rule.
    """
    text = text.strip()
    if not text:
        return float("nan")
    return float(Decimal(text) * Decimal(scale))


def read_table1(path: Path) -> xr.Dataset:
    """Read the Table I transcription into the `table1` group, ordered by increasing pressure."""
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path}: no data rows")

    columns = {}
    for source_name, target_name, scale, units, long_name, provenance in COLUMNS:
        if source_name not in rows[0]:
            raise ValueError(f"{path}: expected a column named {source_name!r}")
        values = np.array([_scaled(row[source_name], scale) for row in rows], dtype="float64")
        columns[target_name] = (values, units, long_name, provenance)

    order = np.argsort(columns["pressure_Pa"][0], kind="stable")
    reordered = not np.array_equal(order, np.arange(order.size))

    dataset = xr.Dataset(
        {
            name: (
                ("level",),
                values[order],
                {
                    "units": units,
                    "long_name": long_name,
                    "provenance": provenance,
                    "value_source": "Table I",
                },
            )
            for name, (values, units, long_name, provenance) in columns.items()
        }
    )
    dataset["pressure_Pa"].attrs.update({"positive": "down", "direction": "increasing"})
    dataset.attrs["source_columns"] = ", ".join(name for name, *_ in COLUMNS)
    dataset.attrs["level_order"] = "increasing pressure"
    dataset.attrs["rows_reordered_from_file"] = "yes" if reordered else "no"
    return dataset


def _attribute_value(value):
    """Render one TOML value as a netCDF attribute, numeric where it is numeric."""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, list):
        if value and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value):
            return np.array(value, dtype="float64")
        return "\n".join(str(v) for v in value)
    return str(value)


def read_scalars(path: Path) -> dict[str, xr.Dataset]:
    """Read the scalars transcription into one group per TOML table, nested tables included."""
    with open(path, "rb") as handle:
        tables = tomllib.load(handle)

    groups: dict[str, xr.Dataset] = {"scalars": xr.Dataset(attrs={
        "content": "every scalar the source states, transcribed; one sub-group per TOML table",
        "source_file": path.name,
    })}

    def add(prefix: str, table: dict) -> None:
        flat = {k: _attribute_value(v) for k, v in table.items() if not isinstance(v, dict)}
        groups[prefix] = xr.Dataset(attrs=flat)
        for key, value in table.items():
            if isinstance(value, dict):
                add(f"{prefix}/{key}", value)

    for name, table in tables.items():
        add(f"scalars/{name}", table)
    return groups


def _repository_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        if (candidate / ".git").exists():
            return candidate
    return start


def build(raw_dir: Path, output: Path) -> Path:
    """Build the raw bundle from the two transcriptions in `raw_dir`."""
    raw_dir = Path(raw_dir).resolve()
    output = Path(output).resolve()
    table1_csv = raw_dir / "lindal_table1.csv"
    scalars_toml = raw_dir / "lindal_scalars.toml"
    notes_md = raw_dir / "notes.md"

    for required in (table1_csv, scalars_toml):
        if not required.exists():
            raise FileNotFoundError(f"{required} is required by {TOOL} and does not exist")

    table1 = read_table1(table1_csv)
    groups = {"table1": table1}
    groups.update(read_scalars(scalars_toml))

    with open(scalars_toml, "rb") as handle:
        tables = tomllib.load(handle)
    citation = tables.get("source", {}).get("citation", "unknown")

    root_dir = _repository_root(raw_dir)
    data_inputs = [table1_csv, scalars_toml]
    record = [*data_inputs] + ([notes_md] if notes_md.exists() else [])

    root = xr.Dataset(
        attrs={
            "title": "Lindal et al. (1985) Voyager 2 ingress, raw transcription bundle",
            "profile_or_run": "lindal",
            "role": "reduction",
            "source": citation,
            "input_hashes": cio.input_hashes(data_inputs, relative_to=root_dir),
            "raw_sources": cio.input_hashes(record, relative_to=root_dir),
        }
    )
    for line in CONVERSIONS:
        cio.history_append(root, f"{TOOL}: {line}")
    if notes_md.exists():
        cio.history_append(
            root, f"{TOOL}: notes.md hashed into raw_sources as documentation, not read"
        )

    return cio.write(output, root, "raw", groups=groups, created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL,
        description="Build the Lindal raw bundle from the two ASCII transcriptions.",
    )
    parser.add_argument(
        "--raw-dir",
        default="occul_data/lindal/raw",
        help="directory holding lindal_table1.csv and lindal_scalars.toml",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="output path; defaults to <raw-dir>/lindal_raw.nc",
    )
    args = parser.parse_args(argv)
    raw_dir = Path(args.raw_dir)
    output = Path(args.output) if args.output else raw_dir / "lindal_raw.nc"
    written = build(raw_dir, output)
    print(f"wrote {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
