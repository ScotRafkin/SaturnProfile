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

**The pressure grid (SPEC_01 v0.23 Step 2, SPEC_03 v0.5 Step 0).** The transcription's
`[pressure_grid]` table declares that the printed pressures are rounded representations of the
grid `10^(k/denominator)` in the declared unit. `table1` then carries the printed value as
`pressure_printed_Pa` and the grid value as `pressure_Pa`. For each row the candidates are every
integer `k` whose grid value, rounded to the decimals printed in that row, equals the printed
value. One candidate fixes the row. Rows with several are resolved one at a time from the fixed
side: each takes the candidate that continues the spacing in `k` of the two nearest fixed rows
on that side, both sides agreeing when both are fixed, and a resolved row counts as fixed for
the next. A row still ambiguous, or with no candidate and not in the exclusion list, is refused.
No spacing is written here; the spacing is a property of the table, which the acceptance checks.
"""

from __future__ import annotations

import argparse
import csv
import math
import tomllib
from decimal import ROUND_HALF_EVEN, Decimal
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


#: The unit Table I prints its pressures in, from the CSV column name `pressure_mbar`. The
#: declared grid must be stated in the same unit, since the rounding test is on the printed text.
PRINTED_PRESSURE_UNIT = "mbar"


class PressureGridError(ValueError):
    """A printed pressure cannot be placed on the declared grid under the declared rule."""


def grid_candidates(text: str, denominator: int) -> list[int]:
    """Every integer `k` whose grid value `10^(k/denominator)`, rounded to the decimals printed
    in `text`, equals the printed value. The grid value is rounded from its exact binary value."""
    printed = Decimal(text.strip())
    quantum = Decimal(1).scaleb(printed.as_tuple().exponent)
    low, high = printed - quantum / 2, printed + quantum / 2
    if low <= 0:
        return []
    # The bracket is widened by one on each side so that float logarithms cannot drop an end.
    k_low = math.floor(denominator * math.log10(float(low))) - 1
    k_high = math.ceil(denominator * math.log10(float(high))) + 1
    return [k for k in range(k_low, k_high + 1)
            if Decimal(10.0 ** (k / denominator)).quantize(quantum, rounding=ROUND_HALF_EVEN)
            == printed]


def resolve_grid(texts, candidates, excluded, names):
    """The grid index of every row under the declared rule, and the order ambiguous rows took.

    `texts`, `candidates` and `excluded` are in increasing pressure. `names` describes each row
    for a refusal. Returns `(k, resolved)`: `k[i]` is None for an excluded row, and `resolved`
    lists `(row, k, side)` in the order the ambiguous rows were settled.
    """
    size = len(texts)
    k = [None] * size
    ambiguous = []
    for i in range(size):
        if excluded[i]:
            continue
        if not candidates[i]:
            raise PressureGridError(
                f"{names[i]}: no value of the declared grid rounds to the printed "
                f"{texts[i]} {PRINTED_PRESSURE_UNIT}, and the row is not in "
                "excluded_printed_values_mbar. Refused rather than kept as printed."
            )
        if len(candidates[i]) == 1:
            k[i] = candidates[i][0]
        else:
            ambiguous.append(i)

    resolved = []
    while ambiguous:
        for i in ambiguous:
            predictions = {}
            for side, step in (("below", 1), ("above", -1)):
                near, far = i + step, i + 2 * step
                if 0 <= far < size and k[near] is not None and k[far] is not None:
                    predictions[side] = 2 * k[near] - k[far]
            if not predictions:
                continue
            values = set(predictions.values())
            if len(values) != 1:
                raise PressureGridError(
                    f"{names[i]}: the two fixed sides continue their spacing to different grid "
                    f"indices {predictions}; candidates {candidates[i]}. Refused."
                )
            value = values.pop()
            if value not in candidates[i]:
                raise PressureGridError(
                    f"{names[i]}: the spacing of the fixed rows {sorted(predictions)} continues "
                    f"to k = {value}, which does not round to the printed {texts[i]} "
                    f"{PRINTED_PRESSURE_UNIT}; candidates {candidates[i]}. Refused."
                )
            k[i] = value
            resolved.append((i, value, " and ".join(sorted(predictions))))
            ambiguous.remove(i)
            break
        else:
            rows = "; ".join(f"{names[i]} candidates {candidates[i]}" for i in ambiguous)
            raise PressureGridError(
                f"rows still ambiguous with no fixed side to resolve them from: {rows}. Refused."
            )
    return k, resolved


def read_table1(path: Path, grid: dict) -> tuple[xr.Dataset, list[str]]:
    """Read the Table I transcription into the `table1` group, ordered by increasing pressure.

    Returns the group and the lines for `history` that record how the grid was applied.
    """
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

    # ---- the declared pressure grid -------------------------------------------------------
    denominator = int(grid["denominator"])
    if grid["unit"] != PRINTED_PRESSURE_UNIT or grid["rule"] != (
            f"10^(k/{denominator}) {PRINTED_PRESSURE_UNIT}"):
        raise PressureGridError(
            f"[pressure_grid] declares rule {grid['rule']!r}, denominator {denominator} and "
            f"unit {grid['unit']!r}; the rule must read '10^(k/{denominator}) "
            f"{PRINTED_PRESSURE_UNIT}', in the unit Table I prints."
        )
    excluded_values = {Decimal(str(v)) for v in grid["excluded_printed_values_mbar"]}
    texts = [rows[i]["pressure_mbar"].strip() for i in order]
    names = [f"row {int(i) + 1} of {path.name} ({rows[i]['pressure_mbar'].strip()} "
             f"{PRINTED_PRESSURE_UNIT})" for i in order]
    excluded = [Decimal(t) in excluded_values for t in texts]
    candidates = [[] if ex else grid_candidates(t, denominator)
                  for t, ex in zip(texts, excluded)]
    k, resolved = resolve_grid(texts, candidates, excluded, names)

    printed = dataset["pressure_Pa"].values
    scale = Decimal(COLUMNS[0][2])
    on_grid = np.array([printed[i] if k[i] is None
                        else float(Decimal(10.0 ** (k[i] / denominator)) * scale)
                        for i in range(len(texts))], dtype="float64")
    applied = np.array([0 if value is None else 1 for value in k], dtype="int8")

    dataset = dataset.rename({"pressure_Pa": "pressure_printed_Pa"})
    dataset["pressure_printed_Pa"].attrs["long_name"] = "pressure as printed"
    dataset["pressure_Pa"] = (("level",), on_grid, {
        "units": "Pa",
        "long_name": "pressure on the grid declared in scalars/pressure_grid",
        "provenance": "derived",
        "value_source": "Table I, placed on the grid declared in scalars/pressure_grid; the "
                        "printed value is pressure_printed_Pa",
        "positive": "down",
        "direction": "increasing",
    })
    dataset["pressure_grid_index"] = (("level",), np.array(
        [np.nan if value is None else float(value) for value in k], dtype="float64"), {
        "units": "1",
        "long_name": f"the integer k of the grid value 10^(k/{denominator}) "
                     f"{PRINTED_PRESSURE_UNIT}; NaN where the row keeps its printed value",
        "provenance": "index",
    })
    dataset["pressure_grid_applied"] = (("level",), applied, {
        "units": "1",
        "long_name": "whether pressure_Pa is the grid value or the printed value",
        "provenance": "index",
        "flag_values": np.array([0, 1], dtype="int8"),
        "flag_meanings": "as_printed snapped_to_grid",
    })
    dataset = dataset[["pressure_printed_Pa", "pressure_Pa", "pressure_grid_index",
                       "pressure_grid_applied", "temperature_K", "nh3_mole_fraction",
                       "height_m"]]
    dataset.attrs["source_columns"] = ", ".join(name for name, *_ in COLUMNS)
    dataset.attrs["level_order"] = "increasing pressure"
    dataset.attrs["rows_reordered_from_file"] = "yes" if reordered else "no"

    single = sum(1 for c, ex in zip(candidates, excluded) if not ex and len(c) == 1)
    lines = [
        f"pressure grid applied under [pressure_grid] rule {grid['rule']!r}: for each row the "
        f"integers k whose grid value, rounded to the decimals printed in the row, equals the "
        f"printed value; one candidate fixes the row; ambiguous rows resolved one at a time "
        f"from the fixed side by continuing the spacing in k of the two nearest fixed rows, "
        f"both sides agreeing when both are fixed; excluded rows keep their printed value",
        f"pressure grid: {single} rows with one candidate, {len(resolved)} resolved by spacing, "
        f"{sum(excluded)} excluded and kept as printed, of {len(texts)}",
    ]
    for i, value, side in resolved:
        lines.append(f"pressure grid: {names[i]} had candidates {candidates[i]}; took k = "
                     f"{value} from the fixed rows {side}")
    lines.append(
        "pressure grid note (SPEC_03 v0.5 Step 0): with k = 0 at 1 mbar the top row is k = -70 "
        "(19.9526 Pa). A nearest integer snap lands 0.32 mbar on k = -49 (0.3236 mbar) instead "
        "of k = -50 (0.3162 mbar), and both round to 0.32; the equal spacing rule settles it, "
        "and the acceptance checks the resulting k sequence"
    )
    return dataset, lines


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


#: SPEC_01 v0.25 Step 2: the keys of `[source]` that date the observation and name its season.
SEASON_KEYS = ("observation_date", "solar_longitude_deg", "subsolar_latitude_deg",
               "solar_longitude_source")


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

    with open(scalars_toml, "rb") as handle:
        tables = tomllib.load(handle)
    citation = tables.get("source", {}).get("citation", "unknown")
    source_table = tables.get("source", {})
    missing = [key for key in SEASON_KEYS if key not in source_table]
    if missing:
        raise ValueError(
            f"{scalars_toml}: [source] lacks {missing}; the raw bundle carries the observation's "
            "date as its epoch and its season with the source of the computation (SPEC_00 v0.17 "
            "section 5, SPEC_01 v0.25 Step 2)."
        )
    if "pressure_grid" not in tables:
        raise PressureGridError(
            f"{scalars_toml}: no [pressure_grid] table. SPEC_01 v0.23 Step 2 places Table I on a "
            "declared grid, and the declaration is data in the transcription, not a default here."
        )

    table1, grid_lines = read_table1(table1_csv, tables["pressure_grid"])
    groups = {"table1": table1}
    groups.update(read_scalars(scalars_toml))

    data_inputs = [table1_csv, scalars_toml]
    record = [*data_inputs] + ([notes_md] if notes_md.exists() else [])

    root = xr.Dataset(
        attrs={
            "title": "Lindal et al. (1985) Voyager 2 ingress, raw transcription bundle",
            "profile_or_run": "lindal",
            # The raw bundle is written from the transcription with no build-file section, so
            # its role is the transcription's own: the reduction of a source.
            "role": "reduction",
            "source": citation,
            "epoch": str(source_table["observation_date"]),
            "solar_longitude_deg": float(source_table["solar_longitude_deg"]),
            "subsolar_latitude_deg": float(source_table["subsolar_latitude_deg"]),
            "solar_longitude_source": str(source_table["solar_longitude_source"]),
            "input_hashes": cio.input_hashes(data_inputs, output),
            "raw_sources": cio.input_hashes(record, output),
        }
    )
    for line in CONVERSIONS:
        cio.history_append(root, f"{TOOL}: {line}")
    for line in grid_lines:
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
