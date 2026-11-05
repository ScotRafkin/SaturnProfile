"""Stage one of the Kliore tool: the transcription and the scalars into the raw bundle.

SPEC_18 v0.3. Reads `occul_data/kliore1980/pioneer11_exit.csv` (Kliore et al. 1980, Table 1, 128
levels, untouched) and `raw/kliore1980_scalars.toml`, and writes `raw/kliore1980_raw.nc`.

Unlike Lindal's bundle this one is not the bare transcription, because Table 1 is not yet a profile
the reduction can take. Two things happen here and nowhere else:

1. **The equatorial mapping is undone.** Table 1's radius is an equatorial radius, mapped from the
   measured one with "a constant oblateness value based on an ellipsoidal approximation" (p. 5861).
   `r_local = r_eq (1 - f) / sqrt((1 - f)^2 cos^2 phi + sin^2 phi)`, with f = 0.096 (p. 5857) and
   phi planetocentric, the reading SPEC_16 ruled.
2. **One level is inserted at exactly 100 mbar (SPEC_18 section 2).** The reduction needs its anchor
   isobar to be a tabulated level, and Table 1 has 99.6 and 101.9 mbar on either side. Temperature and
   radius are linear in ln p between those two. It is the one row of `table1` with no tabulated
   refractivity, and it is identified by its pressure, 100 mbar exactly (SPEC_18 section 5 ruling 5).
   Kind C reads its levels from this table, so the level is here and not only in kind T.

Kliore's uncertainty ranges stay here and go no further (SPEC_18 section 5 ruling 3).

Manuscript equations implemented: none. No physical value appears in this module.
"""

from __future__ import annotations

import argparse
import csv
import math
import tomllib
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.tools.lindal.build_raw import SEASON_KEYS, _scaled, read_scalars

TOOL = "casspian-kliore-raw"

#: The anchor isobar in Pa, the level SPEC_18 section 2 inserts.
INSERTED_PRESSURE_Pa = 1.0e4

#: Column in the CSV, variable in the bundle, decimal scale, units, long name.
COLUMNS = (
    ("pressure_mbar", "pressure_Pa", "100", "Pa", "pressure"),
    ("temperature_K", "temperature_K", "1", "K", "temperature"),
    ("refractivity_n_units", "refractivity_n_units", "1", "1", "refractivity times 1e6, as tabulated"),
    ("equatorial_radius_km", "radius_equatorial_m", "1000", "m", "radius mapped to the equator, as tabulated"),
    ("temperature_uncertainty_minus_K", "temperature_range_minus_K", "1", "K",
     "lower extreme of temperature, Kliore's range, tabulated on every fifth row"),
    ("temperature_uncertainty_plus_K", "temperature_range_plus_K", "1", "K",
     "upper extreme of temperature, Kliore's range, tabulated on every fifth row"),
    ("pressure_uncertainty_minus_mbar", "pressure_range_minus_Pa", "100", "Pa",
     "lower extreme of pressure, Kliore's range, tabulated on every fifth row"),
    ("pressure_uncertainty_plus_mbar", "pressure_range_plus_Pa", "100", "Pa",
     "upper extreme of pressure, Kliore's range, tabulated on every fifth row"),
)


def read_table(path: Path) -> dict:
    """The transcription's columns in SI, in file order (increasing pressure), each cell converted
    once in exact decimal. A row's trailing `#` comment records a corrected misprint."""
    with open(path, encoding="utf-8", newline="") as handle:
        lines = [line.split("#")[0] for line in handle if not line.lstrip().startswith("#")]
    rows = list(csv.DictReader(lines))
    return {target: np.array([_scaled(row[source] or "", scale) for row in rows], dtype="float64")
            for source, target, scale, _, _ in COLUMNS}


def insert_level(table: dict, pressure_Pa: float) -> dict:
    """The table with one row at `pressure_Pa`: temperature and radius linear in ln p between the
    neighbors, every other column NaN there."""
    p = table["pressure_Pa"]
    if np.any(p == pressure_Pa):
        raise ValueError(f"{pressure_Pa} Pa is already a tabulated level")
    k = int(np.searchsorted(p, pressure_Pa))
    weight = math.log(pressure_Pa / p[k - 1]) / math.log(p[k] / p[k - 1])
    out = {}
    for name, column in table.items():
        if name == "pressure_Pa":
            value = pressure_Pa
        elif name in ("temperature_K", "radius_equatorial_m"):
            value = column[k - 1] + weight * (column[k] - column[k - 1])
        else:
            value = np.nan
        out[name] = np.insert(column, k, value)
    return out


def local_radius(radius_equatorial_m, latitude_deg: float, oblateness: float):
    """Kliore's mapping undone, at his planetocentric latitude."""
    phi = math.radians(latitude_deg)
    f = oblateness
    return np.asarray(radius_equatorial_m) * (1.0 - f) / math.sqrt(
        (1.0 - f) ** 2 * math.cos(phi) ** 2 + math.sin(phi) ** 2)


def build(raw_dir, output=None, table_csv=None) -> Path:
    """Build the raw bundle. Returns its path."""
    raw_dir = Path(raw_dir).resolve()
    output = Path(output).resolve() if output else raw_dir / "kliore1980_raw.nc"
    table_csv = Path(table_csv).resolve() if table_csv else raw_dir.parent / "pioneer11_exit.csv"
    scalars_toml = raw_dir / "kliore1980_scalars.toml"
    with open(scalars_toml, "rb") as handle:
        scalars = tomllib.load(handle)
    source = scalars["source"]
    missing = [key for key in SEASON_KEYS if key not in source]
    if missing:
        raise ValueError(f"{scalars_toml}: [source] lacks {missing}")

    table = insert_level(read_table(table_csv), INSERTED_PRESSURE_Pa)
    latitude = float(scalars["latitude"]["planetocentric_deg"])
    oblateness = float(scalars["radius_mapping"]["oblateness"])
    radius = local_radius(table["radius_equatorial_m"], latitude, oblateness)

    variables = {target: (("level",), table[target],
                          {"units": units, "long_name": long_name, "value_source": "Table 1"})
                 for _, target, _, units, long_name in COLUMNS}
    variables["radius_local_m"] = (("level",), radius, {
        "units": "m", "long_name": "local radius, the equatorial mapping undone",
        "value_source": "radius_equatorial_m and scalars/radius_mapping",
    })
    table1 = xr.Dataset(variables, attrs={"level_order": "increasing pressure"})

    groups = {"table1": table1}
    groups.update(read_scalars(scalars_toml))
    data_inputs = [table_csv, scalars_toml]
    root = xr.Dataset(attrs={
        "title": "Kliore et al. (1980) Pioneer 11 exit, raw bundle",
        "source": source["citation"],
        "epoch": str(source["observation_date"]),
        "solar_longitude_deg": float(source["solar_longitude_deg"]),
        "subsolar_latitude_deg": float(source["subsolar_latitude_deg"]),
        "solar_longitude_source": str(source["solar_longitude_source"]),
        "input_hashes": cio.input_hashes(data_inputs, output),
    })
    cio.history_append(root, f"{TOOL}: Table 1 in SI by exact decimal conversion; the equatorial "
                             f"mapping undone at {latitude:g} deg with f = {oblateness:g}; one level "
                             f"inserted at {INSERTED_PRESSURE_Pa:g} Pa, linear in ln p")
    return cio.write(output, root, "raw", groups=groups, created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog=TOOL, description="Build the Kliore raw bundle.")
    parser.add_argument("--raw-dir", default="occul_data/kliore1980/raw",
                        help="directory holding kliore1980_scalars.toml")
    parser.add_argument("--output", default=None, help="defaults to <raw-dir>/kliore1980_raw.nc")
    args = parser.parse_args(argv)
    print(f"wrote {build(args.raw_dir, args.output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
