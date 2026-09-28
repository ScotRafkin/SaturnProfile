"""Acceptance checks for SPEC_01 Step 2, the Lindal raw bundle.

Prints the actual result of each check, then dumps every attribute of every `scalars` group so
that a second reader can check the transcription against the paper.
"""

import sys
from pathlib import Path

import numpy as np

from casspian.lib import io as cio
from casspian.lib.schema import CasspianSchemaError

BUNDLE = Path("occul_data/lindal/raw/lindal_raw.nc")

results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


tree = cio.read(BUNDLE, "raw")
t1 = tree["table1"].dataset

p = t1["pressure_Pa"].values
temp = t1["temperature_K"].values
nh3 = t1["nh3_mole_fraction"].values
h = t1["height_m"].values

# ---------------------------------------------------------------------------
record(1, "66 levels, ordered by increasing pressure",
       p.size == 66 and bool(np.all(np.diff(p) > 0)),
       f"level size = {p.size}; strictly increasing in pressure = {bool(np.all(np.diff(p) > 0))}\n"
       f"level_order attribute = {t1.attrs.get('level_order')!r}; "
       f"rows reordered from file = {t1.attrs.get('rows_reordered_from_file')!r}")

# ---------------------------------------------------------------------------
one_bar = int(np.argmin(np.abs(p - 100000.0)))
got = (p[one_bar], temp[one_bar], h[one_bar])
record(2, "the 1 bar row reads exactly 100000.0 Pa, 134.8 K, 0.0 m",
       got == (100000.0, 134.8, 0.0),
       f"index {one_bar}: pressure_Pa = {got[0]!r}, temperature_K = {got[1]!r}, "
       f"height_m = {got[2]!r}\nexact equality against the specified values: "
       f"{got == (100000.0, 134.8, 0.0)}")

# ---------------------------------------------------------------------------
# SPEC_01 v0.23: pressure_Pa is on the declared grid and pressure_printed_Pa is what Table I
# prints; the top row is 19.9526 Pa on the grid (to six figures) and 20.0 Pa as printed.
pp = t1["pressure_printed_Pa"].values
top = (p[0], temp[0], h[0])
record(3, "the top row reads 19.9526 Pa on the grid (20.0 Pa printed), 138.7 K, 376700.0 m",
       float(f"{top[0]:.6g}") == 19.9526 and pp[0] == 20.0
       and (top[1], top[2]) == (138.7, 376700.0),
       f"pressure_Pa = {top[0]!r}, pressure_printed_Pa = {pp[0]!r}, temperature_K = {top[1]!r}, "
       f"height_m = {top[2]!r}")

bottom = (p[-1], temp[-1], h[-1])
record(4, "the bottom row reads 129848.0 Pa, 146.2 K, -14100.0 m",
       bottom == (129848.0, 146.2, -14100.0),
       f"pressure_Pa = {bottom[0]!r}, temperature_K = {bottom[1]!r}, height_m = {bottom[2]!r}")

# ---------------------------------------------------------------------------
finite = np.isfinite(nh3)
count = int(finite.sum())
first, last = int(np.argmax(finite)), int(len(finite) - 1 - np.argmax(finite[::-1]))
gap_interior = float(nh3[int(np.argmin(np.abs(p - 104713.0)))])
gap_end = float(nh3[int(np.argmin(np.abs(p - 129848.0)))])
ok = (
    count == 9
    and pp[first] == 83176.0 and nh3[first] == 2.6e-6
    and pp[last] == 125893.0 and nh3[last] == 66.9e-6
    and np.isnan(gap_interior) and np.isnan(gap_end)
)
record(5, "exactly nine finite NH3 values, with the two stated gaps left as NaN", ok,
       f"finite count = {count}\n"
       f"first finite: {p[first]!r} Pa ({p[first]/100:g} mbar) = {nh3[first]!r} "
       f"({nh3[first]*1e6:g} ppm)\n"
       f"last  finite: {p[last]!r} Pa ({p[last]/100:g} mbar) = {nh3[last]!r} "
       f"({nh3[last]*1e6:g} ppm)\n"
       f"interior gap at 1047.13 mbar = {gap_interior} (NaN required)\n"
       f"end gap at 1298.48 mbar = {gap_end} (NaN required)\n"
       f"no filled value appears in the raw bundle")

# ---------------------------------------------------------------------------
lat = tree["scalars/latitude"].attrs
swath = np.asarray(lat.get("swath_planetographic_deg"))
record(6, "scalars/latitude carries planetographic_deg = 36.3 and the swath",
       float(lat.get("planetographic_deg")) == 36.3
       and swath.tolist() == [36.3, 36.7],
       f"planetographic_deg = {lat.get('planetographic_deg')!r}\n"
       f"swath_planetographic_deg = {swath.tolist()}\n"
       f"value_source = {lat.get('value_source')!r}")

# ---------------------------------------------------------------------------
record(7, "read as kind raw succeeds", tree is not None,
       f"casspian_kind = {tree.attrs.get('casspian_kind')!r}, "
       f"schema version = {tree.attrs.get('casspian_schema_version')!r}, "
       f"created_by = {tree.attrs.get('created_by')!r}\n"
       f"casspian_git_commit = {tree.attrs.get('casspian_git_commit')!r}")
try:
    cio.read(BUNDLE, "thermo")
    record(8, "read as kind thermo refuses", False, "it did not refuse")
except CasspianSchemaError as exc:
    record(8, "read as kind thermo refuses", True, str(exc))

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("Every attribute of every scalars group, for checking against the paper")
print("=" * 78)
for path in cio.group_paths(BUNDLE):
    if not path.startswith("scalars"):
        continue
    node = tree[path]
    print(f"\n[{path}]")
    for key, value in node.attrs.items():
        rendered = np.asarray(value).tolist() if isinstance(value, np.ndarray) else value
        print(f"  {key} = {rendered!r}")

print()
print("=" * 78)
print("Root globals")
print("=" * 78)
for key, value in tree.attrs.items():
    print(f"  {key} = {value!r}")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
tree.close()
sys.exit(1 if failed else 0)
