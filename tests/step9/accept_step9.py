"""Acceptance checks for SPEC_01 v0.17 Step 9, kinds T and D and the reduction manifest.

Every check reports its measured value.
"""

import sys
import tomllib
from pathlib import Path

import numpy as np

from casspian.lib import io as cio
from casspian.lib.schema import CasspianSchemaError

DIRECTORY = Path("occul_data/lindal")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


KINDS = {
    "lindal_thermo.nc": "thermo",
    "lindal_composition.nc": "composition",
    "lindal_geodesy.nc": "geodesy",
    "lindal_gravity.nc": "gravity",
    "lindal_rotation.nc": "rotation",
    "lindal_wind.nc": "wind",
}

# ---------------------------------------------------------------------------
# 1. All six files validate under their kinds
# ---------------------------------------------------------------------------
lines, ok = [], True
for name, kind in KINDS.items():
    path = DIRECTORY / name
    try:
        handle = cio.read(path, kind)
        commit = str(handle.attrs.get("casspian_git_commit", ""))
        created = str(handle.attrs.get("created_by", ""))
        handle.close()
        lines.append(f"  {name:26s} kind {kind:12s} ok   {created:34s} "
                     f"{'clean' if not commit.endswith('-dirty') else 'DIRTY'}")
    except CasspianSchemaError as exc:
        ok = False
        lines.append(f"  {name:26s} kind {kind:12s} REFUSED: {exc}")
record(1, "all six model input files validate under their declared kinds", ok,
       "\n".join(lines)
       + "\nA file marked DIRTY was written from an uncommitted tree and is to be rebuilt "
         "after the acceptance commit (SPEC_01 section 0).")

# ---------------------------------------------------------------------------
# 2. Kind T carries the profile and nothing else
# ---------------------------------------------------------------------------
thermo = cio.read(DIRECTORY / "lindal_thermo.nc", "thermo")
# SPEC_01 v0.23 Step 9: kind T also carries the printed pressures beside the grid values.
expected = {"pressure_Pa", "temperature_K", "height_m", "pressure_printed_Pa",
            "pressure_uncertainty_Pa", "temperature_uncertainty_K", "height_uncertainty_m"}
present = set(thermo.data_vars)
strays = sorted(v for v in present
                if str(v).startswith("x_") or "radius" in str(v) or "oblate" in str(v))
uncertainties_nan = {
    name: bool(np.all(np.isnan(thermo[name].values)))
    for name in ("pressure_uncertainty_Pa", "temperature_uncertainty_K", "height_uncertainty_m")
}
record(2, "kind T has exactly the three profile variables and their three companions, and no "
       "composition or geodesy content",
       present == expected and not strays and all(uncertainties_nan.values()),
       f"variables: {sorted(present)}\n"
       f"expected : {sorted(expected)}\n"
       f"exactly equal: {present == expected}\n"
       f"x_* or geodesy variables present: {strays if strays else 'none'}\n"
       f"the three companions are present and all NaN: {uncertainties_nan}\n"
       f"thermo_instance = {thermo.attrs.get('thermo_instance')!r} (not written since SPEC_R1), "
       f"vertical_coordinate = {thermo.attrs.get('vertical_coordinate')!r}, "
       f"latitude_planetocentric_absent_meaning = "
       f"{thermo.attrs.get('latitude_planetocentric_absent_meaning')!r}")

# ---------------------------------------------------------------------------
# 3. The anchor latitude and its provenance
# ---------------------------------------------------------------------------
lat = float(thermo.attrs["latitude_planetographic_deg"])
swath = np.asarray(thermo.attrs["latitude_planetographic_deg_swath"])
record(3, "the profile latitude is 36.3 degrees planetographic with its value_source",
       lat == 36.3 and "latitude_planetographic_deg_value_source" in thermo.attrs,
       f"latitude_planetographic_deg = {lat}\n"
       f"value_source = {thermo.attrs['latitude_planetographic_deg_value_source']!r}\n"
       f"swath = {swath.tolist()} degrees, uncertainty "
       f"{thermo.attrs['latitude_planetographic_deg_uncertainty']} "
       f"({thermo.attrs['latitude_planetographic_deg_uncertainty_kind']})\n"
       f"definition = {thermo.attrs['latitude_definition'][:110]}...")
thermo.close()

# ---------------------------------------------------------------------------
# 4. Kind D carries the two fitted surfaces
# ---------------------------------------------------------------------------
geodesy = cio.read(DIRECTORY / "lindal_geodesy.nc", "geodesy")
pressures = geodesy["surface_pressure_Pa"].values
expected_rows = {
    1.0e4: (60367.0, 54438.0, 0.09822),
    1.0e5: (60268.0, 54364.0, 0.09796),
}
rows, ok4 = [], True
for i, p in enumerate(pressures):
    r_eq = float(geodesy["radius_equatorial_m"].values[i]) / 1e3
    r_po = float(geodesy["radius_polar_m"].values[i]) / 1e3
    obl = float(geodesy["oblateness"].values[i])
    want = expected_rows.get(float(p))
    good = want is not None and (abs(r_eq - want[0]) < 1e-6 and abs(r_po - want[1]) < 1e-6
                                 and abs(obl - want[2]) < 1e-12)
    ok4 = ok4 and good
    rows.append(
        f"  {p:8.0f} Pa: r_eq {r_eq:9.3f} +- "
        f"{float(geodesy['radius_equatorial_uncertainty_m'].values[i]) / 1e3:.0f} km, "
        f"r_polar {r_po:9.3f} +- "
        f"{float(geodesy['radius_polar_uncertainty_m'].values[i]) / 1e3:.0f} km, "
        f"oblateness {obl:.5f} +- "
        f"{float(geodesy['oblateness_uncertainty'].values[i]):.5f}"
        + ("" if good else f"   <-- expected {want}")
    )
record(4, "kind D carries the two fitted surfaces with their uncertainties", ok4,
       "\n".join(rows)
       + f"\nfit_residual_m = {geodesy['fit_residual_m'].values.tolist()}, NaN because the "
         f"source states a range across its fits rather than a value per surface; the range "
         f"is in fit_residual_range_m = "
         f"{np.asarray(geodesy.attrs['fit_residual_range_m']).tolist()} m\n"
       f"fit_latitude_convention = {geodesy.attrs['fit_latitude_convention'][:80]}...")
geodesy.close()

# ---------------------------------------------------------------------------
# 5. The manifest parses and every path in it exists
# ---------------------------------------------------------------------------
manifest_path = DIRECTORY / "lindal_reduction.toml"
with open(manifest_path, "rb") as handle:
    manifest = tomllib.load(handle)
named = manifest["inputs"]
missing = [v for v in named.values() if not (DIRECTORY / v).exists()]
# SPEC_R1: the prefix and same-directory rules are retired; the check keeps what code needs.
record(5, "the reduction manifest parses, names six inputs, and every path exists",
       not missing and len(named) == 6,
       f"sections: {sorted(manifest)}\n"
       f"profile slug = {manifest['profile']['slug']!r}\n"
       f"inputs ({len(named)}): {list(named.values())}\n"
       f"missing: {missing if missing else 'none'}\n"
       f"anchor_isobar.pressure_Pa = {manifest['anchor_isobar']['pressure_Pa']}\n"
       f"anchor_latitude = {manifest['anchor_latitude']}\n"
       f"geoid = {manifest['geoid']}\n"
       f"output.product = {manifest['output']['product']!r}, which refrac writes (SPEC_02)")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
