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
       f"thermo_instance = {thermo.attrs['thermo_instance']!r}, "
       f"vertical_coordinate = {thermo.attrs['vertical_coordinate']!r}, "
       f"latitude_planetocentric_absent_meaning = "
       f"{thermo.attrs['latitude_planetocentric_absent_meaning']!r}")

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
unprefixed = [v for v in named.values() if not v.startswith("lindal_")]
outside = [v for v in named.values() if "/" in v or "\\" in v or v.startswith("..")]
record(5, "the reduction manifest parses, names six inputs, and every path exists",
       not missing and not unprefixed and not outside and len(named) == 6,
       f"sections: {sorted(manifest)}\n"
       f"profile slug = {manifest['profile']['slug']!r}\n"
       f"inputs ({len(named)}): {list(named.values())}\n"
       f"missing: {missing if missing else 'none'}; unprefixed: "
       f"{unprefixed if unprefixed else 'none'}; pointing outside the directory: "
       f"{outside if outside else 'none'}\n"
       f"anchor_isobar.pressure_Pa = {manifest['anchor_isobar']['pressure_Pa']}\n"
       f"anchor_latitude = {manifest['anchor_latitude']}\n"
       f"geoid = {manifest['geoid']}\n"
       f"output.product = {manifest['output']['product']!r}, which refrac writes (SPEC_02)")

# ---------------------------------------------------------------------------
# 6. The directory matches SPEC_00 section 2.2
# ---------------------------------------------------------------------------
expected_top = {
    "lindal_build.toml", "lindal_thermo.nc", "lindal_composition.nc", "lindal_geodesy.nc",
    "lindal_gravity.nc", "lindal_rotation.nc", "lindal_wind.nc", "lindal_reduction.toml",
}
expected_raw = {"lindal_table1.csv", "lindal_scalars.toml", "lindal_raw.nc", "notes.md"}
actual_top = {p.name for p in DIRECTORY.iterdir() if p.is_file()}
actual_raw = {p.name for p in (DIRECTORY / "raw").iterdir() if p.is_file()}
not_prefixed = sorted(n for n in actual_top | actual_raw if not n.startswith("lindal_"))
# SPEC_02 Step 4 writes the kind N product into this directory, the last entry of the SPEC_00
# section 2.2 listing, so from then on it is expected here rather than counted as unexpected.
PRODUCT_N = "lindal_refractivity.nc"
record(6, "the directory listing matches SPEC_00 section 2.2",
       actual_top - {PRODUCT_N} == expected_top and actual_raw == expected_raw,
       f"top level ({len(actual_top)}): {sorted(actual_top)}\n"
       f"  expected: {sorted(expected_top)}, plus {PRODUCT_N} once SPEC_02 Step 4 has run\n"
       f"  missing {sorted(expected_top - actual_top) or 'none'}; "
       f"unexpected {sorted(actual_top - expected_top - {PRODUCT_N}) or 'none'}\n"
       f"raw/ ({len(actual_raw)}): {sorted(actual_raw)}\n"
       f"  missing {sorted(expected_raw - actual_raw) or 'none'}; "
       f"unexpected {sorted(actual_raw - expected_raw) or 'none'}\n"
       f"files without the lindal_ prefix: {not_prefixed}. SPEC_00 section 2.2 lists "
       f"notes.md unprefixed in its own directory listing, so that is the specification's "
       f"own exception rather than a defect.\n"
       f"{PRODUCT_N}: "
       f"{'present, written by refrac at SPEC_02 Step 4' if PRODUCT_N in actual_top else 'absent, as before SPEC_02 Step 4'}.")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
