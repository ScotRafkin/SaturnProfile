"""Acceptance checks for SPEC_02 v0.6 Step 4, the kind N product and `casspian-refrac`.

Every check reports its measured value. The product is written by the entry point, as a user
would run it; the refusal cases work on copies under this directory.
"""

import math
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib.schema import CasspianSchemaError
from casspian.refrac.anchor import freeze_anchor
from casspian.refrac.reduce import reduce_profile

HERE = Path("reports/step02_4")
HERE.mkdir(parents=True, exist_ok=True)
PROFILE = Path("occul_data/lindal")
MANIFEST = PROFILE / "lindal_reduction.toml"
PRODUCT = PROFILE / "lindal_refractivity.nc"
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


# ---------------------------------------------------------------------------
# 1. casspian-refrac writes the product, and it reads back as kind N
# ---------------------------------------------------------------------------
run = subprocess.run(["casspian-refrac", str(MANIFEST)], capture_output=True, text=True)
tree = cio.read(PRODUCT, "refractivity")
root = tree.dataset
groups = sorted(node.path for node in tree.subtree if node.path != "/")
record(1, "casspian-refrac writes lindal_refractivity.nc, which reads back as kind N, a DataTree",
       run.returncode == 0 and isinstance(tree, xr.DataTree)
       and root.attrs["casspian_kind"] == "refractivity",
       f"exit {run.returncode}: {run.stdout.strip()} {run.stderr.strip()}\n"
       f"read returns {type(tree).__name__}; casspian_kind {root.attrs['casspian_kind']!r}; "
       f"role {root.attrs['role']!r}; commit {root.attrs['casspian_git_commit']}\n"
       f"level variables: {sorted(v for v in root.data_vars if root[v].dims == ('level',))}\n"
       f"scalars: {sorted(v for v in root.data_vars if root[v].dims == ())}\n"
       f"coordinate: {list(root.coords)}\n"
       f"groups: {groups}\n"
       f"input_hashes:\n  " + str(root.attrs["input_hashes"]).replace(chr(10), chr(10) + "  "))

# ---------------------------------------------------------------------------
# 2. read refuses a copy missing inputs/wind
# ---------------------------------------------------------------------------
missing = HERE / "missing_wind.nc"
nodes = {node.path: node.to_dataset(inherit=False) for node in tree.subtree
         if node.path != "/inputs/wind"}
xr.DataTree.from_dict(nodes).to_netcdf(missing, engine="netcdf4")
try:
    cio.read(missing, "refractivity")
    refused_wind, message_wind = False, "not refused"
except CasspianSchemaError as exc:
    refused_wind, message_wind = True, str(exc)
record(2, "read refuses a copy missing the inputs/wind group", refused_wind,
       f"groups in the copy: {sorted(g for g in cio.group_paths(missing))}\n  {message_wind}")

# ---------------------------------------------------------------------------
# 3. read refuses a copy whose input_hashes entry for the thermo file is edited
# ---------------------------------------------------------------------------
edited = HERE / "edited_hash.nc"
shutil.copy(PRODUCT, edited)
with netCDF4.Dataset(edited, "a") as handle:
    lines = handle.getncattr("input_hashes").splitlines()
    changed = []
    for line in lines:
        if "lindal_thermo.nc" in line:
            head, hexhash = line.rsplit("sha256:", 1)
            line = head + "sha256:" + ("0" if hexhash[0] != "0" else "1") + hexhash[1:]
        changed.append(line)
    handle.setncattr("input_hashes", chr(10).join(changed))
try:
    cio.read(edited, "refractivity")
    refused_hash, message_hash = False, "not refused"
except CasspianSchemaError as exc:
    refused_hash, message_hash = True, str(exc)
record(3, "read refuses a copy whose input_hashes entry for the thermo file has been edited",
       refused_hash and "thermo" in message_hash,
       f"first hex digit of the thermo entry changed\n  {message_hash}")

# ---------------------------------------------------------------------------
# 4. Every embedded input group is identical to the file it came from
# ---------------------------------------------------------------------------
manifest = ctl.read_reduction_manifest(MANIFEST)
inputs = ctl.load_reduction_inputs(manifest)
lines, all_identical = [], True
for key in ctl.INPUT_KINDS:
    source = getattr(inputs, key)
    if isinstance(source, xr.DataTree):
        pairs = [(f"inputs/{key}" + node.path.rstrip("/"), node.to_dataset(inherit=False))
                 for node in source.subtree]
    else:
        pairs = [(f"inputs/{key}", source)]
    for path, original in pairs:
        embedded = tree[path].to_dataset(inherit=False)
        same = bool(embedded.identical(original))
        all_identical = all_identical and same
        lines.append(f"  {path:30s} identical {same}  ({len(original.variables)} variables, "
                     f"{len(original.attrs)} attributes)")
record(4, "every embedded input group is identical in variables and attributes to its file "
       "(xarray.Dataset.identical)", all_identical, "\n".join(lines))

# ---------------------------------------------------------------------------
# 5. The six scalars match Step 2 and Step 3
# ---------------------------------------------------------------------------
anchor = freeze_anchor(inputs, manifest)
red = reduce_profile(inputs, manifest, anchor)
expected = {
    "latitude_planetocentric_deg": math.degrees(anchor.phi_c_rad),
    "psi_deg": math.degrees(anchor.psi_rad),
    "latitude_planetographic_deg": float(inputs.thermo.attrs["latitude_planetographic_deg"]),
    "anchor_isobar_pressure_Pa": float(manifest.anchor_isobar_Pa),
    "anchor_isobar_radius_m": anchor.r0_m,
    "anchor_isobar_height_m": red.h_ref_m,
}
lines, all_equal = [], True
for name, value in expected.items():
    stored = float(root[name].values)
    equal = stored == value
    all_equal = all_equal and equal
    companion = [v for v in root.data_vars if v.startswith(name.rsplit("_", 1)[0] + "_uncertainty")]
    shown = ", ".join(f"{c} = {float(root[c].values):.6g} ({root[c].attrs['uncertainty_kind']})"
                      for c in companion)
    lines.append(f"  {name:30s} {stored!r:>24}  expected {value!r:>24}  equal {equal}; {shown}")
levels_equal = (np.array_equal(root["radius_m"].values, red.radius_m)
                and np.array_equal(root["refractivity"].values, red.refractivity)
                and np.array_equal(root["number_density_m3"].values, red.number_density_m3))
record(5, "the six scalars match Step 2 and Step 3, bit for bit",
       all_equal and levels_equal,
       "\n".join(lines)
       + f"\nper level radius_m, number_density_m3 and refractivity identical to a fresh "
         f"reduction: {levels_equal}")

# ---------------------------------------------------------------------------
# 6. The SPEC_00 section 2.2 listing is complete for the profile
# ---------------------------------------------------------------------------
listing = ["raw/lindal_table1.csv", "raw/lindal_scalars.toml", "raw/lindal_raw.nc",
           "raw/notes.md", "lindal_build.toml", "lindal_thermo.nc", "lindal_composition.nc",
           "lindal_geodesy.nc", "lindal_gravity.nc", "lindal_rotation.nc", "lindal_wind.nc",
           "lindal_reduction.toml", "lindal_refractivity.nc"]
present = {name: (PROFILE / name).exists() for name in listing}
record(6, "lindal_refractivity.nc is in occul_data/lindal/ and the SPEC_00 section 2.2 listing "
       "is complete", all(present.values()),
       "\n".join(f"  {name:28s} {'present' if ok else 'MISSING'}" for name, ok in present.items()))

# ---------------------------------------------------------------------------
# 7. Beyond the specification: the companions and the reduction record
# ---------------------------------------------------------------------------
rec = tree["reduction_record"].attrs
lines = []
for name in sorted(v for v in root.data_vars if "_uncertainty" in v):
    attrs = root[name].attrs
    value = np.atleast_1d(root[name].values)
    shown = "NaN" if np.all(np.isnan(value)) else f"{np.nanmin(value):.6g} to {np.nanmax(value):.6g}"
    lines.append(f"  {name:42s} {shown:26s} {attrs['uncertainty_kind']}; included "
                 f"[{attrs.get('uncertainty_terms_included', '(copied)')}] unstated "
                 f"[{attrs.get('uncertainty_terms_unstated', '')}]"
                 + (f"; conversions [{attrs['uncertainty_kind_conversions']}]"
                    if attrs.get("uncertainty_kind_conversions") else ""))
spread_n = float(rec["anchor_rule_spread_r0_north_pole_m"])
spread_s = float(rec["anchor_rule_spread_r0_south_pole_m"])
r0 = float(root["anchor_isobar_radius_m"].values)
needed = ["fixed_point_iterates_deg", "anchor_rule", "north_polar_start_m",
          "polar_radius_north_m", "polar_radius_south_m", "polar_asymmetry_m",
          "anchoring_residual_m", "nowind_latitude_planetocentric_deg",
          "nowind_anchor_isobar_radius_m", "anchor_rule_spread_r0_north_pole_m",
          "anchor_rule_spread_r0_south_pole_m", "partial_dr0_dr_anchor",
          "partial_dphi_c_dphi_g", "codata_release", "casspian_version", "casspian_git_commit",
          "fixed_point_closure_deg"]
absent = [n for n in needed if n not in rec]
record(7, "beyond the specification: every companion's attributes, and the reduction_record",
       not absent and float(root["radius_uncertainty_m"].values[0]) == float(
           red.companions["radius_uncertainty_m"].value[0]),
       "\n".join(lines)
       + f"\nreduction_record carries {len(rec)} attributes; required ones absent: {absent}\n"
         f"  iterates {np.round(np.asarray(rec['fixed_point_iterates_deg']), 7).tolist()}\n"
         f"  closure {float(rec['fixed_point_closure_deg']):.2e} deg; polar asymmetry "
         f"{float(rec['polar_asymmetry_m']) / 1e3:.3f} km; anchoring residual "
         f"{float(rec['anchoring_residual_m']):.2e} m\n"
         f"  no wind phi_c {float(rec['nowind_latitude_planetocentric_deg']):.6f} deg, r0 "
         f"{float(rec['nowind_anchor_isobar_radius_m']) / 1e3:.3f} km\n"
         f"  anchor rule spread at phi_c: north_pole {spread_n / 1e3:.3f} km "
         f"({(spread_n - r0) / 1e3:+.3f}), south_pole {spread_s / 1e3:.3f} km "
         f"({(spread_s - r0) / 1e3:+.3f}), against mean_polar_radius {r0 / 1e3:.3f} km\n"
         f"  partials: dr0/dr_anchor {float(rec['partial_dr0_dr_anchor']):.6f}, dphi_c/dphi_g "
         f"{float(rec['partial_dphi_c_dphi_g']):.6f}; closure "
         f"{rec['closure_rule']} {rec['closure_species']}")

# ---------------------------------------------------------------------------
# 8. Beyond the specification: no warning while every input on disk matches its recorded hash
# ---------------------------------------------------------------------------
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    cio.read(PRODUCT, "refractivity").close()
disk_warnings = [str(w.message) for w in caught if "input_hashes" in str(w.message)]
record(8, "beyond the specification: reading the product raises no input_hashes warning while "
       "every input on disk matches (SPEC_00 section 8 warns, it does not refuse)",
       not disk_warnings,
       f"input_hashes warnings with every input unchanged: {disk_warnings}")

# ---------------------------------------------------------------------------
# 9. SPEC_02 v0.7: the anchor rule spread by a full fixed point rerun under each rule
# ---------------------------------------------------------------------------
lat_n = float(rec["anchor_rule_spread_latitude_planetocentric_north_pole_deg"])
lat_s = float(rec["anchor_rule_spread_latitude_planetocentric_south_pole_deg"])
phi_c_deg = float(root["latitude_planetocentric_deg"].values)
record(9, "v0.7: the anchor rule spread matches the expected values (r0 58,537.7 and 58,502.1 km "
       "+- 0.1; phi_c about 30.8020 and 30.8079 deg)",
       abs(spread_n / 1e3 - 58537.7) <= 0.1 and abs(spread_s / 1e3 - 58502.1) <= 0.1
       and abs(lat_n - 30.8020) <= 0.0005 and abs(lat_s - 30.8079) <= 0.0005,
       f"north_pole: r0 {spread_n / 1e3:.3f} km ({(spread_n - r0) / 1e3:+.3f}), phi_c {lat_n:.6f} deg "
       f"({lat_n - phi_c_deg:+.6f})\n"
       f"south_pole: r0 {spread_s / 1e3:.3f} km ({(spread_s - r0) / 1e3:+.3f}), phi_c {lat_s:.6f} deg "
       f"({lat_s - phi_c_deg:+.6f})\n"
       f"mean_polar_radius: r0 {r0 / 1e3:.3f} km, phi_c {phi_c_deg:.6f} deg\n"
       f"half range {(spread_n - spread_s) / 2e3:.3f} km; the latitude held v0.6 values were "
       f"+17.506 and -17.475 km\n"
       f"note: {rec['anchor_rule_spread_note']}")

tree.close()
print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
