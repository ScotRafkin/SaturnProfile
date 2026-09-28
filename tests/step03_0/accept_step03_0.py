"""Acceptance checks for SPEC_03 v0.5 Step 0: the pressure grid and Eq. B1 along the local vertical.

Every check prints its measured value. Checks the specification did not ask for are labeled
"beyond the specification". The products before the rebuild were copied to
`reports/step03_0/before/` by the rebuild, so every "was" value below is read from a file, not
quoted. Refusal and exclusion cases are built in copies under this directory; nothing under
`occul_data/` is written by this script.

The candidate kind N (author ruling, 14 September 2026, option 1). The rebuilt inputs carry
`-dirty` until the acceptance commit, and `refrac` refuses such inputs. The candidate product is
therefore built by `refrac.product.build_product` in this process, on copies of the manifest and
the six inputs under `reports/step03_0/candidate/`, with the refusal relaxed in memory by
`relaxed.py`. The code and the rule are unchanged; the registered product comes from the clean
rebuild after acceptance.
"""

import math
import shutil
import sys
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relaxed  # noqa: E402  (applies the in-memory relaxation before refrac is used)

import numpy as np  # noqa: E402

from casspian.lib import control as ctl  # noqa: E402
from casspian.lib import io as cio  # noqa: E402
from casspian.lib.gravity import g_eff_radial  # noqa: E402
from casspian.refrac.anchor import geoid_setup  # noqa: E402
from casspian.refrac.product import build_product  # noqa: E402
from casspian.tools.lindal import build_raw  # noqa: E402
from casspian.tools.lindal.build_raw import PressureGridError, grid_candidates  # noqa: E402
from casspian.tools.plots import describe, render  # noqa: E402

D = Path("occul_data/lindal")
HERE = Path("reports/step03_0")
BEFORE = HERE / "before"
CANDIDATE = HERE / "candidate"
results = []

# The candidate product, built once, before any check reads it.
if CANDIDATE.exists():
    shutil.rmtree(CANDIDATE)
CANDIDATE.mkdir(parents=True)
for name in ("lindal_reduction.toml", "lindal_thermo.nc", "lindal_composition.nc",
             "lindal_geodesy.nc", "lindal_gravity.nc", "lindal_rotation.nc", "lindal_wind.nc"):
    shutil.copy2(D / name, CANDIDATE / name)
CANDIDATE_PRODUCT = build_product(CANDIDATE / "lindal_reduction.toml")
print(f"candidate product {CANDIDATE_PRODUCT}; inputs relaxed in memory: "
      f"{[Path(p).name for p in relaxed.RELAXED]}")
for line in describe(render(CANDIDATE_PRODUCT, CANDIDATE / "figures", "png", 150)):
    print(f"  {line}")
print()


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def root_of(tree):
    return tree.to_dataset(inherit=False)


# ---------------------------------------------------------------------------
# 1. The raw bundle: both pressures, the flag, the k sequence
# ---------------------------------------------------------------------------
raw = cio.read(D / "raw/lindal_raw.nc", "raw")
t1 = raw["table1"].dataset
p = t1["pressure_Pa"].values
pp = t1["pressure_printed_Pa"].values
k = t1["pressure_grid_index"].values
flag = t1["pressure_grid_applied"].values
grid_attrs = dict(raw["scalars/pressure_grid"].attrs)
history = str(raw.attrs["history"])
raw.close()
old_raw = cio.read(BEFORE / "raw/lindal_raw.nc", "raw")
old_t1 = old_raw["table1"].dataset.load()
old_raw.close()

finite_k = k[np.isfinite(k)].astype(int)
steps = np.diff(finite_k)
runs = []
for i, step in enumerate(steps):
    if runs and runs[-1][0] == step:
        runs[-1][1] = int(finite_k[i + 1])
    else:
        runs.append([int(step), int(finite_k[i + 1])])
# SPEC_03 v0.6 and SPEC_01 v0.24 (REVIEW_03_step0 ruling 1): 6 apart to 236, 4 apart to 268.
expected_runs = [[10, 130], [8, 194], [6, 236], [4, 268], [2, 310]]
top_expected = [19.9526, 25.1189, 31.6228, 39.8107, 50.1187, 63.0957, 79.4328, 100.0]
top_got = [float(f"{v:.6g}") for v in p[:8]]
one_bar = int(np.flatnonzero(pp == 100000.0)[0])
unchanged_columns = all(np.array_equal(t1[n].values, old_t1[n].values, equal_nan=True)
                        for n in ("temperature_K", "nh3_mole_fraction", "height_m"))
ok = (finite_k[0] == -70 and runs == expected_runs and int(flag.sum()) == 65
      and flag[-1] == 0 and np.isnan(k[-1]) and p[-1] == 129848.0 and pp[-1] == 129848.0
      and p[one_bar] == 100000.0 and top_got == top_expected
      and np.array_equal(pp, old_t1["pressure_Pa"].values) and unchanged_columns
      and "48 of them" in str(grid_attrs.get("basis", "")))
record(1, "raw bundle: pressure_printed_Pa, pressure_Pa on the grid, the flag, and the k sequence",
       ok,
       f"k starts at {finite_k[0]}; runs of (spacing, last k): {runs}\n"
       f"  expected {expected_runs}\n"
       f"snapped rows {int(flag.sum())} of {flag.size}; last row flag {int(flag[-1])}, k {k[-1]}, "
       f"pressure_Pa {p[-1]!r}, printed {pp[-1]!r}\n"
       f"top eight pressure_Pa to six figures {top_got}\n  expected {top_expected}\n"
       f"1 bar row pressure_Pa {p[one_bar]!r}\n"
       f"pressure_printed_Pa bit for bit the previous bundle's pressure_Pa: "
       f"{np.array_equal(pp, old_t1['pressure_Pa'].values)}\n"
       f"temperature, ammonia and height bit for bit unchanged: {unchanged_columns}\n"
       f"[pressure_grid] basis carries the v0.5 count of 48: {'48 of them' in str(grid_attrs.get('basis'))}\n"
       f"history lines on the grid:\n  "
       + "\n  ".join(l for l in history.splitlines() if "pressure grid" in l))


# ---------------------------------------------------------------------------
# 2 and 3. A printed value no grid value rounds to is refused; the exclusion list is data
# ---------------------------------------------------------------------------
def case_dir(name, old, new):
    case = HERE / name
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    for f in ("lindal_scalars.toml", "notes.md"):
        shutil.copy(D / "raw" / f, case / f)
    text = (D / "raw/lindal_table1.csv").read_text(encoding="utf-8")
    assert text.count(f"\n{old},") == 1, old
    (case / "lindal_table1.csv").write_text(text.replace(f"\n{old},", f"\n{new},"),
                                            encoding="utf-8")
    return case


case = case_dir("case_refused", "1000.00", "1001.00")
try:
    build_raw.build(case, case / "lindal_raw.nc")
    record(2, "a CSV copy with 1000.00 edited to 1001.00 is refused with the row named", False,
           "it was written")
except PressureGridError as exc:
    record(2, "a CSV copy with 1000.00 edited to 1001.00 is refused with the row named",
           "1001.00" in str(exc) and "row" in str(exc),
           f"grid candidates for 1001.00: {grid_candidates('1001.00', 100)}\nmessage: {exc}")

case = case_dir("case_excluded_moved", "1298.48", "1318.26")
build_raw.build(case, case / "lindal_raw.nc")
moved = cio.read(case / "lindal_raw.nc", "raw")
mt = moved["table1"].dataset
record(3, "a copy with the excluded row moved onto the grid (1318.26 mbar) is written with 66 "
          "snapped rows",
       int(mt["pressure_grid_applied"].values.sum()) == 66
       and float(mt["pressure_grid_index"].values[-1]) == 312.0,
       f"grid candidates for 1318.26: {grid_candidates('1318.26', 100)}\n"
       f"snapped rows {int(mt['pressure_grid_applied'].values.sum())}; last row k "
       f"{mt['pressure_grid_index'].values[-1]}, pressure_Pa {mt['pressure_Pa'].values[-1]!r}, "
       f"printed {mt['pressure_printed_Pa'].values[-1]!r}\n"
       f"the exclusion list is unchanged in the copy ({grid_attrs['excluded_printed_values_mbar']}); "
       "the code carries no row")
moved.close()

# ---------------------------------------------------------------------------
# 4. Kind T reads back with the rule and the printed pressures; T and D were rebuilt
# ---------------------------------------------------------------------------
thermo = cio.read(D / "lindal_thermo.nc", "thermo").load()
old_thermo = cio.read(BEFORE / "lindal_thermo.nc", "thermo").load()
geodesy = cio.read(D / "lindal_geodesy.nc", "geodesy").load()
old_geodesy = cio.read(BEFORE / "lindal_geodesy.nc", "geodesy").load()
raw_entry = cio.input_hash_entry(D / "raw/lindal_raw.nc", relative_to=Path.cwd())
ok = ("pressure_grid_rule" in thermo.attrs and "pressure_printed_Pa" in thermo
      and np.array_equal(thermo["pressure_Pa"].values, p)
      and np.array_equal(thermo["pressure_printed_Pa"].values, pp))
record(4, "kind T reads back with pressure_grid_rule and pressure_printed_Pa",
       ok,
       f"pressure_grid_rule = {thermo.attrs.get('pressure_grid_rule')!r}\n"
       f"pressure_Pa equals the bundle's grid values bit for bit: "
       f"{np.array_equal(thermo['pressure_Pa'].values, p)}\n"
       f"pressure_printed_Pa equals the bundle's printed values bit for bit: "
       f"{np.array_equal(thermo['pressure_printed_Pa'].values, pp)}\n"
       f"pressure_printed_Pa attributes {dict(thermo['pressure_printed_Pa'].attrs)}")
record("4b", "beyond the specification: kinds T and D rebuilt, not kept, because the raw bundle "
             "hash they record changed",
       thermo.attrs["created_at"] != old_thermo.attrs["created_at"]
       and geodesy.attrs["created_at"] != old_geodesy.attrs["created_at"]
       and thermo.attrs["raw_bundle"] != old_thermo.attrs["raw_bundle"],
       f"kind T created_at {old_thermo.attrs['created_at']} -> {thermo.attrs['created_at']}\n"
       f"kind D created_at {old_geodesy.attrs['created_at']} -> {geodesy.attrs['created_at']}\n"
       f"kind T raw_bundle {old_thermo.attrs['raw_bundle']}\n"
       f"              -> {thermo.attrs['raw_bundle']}\n"
       f"kind D raw_bundle {old_geodesy.attrs['raw_bundle']}\n"
       f"              -> {geodesy.attrs['raw_bundle']}\n"
       f"the raw bundle on disk: {raw_entry}\n"
       "The keep-by-content rule compares raw_bundle as content, so no change to it was needed "
       "(report finding).")

# ---------------------------------------------------------------------------
# 5. Kind C on kind T's levels; what the grid did to the ammonia fill
# ---------------------------------------------------------------------------
comp = cio.read(D / "lindal_composition.nc", "composition")
old_comp = cio.read(BEFORE / "lindal_composition.nc", "composition")
c_root, oc_root = root_of(comp), root_of(old_comp)
nh3, onh3 = c_root["x_NH3"].values, oc_root["x_NH3"].values
changed = np.flatnonzero(nh3 != onh3)
record(5, "beyond the specification: kind C levels are kind T's bit for bit; ammonia moves only "
          "where the fill interpolates or extrapolates",
       np.array_equal(c_root["pressure_Pa"].values, thermo["pressure_Pa"].values),
       f"levels equal: {np.array_equal(c_root['pressure_Pa'].values, thermo['pressure_Pa'].values)}\n"
       + "\n".join(f"level {i} ({pp[i] / 100:g} mbar): x_NH3 {onh3[i]:.9e} -> {nh3[i]:.9e}, "
                   f"relative {(nh3[i] - onh3[i]) / onh3[i]:+.3e}" for i in changed)
       + f"\nlevels where x_NH3 changed: {changed.tolist()}")
comp.close()
old_comp.close()

# ---------------------------------------------------------------------------
# 6. Kind N: the radius along the local vertical
# ---------------------------------------------------------------------------
N = cio.read(CANDIDATE_PRODUCT, "refractivity")
oN = cio.read(BEFORE / "lindal_refractivity.nc", "refractivity")
root, oroot = root_of(N), root_of(oN)
rec, orec = dict(N["reduction_record"].attrs), dict(oN["reduction_record"].attrs)
r, orad = root["radius_m"].values, oroot["radius_m"].values
dh = root["height_above_anchor_isobar_m"].values
r0 = float(root["anchor_isobar_radius_m"].values)
psi = math.radians(float(root["psi_deg"].values))
ka = int(np.flatnonzero(thermo["pressure_Pa"].values == 1.0e4)[0])
top_formula = r0 + dh[0] * math.cos(psi)
ok = (r[ka] == r0 and abs(r[0] - top_formula) < 1e-6
      and abs(r[0] - 58801571.0) < 0.05 and abs(r[-1] - 58412566.6) < 0.05
      and abs((r[0] - orad[0]) + 1317.2) < 0.05 and abs((r[-1] - orad[-1]) - 478.3) < 0.05
      and f"{math.cos(psi):.6f}" == "0.995406"
      and np.array_equal(dh, oroot["height_above_anchor_isobar_m"].values))
record(6, "radius_m at the anchor level equals r0 exactly and at the top r0 + (h_top - h_ref) cos psi",
       ok,
       f"anchor level {ka}: radius_m {r[ka]!r}, r0 {r0!r}, equal {r[ka] == r0}\n"
       f"top: radius_m {r[0]!r}, r0 + (h_top - h_ref) cos psi {top_formula!r}, difference "
       f"{r[0] - top_formula:.3e} m\n"
       f"top {r[0]:.1f} m (58,801,571.0), was {orad[0]:.1f}, change {r[0] - orad[0]:+.1f} (-1,317.2)\n"
       f"bottom {r[-1]:.1f} m (58,412,566.6), was {orad[-1]:.1f}, change {r[-1] - orad[-1]:+.1f} (+478.3)\n"
       f"cos psi {math.cos(psi):.6f} (0.995406); psi {math.degrees(psi)!r} degrees\n"
       f"height_above_anchor_isobar_m unchanged bit for bit: "
       f"{np.array_equal(dh, oroot['height_above_anchor_isobar_m'].values)}\n"
       f"long_name radius_m: {root['radius_m'].attrs['long_name']!r}\n"
       f"long_name height: {root['height_above_anchor_isobar_m'].attrs['long_name']!r}")

# ---------------------------------------------------------------------------
# 7. Number density and refractivity at the snapped levels
# ---------------------------------------------------------------------------
n, on = root["number_density_m3"].values, oroot["number_density_m3"].values
Nr, oNr = root["refractivity"].values, oroot["refractivity"].values
Rb, oRb = root["mean_refractivity_m3"].values, oroot["mean_refractivity_m3"].values
k1 = int(np.flatnonzero(thermo["pressure_Pa"].values == 1.0e5)[0])
p_ratio = p / pp
same_R = Rb == oRb
N_departure = np.abs(Nr / oNr / p_ratio - 1.0)
n_departure = np.abs(n / on / p_ratio - 1.0)
ok = (abs(n[0] / 1.041934e22 - 1) < 5e-7 and abs(n[1] / 1.270497e22 - 1) < 5e-7
      and n[k1] == on[k1] and float(N_departure[same_R].max()) < 1e-13
      and float(n_departure.max()) < 1e-13)
record(7, "n at the top 1.041934e22 and second level 1.270497e22; 1 bar unchanged; N moves by the "
          "fraction p moves",
       ok,
       f"top n {n[0]:.6e} (1.041934e22), was {on[0]:.6e}\n"
       f"second n {n[1]:.6e} (1.270497e22), was {on[1]:.6e}\n"
       f"1 bar n {n[k1]!r}, was {on[k1]!r}, equal {n[k1] == on[k1]}\n"
       f"largest |n_new / n_old / (p_grid / p_printed) - 1| over all levels: {n_departure.max():.2e}\n"
       f"largest |N_new / N_old / (p_grid / p_printed) - 1| where R_bar is unchanged "
       f"({int(same_R.sum())} levels): {N_departure[same_R].max():.2e}\n"
       + "".join(f"level {i} ({pp[i] / 100:g} mbar): R_bar changed by {Rb[i] / oRb[i] - 1:+.3e}, "
                 f"N by {Nr[i] / oNr[i] - 1:+.3e}, p by {p_ratio[i] - 1:+.3e}\n"
                 for i in np.flatnonzero(~same_R))
       + f"largest |p_grid / p_printed - 1|: {np.abs(p_ratio - 1).max():.3e} at "
         f"{pp[int(np.argmax(np.abs(p_ratio - 1)))] / 100:g} mbar")

# ---------------------------------------------------------------------------
# 8. The record carries the projection rule, its residual and the drift
# ---------------------------------------------------------------------------
drift = np.asarray(rec.get("latitude_drift_neglected_deg", [np.nan, np.nan]), dtype="float64")
residual = float(rec.get("radius_projection_residual_m", np.nan))
ok = (rec.get("radius_projection_rule") == "cos psi at the anchor"
      and abs(residual - 12.0) < 3.0
      and abs(drift[0] - 0.027) < 0.0015 and abs(drift[1] + 0.010) < 0.0015)
record(8, "reduction_record carries radius_projection_rule, radius_projection_residual_m and "
          "latitude_drift_neglected_deg",
       ok,
       f"radius_projection_rule = {rec.get('radius_projection_rule')!r}\n"
       f"radius_projection_residual_m = {residual!r} (about 12 m; the check allows 9 to 15)\n"
       f"latitude_drift_neglected_deg [top, bottom] = {drift.tolist()} (about +0.027 and -0.010; "
       "the check allows 0.0015 degrees)\n"
       f"radius_projection_note = {rec.get('radius_projection_note')!r}")

# ---------------------------------------------------------------------------
# 9. Everything in the anchoring decision is unchanged
# ---------------------------------------------------------------------------
expected_to_change = {"casspian_git_commit"} | {key for key in orec if key.startswith("input_sha256_")}
added = sorted(set(rec) - set(orec))
removed = sorted(set(orec) - set(rec))
differ = []
for key in orec:
    if key in expected_to_change or key not in rec:
        continue
    a, b = np.asarray(orec[key]), np.asarray(rec[key])
    if a.shape != b.shape or not np.array_equal(a, b):
        differ.append(key)
scalars = ["anchor_isobar_radius_m", "latitude_planetocentric_deg", "psi_deg",
           "anchor_isobar_radius_uncertainty_m", "latitude_planetocentric_uncertainty_deg"]
scalar_lines = [f"{s}: {float(oroot[s].values)!r} -> {float(root[s].values)!r}" for s in scalars]
scalars_equal = all(float(oroot[s].values) == float(root[s].values) for s in scalars)
frac = root["refractivity_uncertainty"].values / Nr
record(9, "r0, phi_c, psi, both polar radii, the asymmetry and every other reduction_record number "
          "unchanged bit for bit",
       not differ and not removed and scalars_equal
       and set(added) == {"radius_projection_rule", "radius_projection_residual_m",
                          "latitude_drift_neglected_deg", "radius_projection_note"},
       "\n".join(scalar_lines) + "\n"
       f"polar_radius_north_m {rec['polar_radius_north_m']!r}, south {rec['polar_radius_south_m']!r}, "
       f"asymmetry {rec['polar_asymmetry_m']!r}\n"
       f"record attributes compared: {len(orec) - len(expected_to_change)}; differing: "
       f"{differ if differ else 'none'}\n"
       f"expected to change and not compared: {sorted(expected_to_change)}\n"
       f"added: {added}\nremoved: {removed if removed else 'none'}\n"
       f"beyond the specification: N fractional uncertainty {frac.min():.6f} to {frac.max():.6f}")

# ---------------------------------------------------------------------------
# 10. F2's gravity along the profile
# ---------------------------------------------------------------------------
manifest = ctl.read_reduction_manifest(CANDIDATE / "lindal_reduction.toml")
inputs = ctl.load_reduction_inputs(manifest)
constants = geoid_setup(inputs, manifest).constants
phi_c = math.radians(float(root["latitude_planetocentric_deg"].values))
u = float(rec["u_at_anchor_ms"])
g_new = g_eff_radial(u, r, phi_c, *constants)
g_old = g_eff_radial(u, orad, phi_c, *constants)
change = np.abs(g_new / g_old - 1.0)
# SPEC_03 v0.6 (REVIEW_03_step0 ruling 2): the bound is 6e-5, the Newtonian 4.49e-5 divided by
# g_eff / g_N plus the centrifugal part.
record(10, "the gravity along the profile changes by under 6e-5",
       float(change.max()) < 6e-5,
       f"largest |g(r_new) / g(r_old) - 1| = {change.max():.3e} at level {int(np.argmax(change))} "
       f"({pp[int(np.argmax(change))] / 100:g} mbar); g_eff_radial with u(phi_c) = {u} m/s")
N.close()
oN.close()

# ---------------------------------------------------------------------------
# 11. Wind content unchanged; commits and hashes of every product
# ---------------------------------------------------------------------------
wind = cio.read(D / "lindal_wind.nc", "wind").load()
old_wind = cio.read(BEFORE / "lindal_wind.nc", "wind").load()
drop = {"created_at", "created_by", "casspian_git_commit", "history", "input_hashes",
        "control_file", "raw_bundle"}
w1, w0 = wind.copy(), old_wind.copy()
w1.attrs = {a: v for a, v in wind.attrs.items() if a not in drop}
w0.attrs = {a: v for a, v in old_wind.attrs.items() if a not in drop}
record("11a", "beyond the specification: the wind file, rebuilt because it records the raw bundle "
              "hash, is identical in content",
       w1.identical(w0), f"identical after dropping {sorted(drop)}: {w1.identical(w0)}")

rows = []
for name, kind in [("raw/lindal_raw.nc", "raw"), ("lindal_thermo.nc", "thermo"),
                   ("lindal_geodesy.nc", "geodesy"), ("lindal_gravity.nc", "gravity"),
                   ("lindal_rotation.nc", "rotation"), ("lindal_wind.nc", "wind"),
                   ("lindal_composition.nc", "composition")]:
    handle = cio.read(D / name, kind)
    commit = str(handle.attrs["casspian_git_commit"])
    handle.close()
    rows.append(f"`{Path(name).name}` {kind} {cio.sha256(D / name)} commit {commit}")
handle = cio.read(CANDIDATE_PRODUCT, "refractivity")
rows.append(f"candidate `lindal_refractivity.nc` refractivity {cio.sha256(CANDIDATE_PRODUCT)} "
            f"commit {handle.attrs['casspian_git_commit']}")
handle.close()
rows.append(f"`lindal_reduction.toml` manifest {cio.sha256(D / 'lindal_reduction.toml')} "
            f"(before: {cio.sha256(BEFORE / 'lindal_reduction.toml')})")
record("11b", "beyond the specification: commit and SHA-256 of every file now, provisional until "
              "the sweep", True, "\n".join(rows))

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
