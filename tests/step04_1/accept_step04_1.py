"""Acceptance checks for SPEC_04 v0.10 Step 1: `lib.geoid.through_anchor`, `lib.windfield`, and
the anchors' geopotential under the run's wind.

Every check prints its measured value; checks beyond the specification are labeled so. The anchor
and the run's inputs are the clean swept products of the Step 0 acceptance commit `1c9b310`, so
nothing is relaxed and no `-dirty` refusal is touched.

This step changes the kind `profile` product (it gains `u_column_ms` and its NaN companion), so by
SPEC_03 section 0 the acceptance runs `produce` in memory and compares against the registered
closure product on disk; the registered product is rebuilt on the clean tree at the sweep. Nothing
here writes inside the repository outside `reports/step04_1/`, which git ignores.

The cylinder-extended wind is not built here. Decision P at v0.10 pins the extension on the file's
reference level and builds it on the library's radial columns as a fixed point with the reference
surface, and the columns are Step 2's, so check 9's cylinder half and the first filing's check 14
moved to Step 2 with it (SPEC_04 section 12 ruling 2). Check 13 builds a minimal sheared field
instead, which is not that construction and is used for two measurements only.
"""

import math
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from casspian.forward import production as fp
from casspian.lib import control as ctl
from casspian.lib import geoid
from casspian.lib import io as cio
from casspian.lib import windfield as wf
from casspian.refrac.anchor import wind_of_latitude

HERE = Path("reports/step04_1")
CLOSURE = Path("forward/lindal_closure")
ANCHOR = Path("occul_data/lindal/lindal_refractivity.nc")
GAUGE_ISOBAR_PA = 1.0e4

HERE.mkdir(parents=True, exist_ok=True)
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def identical(a, b):
    """Bit for bit equality of two float64 arrays, not equality of value."""
    a = np.ascontiguousarray(np.asarray(a, dtype="float64"))
    b = np.ascontiguousarray(np.asarray(b, dtype="float64"))
    return a.shape == b.shape and np.array_equal(a.view(np.uint8), b.view(np.uint8))


def departure(a, b):
    """The largest relative difference of `a` from `b` where `b` is nonzero, and the largest in
    units in the last place of `b`. Reported beside the bit-identity, which says only True or
    False (REPORT_05_step0: on a platform other than the one that registered the product, the
    closure rerun differs from it in the last bit)."""
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    d = np.abs(a - b)
    nonzero = b != 0
    rel = float(np.max(d[nonzero] / np.abs(b[nonzero]))) if nonzero.any() else 0.0
    ulp = float(np.max(d / np.spacing(np.abs(b))))
    return rel, ulp


#: REVIEW_05_step0, the ruling on section 5: one bound on every platform.
RELATIVE_BOUND = 1e-14


def within(a, b, bound=RELATIVE_BOUND):
    """`a` agrees with `b` to `bound` relative where `b` is nonzero, and exactly where it is zero."""
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    if a.shape != b.shape:
        return False
    d = np.abs(a - b)
    nonzero = b != 0
    return bool(np.all(d[nonzero] <= bound * np.abs(b[nonzero])) and np.all(d[~nonzero] == 0))


# ---------------------------------------------------------------------------
# The run, the anchor and the geometry every check below shares
# ---------------------------------------------------------------------------
namelist = ctl.read_run_namelist(CLOSURE / "lindal_closure.toml")
inputs = ctl.load_run_inputs(namelist)
anchor = inputs.anchor
anchor_root = anchor.to_dataset(inherit=False)

phi_a = math.radians(float(anchor_root["latitude_planetocentric_deg"].values))
r_a = float(anchor_root["anchor_isobar_radius_m"].values)
r_levels = np.asarray(anchor_root["radius_m"].values, dtype="float64")
h_levels = np.asarray(anchor_root["height_above_anchor_isobar_m"].values, dtype="float64")
p_tab = np.asarray(anchor["inputs"]["thermo"]["pressure_Pa"].values, dtype="float64")
record_group = dict(anchor["reduction_record"].attrs)

constants = (
    float(inputs.rotation["angular_rate_rad_s"]),
    float(inputs.gravity["GM_m3s2"]),
    np.asarray(inputs.gravity["J"].values),
    np.asarray(inputs.gravity["degree"].values),
    float(inputs.gravity["normalization_radius_m"]),
)

closure_field = wf.WindField(inputs.wind)
wind_file_latitude = np.radians(
    np.asarray(inputs.wind["latitude_planetocentric_deg"].values, dtype="float64"))
wind_file_pressure = np.asarray(inputs.wind["pressure_Pa"].values, dtype="float64")
u_reference_file = np.asarray(inputs.wind["u_reference_ms"].values, dtype="float64")
order = np.argsort(wind_file_latitude)
wind_file_latitude = wind_file_latitude[order]
u_reference_file = u_reference_file[order]


def u_at_gauge(phi):
    """`u_total` at the gauge isobar's pressure, the wind the march is supplied with."""
    return closure_field.wind_at(np.asarray(phi, dtype="float64"), GAUGE_ISOBAR_PA)


WANTED_DEG = np.array([0.0, 10.0, 20.0, 45.0, 60.0])
surface = geoid.through_anchor(np.radians(WANTED_DEG), phi_a, r_a, u_at_gauge, *constants)

# ---------------------------------------------------------------------------
# 1. The equatorial radius
# ---------------------------------------------------------------------------
equator_expected = 60367000.0
d_equator = surface.equator_radius_m - equator_expected
record(1, "the surface through the anchor reproduces the equatorial radius to 0.01 m",
       abs(d_equator) <= 0.01,
       f"through ({math.degrees(phi_a):.6f} deg, {r_a:.2f} m): equator {surface.equator_radius_m:.4f} m, "
       f"expected {equator_expected:.1f}, difference {d_equator:+.2e} m\n"
       f"the reduction's own anchoring residual, its record: "
       f"{float(record_group['anchoring_residual_m']):+.2e} m")

# ---------------------------------------------------------------------------
# 2. The polar radii and the asymmetry
# ---------------------------------------------------------------------------
north_expected, south_expected = 54420643.5, 54449383.8
asymmetry_expected = 28740.37
d_north = surface.polar_north_m - north_expected
d_south = surface.polar_south_m - south_expected
d_asymmetry = surface.polar_asymmetry_m - asymmetry_expected
record(2, "the polar radii and the asymmetry (28,740.37 m) to 0.1 m",
       max(abs(d_north), abs(d_south), abs(d_asymmetry)) <= 0.1,
       f"north {surface.polar_north_m:.4f} m (expected {north_expected}), difference {d_north:+.4f}\n"
       f"south {surface.polar_south_m:.4f} m (expected {south_expected}), difference {d_south:+.4f}\n"
       f"asymmetry {surface.polar_asymmetry_m:.4f} m (expected {asymmetry_expected}), "
       f"difference {d_asymmetry:+.4f}\n"
       f"the reduction's record: north {float(record_group['polar_radius_north_m']):.4f}, "
       f"south {float(record_group['polar_radius_south_m']):.4f}, "
       f"asymmetry {float(record_group['polar_asymmetry_m']):.4f} m")

# ---------------------------------------------------------------------------
# 3. The stated latitudes
# ---------------------------------------------------------------------------
expected_at = {10.0: 60128613.0, 20.0: 59492075.7, 45.0: 57053249.5, 60.0: 55675128.9}
lines, worst = [], 0.0
for latitude, want in expected_at.items():
    got = float(surface.radius[int(np.flatnonzero(WANTED_DEG == latitude)[0])])
    worst = max(worst, abs(got - want))
    lines.append(f"{latitude:4.0f} N: {got:.4f} m, expected {want:.1f}, difference {got - want:+.4f}")
record(3, "the values at 10, 20, 45 and 60 N to 1 m", worst <= 1.0,
       "\n".join(lines) + f"\nlargest difference {worst:.4f} m")

# ---------------------------------------------------------------------------
# 4. Against the reduction's own march
# ---------------------------------------------------------------------------
reduction = geoid.wind_geoid(surface.march_latitude_rad, equator_expected, "equatorial_radius",
                             u_at_gauge, *constants,
                             tol_m=float(record_group.get("convergence_m", 1.0e-6)))
difference = np.asarray(reduction.radius, dtype="float64") - surface.march_radius_m
worst_node = int(np.argmax(np.abs(difference)))
record(4, "agrees with the reduction's own march (equatorial_radius rule) at every dense node to 0.01 m",
       float(np.max(np.abs(difference))) <= 0.01,
       f"{difference.size} dense nodes, step 0.05 deg; largest |difference| "
       f"{float(np.max(np.abs(difference))):.6e} m at "
       f"{math.degrees(surface.march_latitude_rad[worst_node]):+.3f} deg\n"
       f"the difference runs from {float(difference.min()):+.6e} to {float(difference.max()):+.6e} m, "
       f"one sign throughout: the two marches take the same steps and differ only in the constant\n"
       f"at the equator the difference is "
       f"{float(difference[int(np.flatnonzero(surface.march_latitude_rad == 0.0)[0])]):+.6e} m, which is "
       f"the anchoring residual the reduction recorded, "
       f"{-float(record_group['anchoring_residual_m']):+.6e} m: the anchor radius this surface is "
       f"marched through was read off that march, so it carries that residual, while this rerun of "
       f"wind_geoid converges to {float(reduction.anchor_residual_m):+.2e} m at the equator")

# ---------------------------------------------------------------------------
# 5. A second anchor entry on the same file
# ---------------------------------------------------------------------------
same = float(surface.residual(phi_a, r_a))
record(5, "a second anchor entry pointing at the same file gives a residual of zero at its latitude",
       same == 0.0,
       f"the second entry's measured radius is the first's, {r_a:.6f} m at "
       f"{math.degrees(phi_a):.6f} deg; residual {same:+.6e} m\n"
       f"the anchor latitude is a march node, so the surface is marched there and not interpolated")

# ---------------------------------------------------------------------------
# 6. A synthetic anchor moved by 1 km
# ---------------------------------------------------------------------------
moved = float(surface.residual(phi_a, r_a + 1000.0))
record(6, "a synthetic anchor entry with its radius moved by 1 km gives a residual of 1 km",
       abs(moved - 1000.0) <= 1.0e-6,
       f"measured radius {r_a + 1000.0:.6f} m at the same latitude; residual {moved:.9f} m, "
       f"departure from 1000 m {moved - 1000.0:+.3e}")

# ---------------------------------------------------------------------------
# 7. wind_at on the file's own nodes
# ---------------------------------------------------------------------------
u_total_file = np.asarray(inputs.wind["u_total_ms"].values, dtype="float64")[order]
node_phi, node_p = np.meshgrid(wind_file_latitude, wind_file_pressure, indexing="ij")
on_nodes = np.asarray(closure_field.wind_at(node_phi, node_p), dtype="float64")
exact = identical(on_nodes, u_total_file)
record(7, "wind_at reproduces the file's nodes exactly",
       exact,
       f"{u_total_file.size} nodes ({wind_file_latitude.size} latitudes by "
       f"{wind_file_pressure.size} pressures); bit for bit identical: {exact}; "
       f"largest |difference| {float(np.max(np.abs(on_nodes - u_total_file))):.3e} m/s\n"
       f"reference level {closure_field.reference_pressure_Pa:.1f} Pa: wind_at there equals "
       f"refrac.anchor.wind_of_latitude at the anchor bit for bit: "
       f"{identical(closure_field.reference_wind(np.array([phi_a])), wind_of_latitude(inputs.wind)(np.array([phi_a])))}")

# ---------------------------------------------------------------------------
# 8. The refusals outside coverage
# ---------------------------------------------------------------------------
low_p, high_p = closure_field.pressure_bounds_Pa
cases = [("latitude above the grid", math.radians(90.5), 1.0e4),
         ("latitude below the grid", math.radians(-90.5), 1.0e4),
         ("pressure above the grid", phi_a, high_p * 1.000001),
         ("pressure below the grid", phi_a, low_p * 0.999999)]
lines, refused = [], 0
for name, phi, p in cases:
    try:
        closure_field.wind_at(phi, p)
    except ValueError as error:
        refused += 1
        lines.append(f"{name}: refused, {str(error).splitlines()[0][:110]}")
    else:
        lines.append(f"{name}: NOT refused")
for name, phi, p in [("on the grid's corner, north pole at the deepest node",
                      float(wind_file_latitude[-1]), high_p)]:
    value = float(closure_field.wind_at(phi, p))
    lines.append(f"{name}: accepted, u = {value:.6f} m/s (the corner node itself)")
record(8, "wind_at refuses a point outside the file's coverage and accepts the grid's own corner",
       refused == len(cases), "\n".join(lines))

# ---------------------------------------------------------------------------
# 9. The wind along the anchor's column, the closure wind
# ---------------------------------------------------------------------------
# The cylinder-extended wind's half of this check, and check 14 of the first filing, moved to
# Step 2 with the construction they test: decision P at SPEC_04 v0.10 pins the extension on the
# file's reference level and builds it on the library's radial columns as a fixed point with the
# reference surface, and the columns are Step 2's (SPEC_04 section 12 ruling 2). Nothing was
# loosened: the closure half is unchanged and still bounded at 1e-3 m/s.
closure_column = np.asarray(
    closure_field.wind_at(np.full(p_tab.shape, phi_a), p_tab), dtype="float64")
gauge_level = int(np.flatnonzero(p_tab == namelist.gauge_isobar_Pa)[0])
closure_expected = 2.167
closure_worst = float(np.max(np.abs(closure_column - closure_expected)))
record(9, "the wind along the anchor's column under the closure wind matches the stated value to "
          "1e-3 m/s",
       closure_worst <= 1.0e-3,
       f"closure wind: {closure_column.min():.6f} to {closure_column.max():.6f} m/s over the "
       f"{p_tab.size} levels, expected {closure_expected} at every level; largest departure "
       f"{closure_worst:.3e} m/s\n"
       f"the file carries no shear, so the column is one value: max |u_shear| "
       f"{float(np.max(np.abs(inputs.wind['u_shear_ms'].values))):.1e} m/s\n"
       f"at the gauge isobar it equals the reduction's u_at_anchor_ms "
       f"{float(record_group['u_at_anchor_ms']):.11f} bit for bit: "
       f"{identical(closure_column[gauge_level], float(record_group['u_at_anchor_ms']))}\n"
       f"the cylinder-extended wind's half of this check is Step 2's at v0.10, by decision P")

# ---------------------------------------------------------------------------
# 10. The closure rerun through produce with the wind array
# ---------------------------------------------------------------------------
profile = fp.profile_from_anchor(anchor)
p_b = float(p_tab[0])
scalar = fp.produce(profile, inputs, gauge_level, p_b)
with_array = fp.produce(profile, inputs, gauge_level, p_b, u_column=closure_column)

fields = ("mean_refractivity_m3", "mean_molar_mass_kg_mol", "number_density_m3", "density_kg_m3",
          "layer_mass_Pa", "pressure_Pa", "temperature_K")
same_as_scalar = all(identical(getattr(scalar, name), getattr(with_array, name)) for name in fields)
same_geopotential = all(
    identical(getattr(scalar.geopotential, name), getattr(with_array.geopotential, name))
    for name in ("geopotential_m2s2", "layer_increments_m2s2", "g_magnitude_ms2", "psi_rad"))

product = cio.read(CLOSURE / "output/lindal_closure_profile.nc", "profile").to_dataset(inherit=False)
on_disk = {"pressure_Pa": with_array.pressure_Pa,
           "temperature_K": with_array.temperature_K,
           "number_density_m3": with_array.number_density_m3,
           "mean_refractivity_m3": with_array.mean_refractivity_m3,
           "mean_molar_mass_kg_mol": with_array.mean_molar_mass_kg_mol,
           "geopotential_m2s2": with_array.geopotential.geopotential_m2s2}
disk_same = {name: identical(product[name].values, value) for name, value in on_disk.items()}
disk_departure = {name: departure(value, product[name].values)
                  for name, value in on_disk.items()}
disk_within = {name: within(value, product[name].values) for name, value in on_disk.items()}
record(10, "the SPEC_03 closure rerun through produce with the wind array is bit-identical to the "
           "scalar path and agrees with the closure product to 1e-14 relative (REVIEW_05_step0)",
       same_as_scalar and same_geopotential and all(disk_within.values()),
       f"the wind array is the run's own wind along the column, {closure_column[0]:.11f} m/s at "
       f"every level, against the scalar {scalar.u_ms:.11f}\n"
       f"against the scalar path, bit for bit: {sorted(fields)} {same_as_scalar}; the "
       f"geopotential, its increments, |g_eff| and psi {same_geopotential}\n"
       f"against the registered product on disk (commit "
       f"{product.attrs.get('casspian_git_commit', '')[:12]}), bit for bit: "
       + "; ".join(f"{name} {value}" for name, value in disk_same.items()) + "\n"
       + f"within {RELATIVE_BOUND:g} relative of the registered product: {all(disk_within.values())}; "
       + "largest departure: "
       + "; ".join(f"{name} {rel:.3e} relative ({ulp:g} ulp)"
                   for name, (rel, ulp) in disk_departure.items()) + "\n"
       f"u_at_phi_c_ms keeps its meaning, the reference-level value {with_array.u_ms:.11f} m/s, and "
       f"u_column_ms holds {with_array.u_column_ms.size} levels")

# ---------------------------------------------------------------------------
# 11. Beyond the specification: the surface passes through the stated point exactly
# ---------------------------------------------------------------------------
at_anchor = float(geoid.radius_at(phi_a, surface.march_latitude_rad, surface.march_radius_m))
record(11, "beyond the specification: the surface passes through the stated point exactly, with no "
           "root find and no declared tolerance",
       at_anchor - r_a == 0.0,
       f"r0 at the anchor latitude {at_anchor:.6f} m, stated {r_a:.6f} m, difference "
       f"{at_anchor - r_a:+.1e} m\n"
       f"the equatorial_radius rule reaches its anchored radius by a secant on the polar start and "
       f"stops at a residual it reports, {float(reduction.anchor_residual_m):+.2e} m here and "
       f"{float(record_group['anchoring_residual_m']):+.2e} m in the reduction of record; this rule "
       f"has no residual to report")

# ---------------------------------------------------------------------------
# 12. Beyond the specification: the derivatives are the interpolant's
# ---------------------------------------------------------------------------
j = int(np.searchsorted(wind_file_latitude, phi_a, side="right") - 1)
step = 1.0e-7
mid_phi = 0.5 * (wind_file_latitude[j] + wind_file_latitude[j + 1])
d_phi, d_ln_p = closure_field.wind_derivatives(mid_phi, 1.0e4)
secant = float((closure_field.wind_at(mid_phi + step, 1.0e4)
                - closure_field.wind_at(mid_phi - step, 1.0e4)) / (2.0 * step))
# SPEC_07 v0.3 (decision L2): the slope in latitude is continuous across a node, so the two one
# sided values agree, and both are scipy's PCHIP derivative of the reference row there.
node_north, _ = closure_field.wind_derivatives(wind_file_latitude[j] + 1.0e-12, 1.0e4)
node_south, _ = closure_field.wind_derivatives(wind_file_latitude[j] - 1.0e-12, 1.0e4)
scipy_slope = float(PchipInterpolator(wind_file_latitude, u_reference_file).derivative()(
    wind_file_latitude[j]))
interval_north = interval_south = scipy_slope
record(12, "beyond the specification: wind_derivatives returns the interpolant's own slopes, and at "
           "a node the one sided slopes agree with each other and with scipy's PCHIP derivative "
           "(decision L2, SPEC_07)",
       abs(float(d_phi) - secant) <= 1.0e-6 * abs(secant)
       and abs(float(node_north) - scipy_slope) <= 1.0e-9 * abs(scipy_slope)
       and abs(float(node_south) - scipy_slope) <= 1.0e-9 * abs(scipy_slope)
       and float(d_ln_p) == 0.0,
       f"(du/dphi)_p inside the cell bracketing the anchor: {float(d_phi):.6f} m/s per radian, "
       f"central difference of wind_at {secant:.6f}\n"
       f"in m/s per degree that is {math.radians(float(d_phi)):.5f}, the "
       f"{math.degrees(wind_file_latitude[j]):.1f} to {math.degrees(wind_file_latitude[j + 1]):.1f} deg "
       f"interval where the anchor sits\n"
       f"at the {math.degrees(wind_file_latitude[j]):.1f} deg node: just north {float(node_north):.6f} "
       f"per radian, just south {float(node_south):.6f}, scipy's PCHIP derivative {interval_north:.6f}\n"
       f"(du/dln p)_phi is {float(d_ln_p):.1e} per unit ln p: the closure wind has no shear, "
       f"max |u_shear| {float(np.max(np.abs(inputs.wind['u_shear_ms'].values))):.1e} m/s")

# 13. Beyond the specification: wind_on_mesh and the pressure derivative, on a sheared field
# ---------------------------------------------------------------------------
# Both wind files in hand carry no shear, so `(du/dln p)_phi` is exactly zero on them and
# `wind_on_mesh` would be exercised on a field uniform in pressure. A minimal sheared field is
# built here for those two, as SPEC_04 section 1 allows ("synthetic files built inside the
# acceptance scripts"): `u_total` plus `a (ln p - ln p_ref) cos(phi)`, which keeps the poles at
# zero and is exactly linear in `ln p`. It is not the cylinder-extended wind, which decision P
# places at Step 2; it is the smallest field that makes these two checks mean something.
SHEAR_PER_LN_P = 1.5
sheared_grid = (np.asarray(inputs.wind["u_total_ms"].values, dtype="float64")[order]
                + SHEAR_PER_LN_P
                * (np.log(wind_file_pressure) - math.log(closure_field.reference_pressure_Pa))
                * np.cos(wind_file_latitude)[:, None])
sheared = inputs.wind.copy(deep=True)
sheared["latitude_planetocentric_deg"] = (
    ("latitude_planetocentric",), np.degrees(wind_file_latitude))
sheared["u_total_ms"] = (("latitude_planetocentric", "pressure"), sheared_grid)
sheared_field = wf.WindField(sheared)

columns = np.radians(np.array([0.0, 10.0, 30.805568, 60.0]))
pressure_map = np.tile(p_tab, (columns.size, 1))
on_mesh = np.asarray(sheared_field.wind_on_mesh(columns, pressure_map), dtype="float64")
one_by_one = np.array([[float(sheared_field.wind_at(phi, p)) for p in p_tab] for phi in columns])

_, sheared_d_ln_p = sheared_field.wind_derivatives(mid_phi, 1.0e4)
ln_step = 1.0e-6
secant_ln_p = float((sheared_field.wind_at(mid_phi, math.exp(math.log(1.0e4) + ln_step))
                     - sheared_field.wind_at(mid_phi, math.exp(math.log(1.0e4) - ln_step)))
                    / (2.0 * ln_step))
# SPEC_07 v0.3 (decision L2): the interpolant's own du/dln p, formed independently by scipy:
# each pressure row by PCHIP in latitude at mid_phi, then PCHIP in ln p. (PCHIP is not linear in
# the data, so the shear term's PCHIP alone is not the answer.)
_rows_mid = PchipInterpolator(wind_file_latitude, sheared_grid, axis=0)(mid_phi)
analytic = float(PchipInterpolator(np.log(wind_file_pressure), _rows_mid).derivative()(
    math.log(1.0e4)))
record(13, "beyond the specification: wind_on_mesh gives every node what wind_at gives it, and the "
           "pressure derivative is the interpolant's on a sheared field",
       on_mesh.shape == pressure_map.shape and identical(on_mesh, one_by_one)
       and abs(float(sheared_d_ln_p) - secant_ln_p) <= 1.0e-6 * abs(secant_ln_p)
       and abs(float(sheared_d_ln_p) - analytic) <= 1.0e-9 * abs(analytic),
       f"{columns.size} columns by {p_tab.size} levels, shape {on_mesh.shape}; bit for bit equal to "
       f"wind_at node by node: {identical(on_mesh, one_by_one)}\n"
       f"the mesh object arrives at Step 2; this takes the column latitudes and the pressure map, "
       f"which is what it needs of one\n"
       f"on the synthetic field, shear {SHEAR_PER_LN_P} m/s per unit ln p times cos(phi): "
       f"(du/dln p)_phi is {float(sheared_d_ln_p):.9f}, the central difference of wind_at "
       f"{secant_ln_p:.9f}, the field's own value {analytic:.9f}\n"
       f"the sum identity of the synthetic file is not rebuilt, so it is never written and never "
       f"read through the schema; it exists in memory for these two measurements only")

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
print()
passed = sum(1 for _, _, ok, _ in results if ok)
print(f"{passed} of {len(results)} checks pass")
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n" + "\n".join("        " + line
              for line in str(detail).splitlines())
              for n, d, ok, detail in results) + f"\n\n{passed} of {len(results)} checks pass\n",
    encoding="utf-8")
sys.exit(0 if passed == len(results) else 1)
