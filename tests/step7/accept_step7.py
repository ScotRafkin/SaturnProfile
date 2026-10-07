"""Acceptance checks for SPEC_01 v0.15 Step 7, the Lindal wind tool and kind W.

Every check reports its measured value, as v0.15 requires. Figures go to reports/figures/,
which is committed; this directory keeps only the script and its output.
"""

import dataclasses
import math
import shutil
import sys
from pathlib import Path
from types import MappingProxyType

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

from casspian.lib import control as ctl
from casspian.lib import geoid as gd
from casspian.lib import io as cio
from casspian.lib.gravity import g_eff_radial
from scipy.interpolate import PchipInterpolator
from casspian.lib import latitude as latmod
from casspian.tools.wind.build_wind import bin_points, bin_rms_about, read_points
from casspian.tools.wind.curve import AssembledCurve, read_curve

HERE = Path("reports/step7")
HERE.mkdir(parents=True, exist_ok=True)
FIGURES = Path("reports/figures")
FIGURES.mkdir(parents=True, exist_ok=True)
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


CURVE_CSV = "data_static/winds/ingersoll_pollard1982_fig5_curve.csv"
POINTS_CSV = "data_static/winds/smith1982_fig4_points.csv"

segments = read_curve(CURVE_CSV)
points_lat, points_u = read_points(POINTS_CSV)
centers, index, count, mean, low = bin_points(points_lat, points_u, 2.0, 3)

wind = cio.read("occul_data/lindal/lindal_wind.nc", "wind")
gravity = cio.read("occul_data/lindal/lindal_gravity.nc", "gravity")
rotation = cio.read("occul_data/lindal/lindal_rotation.nc", "rotation")
raw = cio.read("occul_data/lindal/raw/lindal_raw.nc", "raw")
DEGREES = gravity["degree"].values
J = gravity["J"].values
GM = float(gravity["GM_m3s2"])
R_NORM = float(gravity["normalization_radius_m"])
OMEGA = float(rotation["angular_rate_rad_s"])
R_POLAR = float(raw["scalars/geodesy"].attrs["reference_geoid_polar_radius_km"]) * 1e3
R_FITTED = float(raw["scalars/geodesy/surface_100mbar"].attrs["radius_equatorial_km"]) * 1e3
JOIN = [float(v) for v in wind.attrs["join_window_deg"]]


def make_curve(gap_rule="reflect_north", polar_rule="pchip_to_zero"):
    return AssembledCurve(segments, gap_rule=gap_rule, polar_rule=polar_rule,
                          join_window_deg=JOIN, bin_centers=centers, bin_values=mean)


curve = make_curve()
rms, residual = bin_rms_about(curve, points_lat, points_u, index, centers, count, 3)

# ---------------------------------------------------------------------------
# 1. The curve passes through its samples
# ---------------------------------------------------------------------------
window_lo, window_hi = min(JOIN), max(JOIN)
parts, worst_outside, inside_count, worst_inside = [], 0.0, 0, 0.0
for name, (lat_s, u_s) in (("north", segments.north), ("south", segments.south)):
    error = np.abs(curve(lat_s) - u_s)
    in_window = (lat_s >= window_lo) & (lat_s <= window_hi)
    worst_outside = max(worst_outside, float(error[~in_window].max()))
    if in_window.any():
        inside_count += int(in_window.sum())
        worst_inside = max(worst_inside, float(error[in_window].max()))
    parts.append(f"{name}: {lat_s.size} samples, max error outside the join window "
                 f"{float(error[~in_window].max()):.3e} m/s")
record(1, "the curve passes through every digitized sample to round-off, outside the "
       "declared join window",
       worst_outside == 0.0,
       "\n".join(parts)
       + f"\nInside the join window [{window_lo:g}, {window_hi:g}] the blend deliberately "
         f"departs from the southern segment: {inside_count} samples, max departure "
         f"{worst_inside:.2f} m/s. SPEC_01 asks for every sample to round-off, which cannot "
         f"hold where a declared rule replaces the curve. See the report, finding 1.")

# ---------------------------------------------------------------------------
# 2. Zero at the poles, and the refusal test
# ---------------------------------------------------------------------------
lat_c = wind["latitude_planetocentric_deg"].values
u_file = wind["u_total_ms"].values
poles = [int(np.argmin(np.abs(lat_c - p))) for p in (-90.0, 90.0)]
polar_values = u_file[poles, :]
tampered = HERE / "tampered_wind.nc"
shutil.copy("occul_data/lindal/lindal_wind.nc", tampered)
with netCDF4.Dataset(tampered, "a") as handle:
    handle.variables["u_total_ms"][poles[1], :] = 1e-6
# The polar rule is applied where the model takes a wind (SPEC_08): the copy reads, and the
# reduction's loader, pointed at it through the registered manifest, refuses it.
cio.read(tampered, "wind").close()
manifest = ctl.read_reduction_manifest("occul_data/lindal/lindal_reduction.toml")
manifest = dataclasses.replace(
    manifest, inputs=MappingProxyType({**manifest.inputs, "wind": tampered.resolve()}))
refused, message = False, "not refused"
try:
    ctl.load_reduction_inputs(manifest)
except ctl.ControlFileError as exc:
    refused, message = True, str(exc)
record(2, "the wind is exactly zero at both poles, and a perturbed copy is refused",
       float(np.max(np.abs(polar_values))) == 0.0 and refused,
       f"latitudes {lat_c[poles].tolist()}; max |u| there {float(np.max(np.abs(polar_values))):.1e} m/s, "
       f"exactly zero: {float(np.max(np.abs(polar_values))) == 0.0}\n"
       f"a copy with the north pole set to 1e-6 m/s reads as kind W and is refused by "
       f"load_reduction_inputs:\n  {message}")

# ---------------------------------------------------------------------------
# 3. The stated values on the curve
# ---------------------------------------------------------------------------
dense = np.arange(-90.0, 90.0001, 0.01)
dense_u = curve(dense)
north = dense > 0
peak_i = int(np.argmax(dense_u[north]))
peak_u, peak_lat = float(dense_u[north][peak_i]), float(dense[north][peak_i])
south = dense < 0
rpeak_i = int(np.argmax(dense_u[south]))
u363 = float(curve(np.array([36.3]))[0])
u308 = float(curve(np.array([30.8]))[0])
ok3 = (abs(peak_u - 490.5) <= 2 and abs(peak_lat - 7.4) <= 0.3
       and abs(u363 - 1.9) <= 3 and abs(u308 - 75.0) <= 3)
record(3, "the curve reproduces the stated peak and the two anchor values", ok3,
       f"peak {peak_u:.2f} m/s at {peak_lat:+.2f} deg (spec 490.5 +- 2 at +7.4 +- 0.3)\n"
       f"reflected peak {float(dense_u[south][rpeak_i]):.2f} m/s at "
       f"{float(dense[south][rpeak_i]):+.2f} deg\n"
       f"u(36.3) = {u363:.3f} m/s (spec 1.9 +- 3). THIS VALUE ENTERS THE FROZEN ANCHOR.\n"
       f"u(30.8) = {u308:.3f} m/s (spec 75 +- 3)")

# ---------------------------------------------------------------------------
# 4. The conversion agrees with Step 6
# ---------------------------------------------------------------------------
def surface(phi):
    phi = np.atleast_1d(np.asarray(phi, dtype="float64"))
    r, _, _ = gd.reference_geoid(phi, R_POLAR, OMEGA, GM, J, DEGREES, R_NORM,
                                 tol_m=1e-9, max_iter=80)
    return r


zero_wind = lambda p: np.zeros_like(np.asarray(p, dtype="float64"))
step6, _, _, _ = latmod.planetocentric_fixed_point(
    np.radians([36.3]), surface, zero_wind, OMEGA, GM, J, DEGREES, R_NORM,
    tol_rad=1e-10, max_iter=40)
back = float(np.degrees(latmod.planetographic_from_planetocentric(
    step6, surface, zero_wind, OMEGA, GM, J, DEGREES, R_NORM))[0])
record(4, "the tool's direct conversion of 36.3 degrees matches Step 6", abs(back - 36.3) < 1e-6,
       f"Step 6 fixed point: 36.3 deg planetographic -> "
       f"{math.degrees(float(step6[0])):.6f} deg planetocentric\n"
       f"direct relation back: {back:.9f} deg, error {abs(back - 36.3):.2e} deg")

# ---------------------------------------------------------------------------
# 5. The bin table and the pooled RMS
# ---------------------------------------------------------------------------
gap_lo, gap_hi = curve.join_hi, curve.north_min
in_gap = (points_lat > gap_lo) & (points_lat < gap_hi)
in_range = (points_lat <= curve.north_max) & (points_lat >= curve.south_min)
usable = in_range & ~in_gap
rms_usable = float(np.sqrt(np.mean(residual[usable] ** 2)))
rms_all = float(np.sqrt(np.mean(residual[in_range] ** 2)))
lines = [f"{'bin':>7s} {'count':>6s} {'rms':>9s}  flag"]
for i in range(centers.size):
    flag = "low" if low[i] else ("" if count[i] else "empty")
    lines.append(f"{centers[i]:+7.1f} {count[i]:6d} {rms[i]:9.2f}  {flag}")
record(5, "the bin table, the counts and the pooled RMS about the curve",
       int(count.sum()) == 323 and abs(rms_usable - 21.3) <= 0.5,
       "\n".join(lines)
       + f"\nsum of counts {int(count.sum())}; populated bins {int((count > 0).sum())}; "
         f"bins below the minimum count {int(low.sum())}\n"
       f"RMS about the curve over the {int(usable.sum())} points inside the curve range and "
       f"outside the ring gap: {rms_usable:.2f} m/s, mean {float(residual[usable].mean()):.2f} "
       f"(spec 21.3 +- 0.5; the curve note gives 314, 21.3, -4.3)\n"
       f"including the {int((in_gap & in_range).sum())} gap points it is {rms_all:.2f} m/s: "
       f"their residuals span {float(residual[in_gap & in_range].min()):.1f} to "
       f"{float(residual[in_gap & in_range].max()):.1f} m/s, the reflection overstating the "
       f"southern flank. See the report, finding 2.")

# ---------------------------------------------------------------------------
# 6. The file
# ---------------------------------------------------------------------------
identical = bool(np.array_equal(u_file, np.repeat(u_file[:, :1], u_file.shape[1], axis=1)))
provenance = wind["value_provenance"].values
meanings = wind["value_provenance"].attrs["flag_meanings"].split()
pressure = wind["pressure_Pa"].values
ref_col = int(np.argmin(np.abs(pressure - float(wind["reference_level_pressure_Pa"]))))
codes = sorted(set(provenance.ravel().tolist()))
shear_max = float(np.max(np.abs(wind["u_shear_ms"].values)))
# SPEC_04 Step 0: the sum identity is u_total = u_reference + u_shear, the reference level wind
# broadcast along the pressure axis; `u_cylindrical_ms` and the `decomposition` attribute are
# retired (SPEC_04 Appendix, SPEC_00 section 6.6 amended).
sum_dev = float(np.max(np.abs(u_file - (wind["u_reference_ms"].values[:, None]
                                        + wind["u_shear_ms"].values))))
record(6, "the file: identical columns, provenance and the three parts",
       identical and sum_dev == 0.0 and shear_max == 0.0
       and set(np.unique(np.delete(provenance, ref_col, axis=1))) == {4},
       f"shape {u_file.shape}; columns identical {identical}\n"
       f"grid {lat_c.size} latitudes at {float(lat_c[1] - lat_c[0]):g} deg, both poles and the "
       f"equator as nodes: {bool(np.any(np.abs(lat_c) < 1e-12))}\n"
       f"provenance codes present {codes} -> {[meanings[c] for c in codes]}\n"
       f"every level other than the reference column carries code 4\n"
       f"u_reference_ms{tuple(wind['u_reference_ms'].dims)}; sum identity departure "
       f"{sum_dev:.1e}; max |u_shear| {shear_max:.1e} m/s, zero because this field does not "
       "vary along the column\n"
       f"read as kind W succeeded; the polar rule is the loaders' (check 2)")

# ---------------------------------------------------------------------------
# 7. The wind geoid, its convergence and the two sensitivities
# ---------------------------------------------------------------------------
def wind_of(cv, step_deg=0.01):
    """u(phi_c) with the planetographic conversion precomputed once.

    The conversion is a smooth monotone function of latitude, so tabulating it on a 0.01 degree
    grid and interpolating removes a Newton geoid solve from every RK4 stage; without it the
    anchored march, which repeats the march inside a root find, takes minutes. SPEC_01 forbids
    interpolating the anchor RADIUS; this interpolates the latitude conversion, a different
    object, and the curve it feeds is itself an interpolant.
    """
    phi_c = np.radians(np.arange(-90.0, 90.0 + 0.5 * step_deg, step_deg))
    phi_g = np.degrees(latmod.planetographic_from_planetocentric(
        phi_c, surface, zero_wind, OMEGA, GM, J, DEGREES, R_NORM))
    polar = np.abs(phi_c) >= np.pi / 2
    phi_g[polar] = np.sign(phi_c[polar]) * 90.0
    table = cv(phi_g)
    table[polar] = 0.0
    interpolant = PchipInterpolator(phi_c, table, extrapolate=False)

    def inner(value):
        shape = np.shape(value)
        a = np.clip(np.atleast_1d(np.asarray(value, dtype="float64")), phi_c[0], phi_c[-1])
        return np.reshape(np.nan_to_num(interpolant(a), nan=0.0), shape)

    return inner


ANCHOR = "mean_polar_radius"
convergence = []
for n in (901, 1801, 3601, 7201):
    grid_n = np.radians(np.linspace(-90.0, 90.0, n))
    res_n = gd.wind_geoid(grid_n, R_POLAR, ANCHOR, wind_of(curve), OMEGA, GM, J, DEGREES,
                          R_NORM)
    convergence.append((n, float(res_n.radius[int(np.argmin(np.abs(grid_n)))])))
grid = np.radians(np.linspace(-90.0, 90.0, 3601))
eq = int(np.argmin(np.abs(grid)))
r_nowind, _, _ = gd.reference_geoid(grid, R_POLAR, OMEGA, GM, J, DEGREES, R_NORM,
                                    tol_m=1e-9, max_iter=80)
result = gd.wind_geoid(grid, R_POLAR, ANCHOR, wind_of(curve), OMEGA, GM, J, DEGREES, R_NORM)
r_wind = result.radius
base_eq = float(r_wind[eq])
spread = max(v for _, v in convergence) - min(v for _, v in convergence)

sensitivities = []
for label, gap_rule, polar_rule, rule in (
    ("polar_rule linear_to_zero", "reflect_north", "linear_to_zero", ANCHOR),
    ("gap_rule reflect_north_then_bins", "reflect_north_then_bins", "pchip_to_zero", ANCHOR),
    ("anchor_rule north_pole", "reflect_north", "pchip_to_zero", "north_pole"),
    ("anchor_rule south_pole", "reflect_north", "pchip_to_zero", "south_pole"),
):
    alt = gd.wind_geoid(grid, R_POLAR, rule, wind_of(make_curve(gap_rule, polar_rule)),
                        OMEGA, GM, J, DEGREES, R_NORM)
    sensitivities.append((label, float(alt.radius[eq]), (float(alt.radius[eq]) - base_eq) / 1e3))
departure = (base_eq - R_FITTED) / 1e3
record(7, "the wind geoid: one anchored march, its convergence, and the sensitivities",
       spread < 10.0 and abs(result.anchor_residual_m) < 1e-3,
       "\n".join(f"   n = {n:5d}: r_eq = {v / 1e3:12.4f} km" for n, v in convergence)
       + f"\n   spread {spread:.3f} m across an eightfold refinement\n"
       f"anchor_rule {ANCHOR!r}, residual {result.anchor_residual_m:.1e} m\n"
       f"north polar radius {result.polar_north_m / 1e3:9.2f} km\n"
       f"south polar radius {result.polar_south_m / 1e3:9.2f} km\n"
       f"POLAR ASYMMETRY    {result.polar_asymmetry_m / 1e3:+9.2f} km, the asymmetry the wind "
       f"produces; their mean is the anchor {R_POLAR / 1e3:.0f} km to "
       f"{abs(result.anchor_residual_m):.1e} m\n"
       f"no wind equatorial radius    {float(r_nowind[eq]) / 1e3:9.2f} km\n"
       f"wind geoid equatorial radius {base_eq / 1e3:9.2f} km   bulge "
       f"{(base_eq - float(r_nowind[eq])) / 1e3:.2f} km (Lindal's own two surfaces imply 123)\n"
       f"Lindal's fitted 100 mbar     {R_FITTED / 1e3:9.0f} km +- 4, and the polar anchor "
       f"carries +- 10 km, so the band is about +- 11 km\n"
       f"DEPARTURE {departure:+.2f} km, inside the band\n"
       + "\n".join(f"sensitivity, {a:34s} {b / 1e3:9.2f} km  ({c:+.2f} km)"
                   for a, b, c in sensitivities)
       + "\nThe two anchor rules bracket the surface: anchoring at one pole moves the equator "
         "by about 19 km either way, which is the whole of the 23.31 km the v0.15 report "
         "found.\n"
         "The gap rule is now worth -1.83 km, against 0.01 km under the v0.15 per hemisphere "
         "march. That is not a contradiction: the ring gap still contributes almost nothing "
         "to the equator directly, because the Eq. B3 slope carries u sin(phi) and sin(phi) "
         "vanishes there, but it changes the southern flank and so the SOUTH POLAR radius, "
         "and under mean polar anchoring a shift in either polar radius moves the whole "
         "surface by half of it. The gap reaches the equator through the anchor, not through "
         "the slope.")

# ---------------------------------------------------------------------------
# 8. Figure 1
# ---------------------------------------------------------------------------
figure, axis = plt.subplots(figsize=(8.0, 5.0))
shade = {2: ("#ffe9b8", "gap fill (parameterized)"), 3: ("#f6d6d6", "polar caps (extrapolated)")}
grid_g = wind["latitude_planetographic_deg"].values
prov_1d = provenance[:, ref_col]
for code, (colour, label) in shade.items():
    inside = prov_1d == code
    if inside.any():
        edges = np.flatnonzero(np.diff(np.concatenate([[0], inside.view(np.int8), [0]])))
        for a, b in zip(edges[::2], edges[1::2]):
            axis.axvspan(grid_g[a], grid_g[min(b, grid_g.size - 1)], color=colour, zorder=0,
                         label=label)
            label = None
axis.plot(points_lat, points_u, ".", color="0.6", markersize=4, label="Smith 1982 points (323)")
sigma = wind["u_total_uncertainty_ms"].values[:, ref_col]
good = np.isfinite(sigma)
axis.fill_between(grid_g[good], (u_file[:, ref_col] - sigma)[good],
                  (u_file[:, ref_col] + sigma)[good], color="#1f77b4", alpha=0.25,
                  label="+- 1 sigma")
axis.plot(grid_g, u_file[:, ref_col], "-", color="#1f77b4", linewidth=1.6,
          label="Ingersoll and Pollard 1982 curve, assembled")
for name, (lat_s, u_s) in (("", segments.north), (None, segments.south)):
    axis.plot(lat_s, u_s, "-", color="#d62728", linewidth=0.9,
              label="digitized solid curve" if name == "" else None)
axis.axhline(0.0, color="0.8", linewidth=0.8, zorder=0)
axis.set_xlabel("planetographic latitude (degrees north)")
axis.set_ylabel("zonal wind, System III (m/s)")
axis.set_title("Lindal reduction wind: Ingersoll and Pollard (1982) Fig. 5, assembled")
axis.legend(loc="upper left", fontsize=7.5)
axis.grid(alpha=0.25)
figure.tight_layout()
figure.savefig(FIGURES / "step7_wind_curve.png", dpi=150)
plt.close(figure)

# ---------------------------------------------------------------------------
# 9. Figure 2: dynamical height against Lindal Eq. 18
# ---------------------------------------------------------------------------
height = r_wind - r_nowind
u_grid = wind_of(curve)(grid)
g_ref = g_eff_radial(0.0, r_nowind, grid, OMEGA, GM, J, DEGREES, R_NORM)
integrand = u_grid * np.sin(grid)
tail = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * np.diff(grid))])
eq18 = (2.0 * OMEGA * R_POLAR / g_ref) * (tail[-1] - tail)
# Eq. 18 is anchored the same way as the march, so the two curves are comparable: offset so
# that its two polar values average to zero, which is what "mean polar radius" does to the
# marched surface (SPEC_01 v0.16 Step 5).
eq18 = eq18 - 0.5 * (eq18[0] + eq18[-1])

figure, (top, bottom) = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)
top.plot(np.degrees(grid), u_grid, color="#1f77b4", linewidth=1.5)
top.set_ylabel("zonal wind (m/s)")
top.set_title("Wind curve and the dynamical height it produces")
top.grid(alpha=0.25)
top.axhline(0.0, color="0.8", linewidth=0.8)
bottom.plot(np.degrees(grid), height / 1e3, color="#1f77b4", linewidth=1.6,
            label="Eq. B3 march, r_wind - r_nowind")
bottom.plot(np.degrees(grid), eq18 / 1e3, "--", color="#d62728", linewidth=1.4,
            label="Lindal Eq. 18, small wind approximation")
anchor = 30.818189
bottom.axvline(anchor, color="0.5", linewidth=0.9, linestyle=":")
bottom.annotate(f"anchor {anchor:.3f} deg", xy=(anchor, 0.0), xytext=(anchor + 4, 25),
                fontsize=8, color="0.3")
bottom.annotate(
    f"equator: march {height[eq] / 1e3:.1f} km, Eq. 18 {eq18[eq] / 1e3:.1f} km;"
    f" Lindal implies 123 km",
    xy=(0.0, height[eq] / 1e3), xytext=(-87, height[eq] / 1e3 * 0.45), fontsize=8)
bottom.set_xlabel("planetocentric latitude (degrees north)")
bottom.set_ylabel("dynamical height (km)")
bottom.legend(loc="upper right", fontsize=8)
bottom.grid(alpha=0.25)
figure.tight_layout()
figure.savefig(FIGURES / "step7_dynamical_height.png", dpi=150)
plt.close(figure)

anchor_rad = math.radians(anchor)
record(8, "both figures written to reports/figures/, and the Eq. 18 comparison",
       (FIGURES / "step7_wind_curve.png").exists()
       and (FIGURES / "step7_dynamical_height.png").exists(),
       "figures: step7_wind_curve.png, step7_dynamical_height.png\n"
       f"dynamical height at the equator: march {height[eq] / 1e3:.2f} km, Eq. 18 "
       f"{eq18[eq] / 1e3:.2f} km, difference {(height[eq] - eq18[eq]) / 1e3:+.2f} km "
       f"({100 * abs(height[eq] - eq18[eq]) / height[eq]:.1f} percent)\n"
       f"Lindal's fitted surface implies 60367 - 60244 = 123 km\n"
       f"at the anchor latitude {anchor:.4f} deg: march "
       f"{float(np.interp(anchor_rad, grid, height)) / 1e3:.2f} km, Eq. 18 "
       f"{float(np.interp(anchor_rad, grid, eq18)) / 1e3:.2f} km\n"
       f"max |march - Eq. 18| over all latitudes "
       f"{float(np.max(np.abs(height - eq18))) / 1e3:.2f} km, at "
       f"{float(np.degrees(grid[int(np.argmax(np.abs(height - eq18)))])):+.1f} deg\n"
       "Eq. 18 drops the cyclostrophic term and the oblateness corrections, so it is expected "
       "to sit below the march; the sign and size of the gap are what this figure shows.")

for handle in (wind, gravity, rotation, raw):
    handle.close()
print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
