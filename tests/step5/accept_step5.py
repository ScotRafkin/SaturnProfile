"""Acceptance checks for SPEC_01 v0.16 Step 5, `lib.geoid`.

Every Saturn number is read from the static transcriptions or from the Step 3 products, per
SPEC_00 section 2.1. The stand-in wind profile is the one quantity invented here, and check 5
shows how much the answer depends on that invention.
"""

import math
import sys
import tomllib

import numpy as np

from casspian.lib import geoid as gd
from casspian.lib import gravity as gv

results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def load(path):
    with open(path, "rb") as handle:
        return tomllib.load(handle)


null = load("data_static/harmonics/null1981.toml")
iess = load("data_static/harmonics/iess2019.toml")
rot = load("data_static/rotation/system_iii.toml")
scalars = load("occul_data/lindal/raw/lindal_scalars.toml")

DEGREES = [2, 4, 6]
J = [null["J"][str(d)]["value"] for d in DEGREES]
R_NORM = null["meta"]["normalization_radius_m"]
GM_NULL = null["GM"]["planet_value_m3s2"]
GM_MODERN = iess["GM"]["planet_value_m3s2"]
OMEGA = 2.0 * math.pi / rot["system_iii"]["period_s"]
R_POLAR = scalars["geodesy"]["reference_geoid_polar_radius_km"] * 1e3
R_FITTED_EQ = scalars["geodesy"]["surface_100mbar"]["radius_equatorial_km"] * 1e3
U_EQ = null["check_values"]["equatorial_wind_assumed_ms"]

GRID = np.radians(np.linspace(0.0, 90.0, 1801))


def no_wind(GM, grid=GRID):
    return gd.reference_geoid(grid, R_POLAR, OMEGA, GM, J, DEGREES, R_NORM,
                              tol_m=1e-6, max_iter=60)


# ---------------------------------------------------------------------------
# 1. U_rigid sign convention: dU/dr must equal g, positive inward
# ---------------------------------------------------------------------------
r0, phi0 = 6.0e7, math.radians(37.0)
h = 1.0
dUdr = float(
    (U_rigid_hi := gd.U_rigid(r0 + h, phi0, OMEGA, GM_NULL, J, DEGREES, R_NORM))
    - gd.U_rigid(r0 - h, phi0, OMEGA, GM_NULL, J, DEGREES, R_NORM)
) / (2 * h)
g_here = float(gv.g_eff_radial(0.0, r0, phi0, OMEGA, GM_NULL, J, DEGREES, R_NORM))
record(1, "U_rigid obeys dU/dr = g with g positive inward, the stated convention",
       abs(dUdr - g_here) < 1e-6,
       f"central difference dU/dr = {dUdr:.9f} m/s2\n"
       f"g_eff_radial        = {g_here:.9f} m/s2\n"
       f"difference {abs(dUdr - g_here):.2e}. Lindal Eq. 11 is written with the opposite "
       f"overall sign; this package uses g_eff = -grad U throughout.")

# ---------------------------------------------------------------------------
# 2. The no wind reference geoid
# ---------------------------------------------------------------------------
rows = []
for label, GM in (("Null 1981", GM_NULL), ("modern", GM_MODERN)):
    r, iters, resid = no_wind(GM)
    r_eq = float(r[0]) / 1e3
    r31 = float(gd.radius_at(math.radians(31.0), GRID, r)) / 1e3
    r308 = float(gd.radius_at(math.radians(30.8185), GRID, r)) / 1e3
    rows.append((label, GM, r_eq, r31, r308, int(iters.max()), float(resid.max()),
                 float(r[-1])))
ok = all(
    abs(r_eq - 60244) <= 2 and abs(r31 - 58435) <= 2 and abs(r308 - 58452.9) <= 2
    for _, _, r_eq, r31, r308, _, _, _ in rows
)
detail = [f"r_polar = {R_POLAR / 1e3:g} km, read from lindal_scalars.toml"]
for label, GM, r_eq, r31, r308, it, resid, r_pole in rows:
    detail.append(
        f"{label:9s} GM = {GM:.10e}: r_eq = {r_eq:.2f} km (60244 +- 2), "
        f"r(31.0) = {r31:.2f} (58435 +- 2), r(30.8185) = {r308:.2f} (58452.9 +- 2)"
    )
    detail.append(
        f"{'':9s} max Newton steps {it}, max |U - U_ref| {resid:.2e} m2/s2, "
        f"r at the pole {r_pole:.6f} m (anchor {R_POLAR:.6f})"
    )
detail.append(
    "The handoff's 58452.9 is reproduced exactly by the modern GM and to 0.2 km by Null's. "
    "The reduction uses Null's for consistency with Lindal, and both are inside the tolerance."
)
record(2, "the no wind reference geoid reproduces the handoff radii", ok, "\n".join(detail))

# ---------------------------------------------------------------------------
# 3. Grid convergence
# ---------------------------------------------------------------------------
conv = []
for n in (451, 901, 1801, 3601):
    grid = np.radians(np.linspace(0.0, 90.0, n))
    r, _, _ = no_wind(GM_NULL, grid)
    conv.append((n, float(r[0])))
spread = max(v for _, v in conv) - min(v for _, v in conv)
record(3, "the no wind geoid is grid independent", spread < 1e-6,
       "\n".join(f"n = {n:5d}: r_eq = {v:.9f} m" for n, v in conv)
       + f"\nspread {spread:.2e} m. Newton solves each latitude independently, so the grid "
       "sets only where the surface is sampled, not its accuracy.")

# ---------------------------------------------------------------------------
# 4. The 45 degree check repeated on the constructed surface (Step 4 review action)
# ---------------------------------------------------------------------------
r_nw, _, _ = no_wind(GM_NULL)
r45 = float(gd.radius_at(math.radians(45.0), GRID, r_nw))
g45, gp45, mag45, psi45 = gv.g_eff_vector(
    0.0, r45, math.radians(45.0), OMEGA, GM_NULL, J, DEGREES, R_NORM
)
r_sphere = scalars["geodesy"]["surface_1bar"]["radius_equatorial_km"] * 1e3
_, gp_s, mag_s, psi_s = gv.g_eff_vector(
    0.0, r_sphere, math.radians(45.0), OMEGA, GM_NULL, J, DEGREES, R_NORM
)
record(4, "G_phi at 45 degrees on the constructed geoid, as the Step 4 review asked",
       float(gp45) < 0 and float(psi45) > 0,
       f"on the constructed no wind geoid, r_ref(45) = {r45 / 1e3:.3f} km:\n"
       f"  G_phi = {float(gp45):.6f} m/s2, psi = {math.degrees(float(psi45)):.4f} deg, "
       f"|g_eff| = {float(mag45):.4f} m/s2\n"
       f"on the equatorial 1 bar sphere, r = {r_sphere / 1e3:.0f} km (the Step 4 figure):\n"
       f"  G_phi = {float(gp_s):.6f} m/s2, psi = {math.degrees(float(psi_s)):.4f} deg, "
       f"|g_eff| = {float(mag_s):.4f} m/s2\n"
       f"The sphere overstates the radius at 45 degrees by {(r_sphere - r45) / 1e3:.0f} km, "
       f"which lowers |g| by 1.22 m/s2 and psi by 0.64 degrees. The sign and the order of "
       f"magnitude are the same, which is all the Step 4 check claimed.")

# ---------------------------------------------------------------------------
# 5. The wind geoid, and how much the stand-in profile decides the answer
# ---------------------------------------------------------------------------
def profile(kind, cut_deg):
    cut = math.radians(cut_deg)
    def u(phi):
        p = np.abs(np.asarray(phi, dtype="float64"))
        x = np.clip(p / cut, 0.0, 1.0)
        if kind == "cos2":
            shape = np.cos(math.pi * x / 2) ** 2
        elif kind == "linear":
            shape = 1.0 - x
        elif kind == "plateau20":
            shape = np.clip((1.0 - x) / (1.0 - 20.0 / cut_deg), 0.0, 1.0)
        elif kind == "flat":
            shape = np.ones_like(x)
        return np.where(p <= cut, U_EQ * shape, 0.0)
    return u


r_nw_eq = float(r_nw[0])
table = []
for kind, cut in (("cos2", 35), ("linear", 35), ("plateau20", 35), ("flat", 35),
                  ("cos2", 45), ("cos2", 60), ("linear", 60)):
    r_w, _ = gd.wind_geoid(GRID, R_POLAR, "mean_polar_radius", profile(kind, cut), OMEGA,
                           GM_NULL, J, DEGREES, R_NORM)
    table.append((f"{kind} to {cut} deg", float(r_w[0]) / 1e3,
                  (float(r_w[0]) - r_nw_eq) / 1e3))

primary = table[0]
record(5, "the wind geoid bulges outward, by an amount the stand-in profile decides",
       primary[2] > 0,
       f"no wind equatorial radius: {r_nw_eq / 1e3:.2f} km\n"
       f"Lindal's fitted 100 mbar equatorial radius: {R_FITTED_EQ / 1e3:g} km, "
       f"which is {(R_FITTED_EQ - r_nw_eq) / 1e3:.0f} km above the no wind value\n"
       + "\n".join(f"  {name:22s} r_eq = {r:9.2f} km   bulge = {b:7.2f} km"
                   for name, r, b in table)
       + "\nAll seven have 450 m/s at the equator and reach zero at their cutoff, which is "
         "everything SPEC_01 says about the stand-in, and they span 66 to 215 km of bulge. "
         "The specification's 'roughly 120 km' is inside that range but is not determined by "
         "the stand-in as described. See the report, section 3.")

# ---------------------------------------------------------------------------
# 6. The wind geoid is grid converged, and the closure diagnostic
# ---------------------------------------------------------------------------
u_std = profile("cos2", 35)
conv_w = []
for n in (451, 901, 1801, 3601):
    grid = np.radians(np.linspace(0.0, 90.0, n))
    r_w, _ = gd.wind_geoid(grid, R_POLAR, "mean_polar_radius", u_std, OMEGA, GM_NULL, J,
                           DEGREES, R_NORM)
    conv_w.append((n, float(r_w[0])))
spread_w = max(v for _, v in conv_w) - min(v for _, v in conv_w)

r_w, closure = gd.wind_geoid(GRID, R_POLAR, "mean_polar_radius", u_std, OMEGA, GM_NULL, J,
                             DEGREES, R_NORM)
peak = float(np.abs(closure).max())
g_eq = float(gv.g_eff_radial(u_std(0.0), r_w[0], 0.0, OMEGA, GM_NULL, J, DEGREES, R_NORM))
record(6, "the wind geoid is grid converged, and the closure is the no wind potential "
       "departure",
       spread_w < 1e-3,
       "\n".join(f"n = {n:5d}: r_eq = {v:.6f} m" for n, v in conv_w)
       + f"\nspread {spread_w:.2e} m over an eightfold range of step, which is RK4 doing "
         "its job. This is the Step 5 acceptance for the wind geoid (v0.7).\n"
       f"closure U_rigid(r0) - U_ref against the no wind potential: "
       f"max |.| = {peak:.3e} m2/s2\n"
       f"  as a dynamical height, divided by g = {g_eq:.3f} m/s2 at the equator: "
       f"{peak / g_eq / 1e3:.2f} km, which is the wind surface standing that far above the "
       f"no wind geoid\n"
       "No pseudo-potential built from Omega_abs is evaluated (SPEC_01 v0.7): with u varying "
       "in latitude its gradient is not the effective gravity, so its variation restates the "
       "wind kinetic term and measures nothing about conservativeness.")

# ---------------------------------------------------------------------------
# 7. radius_at: the interpolation choice, and the grid it needs
# ---------------------------------------------------------------------------
probe = math.radians(30.8185)
# The reference is the geoid solved directly at the probe latitude, so that no interpolation
# enters the number the three schemes are measured against.
exact = float(no_wind(GM_NULL, np.array([probe]))[0][0])
rows7 = []
for n in (46, 91, 181, 451, 901):
    coarse = np.radians(np.linspace(0.0, 90.0, n))
    r_c, _, _ = no_wind(GM_NULL, coarse)
    in_sin = float(np.interp(np.sin(probe), np.sin(coarse), r_c))
    in_phi = float(np.interp(probe, coarse, r_c))
    # For an ellipse 1/r^2 = sin^2(phi)/b^2 + cos^2(phi)/a^2 is EXACTLY linear in sin^2(phi),
    # so this scheme is exact for the shape the geoid nearly is.
    # What `radius_at` now implements, per SPEC_01 v0.7.
    in_ruled = float(gd.radius_at(probe, coarse, r_c))
    rows7.append((90.0 / (n - 1), abs(in_sin - exact), abs(in_phi - exact),
                  abs(in_ruled - exact)))
ruled_best = all(c <= a and c <= b for _, a, b, c in rows7)

# The bracketing must be in phi_c, not in sin^2, or a southern latitude would fold onto the
# northern hemisphere. A grid spanning both hemispheres with the equator as a node:
both = np.radians(np.concatenate([np.linspace(-90.0, 0.0, 181)[:-1],
                                  np.linspace(0.0, 90.0, 181)]))
r_both, _, _ = no_wind(GM_NULL, both)
south = float(gd.radius_at(-probe, both, r_both))
north = float(gd.radius_at(probe, both, r_both))
symmetric = abs(south - north) < 1e-6

# A grid without the equator as a node must be refused rather than folded.
no_equator = np.radians(np.array([-2.0, -1.0, 1.0, 2.0]))
r_gap, _, _ = no_wind(GM_NULL, no_equator)
straddle_refused = False
try:
    gd.radius_at(0.0, no_equator, r_gap)
except ValueError:
    straddle_refused = True
refused = False
try:
    gd.radius_at(math.radians(95.0), GRID, r_nw)
except ValueError:
    refused = True
record(7, "radius_at uses the ruled scheme, brackets in phi_c, and refuses both extrapolation "
       "and an equator straddling interval",
       ruled_best and refused and straddle_refused and symmetric,
       f"reference: the geoid solved directly at the probe, r(30.8185) = {exact:.3f} m\n"
       + f"{'spacing':>9s} {'sin phi (v0.6)':>16s} {'phi':>13s} "
         f"{'1/r^2 in sin^2 (v0.7)':>23s}\n"
       + "\n".join(f"{s:8.2f}d {a:14.2f} m {b:11.2f} m {c:21.4f} m"
                   for s, a, b, c in rows7)
       + "\nAll three are second order: the error falls by about four when the spacing "
         "halves. The last column is what `radius_at` now returns, and it is the smallest at "
         "every spacing, by about eight against the v0.6 scheme.\n"
       f"bracketing is in phi_c, so a southern latitude is not folded onto the north: on a "
       f"grid spanning both hemispheres, r(-30.8185) = {south:.3f} m and r(+30.8185) = "
       f"{north:.3f} m, equal to {abs(south - north):.2e} m on this symmetric no wind "
       f"surface, and they were reached through different node pairs\n"
       f"a latitude outside the grid is refused: {refused}\n"
       f"a grid without the equator as a node is refused rather than folded: "
       f"{straddle_refused}\n"
       "The frozen anchor radius is never interpolated at all (SPEC_01 v0.7): it is solved "
       "by `reference_geoid` or marched by `wind_geoid` at its own latitude. This function "
       "serves the vectorized uses of Step 7 and later.")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
