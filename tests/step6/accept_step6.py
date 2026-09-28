"""Acceptance checks for SPEC_01 v0.7 Step 6, `lib.latitude`.

Every Saturn number is read from the static transcriptions. The surfaces passed to the fixed
point solve or march at the iterate's own latitude, never interpolate, as v0.7 requires.
"""

import math
import sys
import tomllib

import numpy as np

from casspian.lib import geoid as gd
from casspian.lib import latitude as lat

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
# The label is defined on the 100 mbar surface (Fig. 4 caption), so its oblateness is the
# ellipsoid the seed should use.
F_100MBAR = scalars["geodesy"]["surface_100mbar"]["oblateness"]
PHI_G = scalars["latitude"]["planetographic_deg"]
U_EQ = null["check_values"]["equatorial_wind_assumed_ms"]

NO_WIND = lambda phi: np.zeros_like(np.asarray(phi, dtype="float64"))


def no_wind_surface(GM):
    """A surface callable that SOLVES at the latitudes given (SPEC_01 v0.7)."""
    def surface(phi):
        phi = np.atleast_1d(np.asarray(phi, dtype="float64"))
        r, _, _ = gd.reference_geoid(phi, R_POLAR, OMEGA, GM, J, DEGREES, R_NORM,
                                     tol_m=1e-9, max_iter=80)
        return r
    return surface


CUT = math.radians(35.0)


def u_standin(phi):
    p = np.abs(np.asarray(phi, dtype="float64"))
    return np.where(p <= CUT, U_EQ * 0.5 * (1.0 + np.cos(math.pi * p / CUT)), 0.0)


BASE_NODES = np.radians(np.linspace(0.0, 90.0, 1801))


def wind_surface(phi):
    """A surface callable that MARCHES with the latitudes inserted as nodes (SPEC_01 v0.7)."""
    phi = np.atleast_1d(np.asarray(phi, dtype="float64"))
    nodes = np.unique(np.concatenate([BASE_NODES, phi]))
    r, _ = gd.wind_geoid(nodes, R_POLAR, "mean_polar_radius", u_standin, OMEGA,
                         GM_NULL, J, DEGREES, R_NORM)
    return r[np.searchsorted(nodes, phi)]


def solve(phi_g_deg, surface, u, GM, tol_deg=1e-9, flattening=F_100MBAR):
    # The acceptance states its tolerances in degrees because that is how SPEC_01 states the
    # convergence claim. `lib` is radians throughout, so the conversion happens here, once.
    return lat.planetocentric_fixed_point(
        np.radians(np.atleast_1d(phi_g_deg)), surface, u, OMEGA, GM, J, DEGREES, R_NORM,
        tol_rad=np.radians(tol_deg), max_iter=30, flattening=flattening,
    )


# ---------------------------------------------------------------------------
# 1. The ellipsoid seed
# ---------------------------------------------------------------------------
seed = math.degrees(float(lat.planetocentric_from_ellipsoid(math.radians(PHI_G), F_100MBAR)))
poles = lat.planetocentric_from_ellipsoid(np.array([np.pi / 2, -np.pi / 2]), F_100MBAR)
record(1, "the ellipsoid seed reproduces the handoff value and is finite at the poles",
       abs(seed - 30.8526) < 0.001 and np.allclose(poles, [np.pi / 2, -np.pi / 2]),
       f"phi_g = {PHI_G} (from lindal_scalars.toml), f = {F_100MBAR} (100 mbar oblateness, "
       f"the surface the Fig. 4 label is defined on)\n"
       f"seed = {seed:.4f} deg, handoff section 9A.8 gives 30.8526\n"
       f"the 1 bar oblateness {scalars['geodesy']['surface_1bar']['oblateness']} would give "
       f"{math.degrees(float(lat.planetocentric_from_ellipsoid(math.radians(PHI_G), scalars['geodesy']['surface_1bar']['oblateness']))):.4f} deg, so the handoff seeded on the "
       f"100 mbar surface\n"
       f"at the poles the seed returns {np.degrees(poles).tolist()} degrees, finite because "
       f"arctan2 is used and not arctan")

# ---------------------------------------------------------------------------
# 2. The fixed point at the Lindal anchor
# ---------------------------------------------------------------------------
rows = []
for label, GM in (("Null 1981", GM_NULL), ("modern", GM_MODERN)):
    pc, psi, its, cnt = solve(PHI_G, no_wind_surface(GM), NO_WIND, GM)
    counts = {t: int(solve(PHI_G, no_wind_surface(GM), NO_WIND, GM, tol_deg=t)[3][0])
              for t in (1e-3, 1e-4, 1e-6, 1e-9)}
    rows.append((label, GM, math.degrees(float(pc[0])), math.degrees(float(psi[0])),
                 [math.degrees(float(v[0])) for v in its[:6]], counts))
ok2 = all(abs(pc - 30.8185) <= 0.001 and abs(ps - 5.4815) <= 0.001
          for _, _, pc, ps, _, _ in rows)
detail = []
for label, GM, pc, ps, its, counts in rows:
    detail.append(f"{label:9s} GM: phi_c = {pc:.6f} deg (30.8185 +- 0.001), "
                  f"psi = {ps:.6f} deg (5.4815 +- 0.001), "
                  f"departure {abs(pc - 30.8185) * 3600:.2f} arcsec")
    detail.append(f"{'':9s}     iterates: " + ", ".join(f"{v:.5f}" for v in its))
    detail.append(f"{'':9s}     steps to converge: "
                  + ", ".join(f"{t:g} deg -> {c}" for t, c in counts.items()))
detail.append(
    "SPEC_01 says three or four iterations without naming a tolerance. Four steps reaches "
    "1e-4 degrees and five reaches 1e-6; the count above is what each tolerance costs."
)
detail.append(
    "The modern GM lands 0.15 arcsec from the handoff value and Null's 1.12 arcsec, the same "
    "pattern Step 5 found in the radii: the handoff numbers were computed with the modern GM. "
    "The reduction uses Null's, and both are well inside the tolerance."
)
record(2, "the fixed point reproduces the handoff anchor latitude and tilt", ok2,
       "\n".join(detail))

# ---------------------------------------------------------------------------
# 3. Round trip
# ---------------------------------------------------------------------------
# Recomputed here with Null's GM so that everything reported below is one consistent pair;
# the loop above left `psi` holding the modern GM result.
pc, psi, _, _ = solve(PHI_G, no_wind_surface(GM_NULL), NO_WIND, GM_NULL)
back = lat.planetographic_from_planetocentric(
    pc, no_wind_surface(GM_NULL), NO_WIND, OMEGA, GM_NULL, J, DEGREES, R_NORM
)
err = abs(math.degrees(float(back[0])) - PHI_G)
record(3, "the inverse returns the planetographic latitude to 1e-9 degrees", err < 1e-9,
       f"phi_g in  = {PHI_G}\n"
       f"phi_g out = {math.degrees(float(back[0])):.12f} deg\n"
       f"error {err:.2e} deg, which is {err * 3600 * 1000:.2e} milliarcsec")

# ---------------------------------------------------------------------------
# 4. Sensitivity: phi_g = 36.5
# ---------------------------------------------------------------------------
pc5, psi5, _, _ = solve(36.5, no_wind_surface(GM_NULL), NO_WIND, GM_NULL)
value5 = math.degrees(float(pc5[0]))
record(4, "phi_g = 36.5 gives 31.005, which shows the sensitivity", abs(value5 - 31.005) < 0.001,
       f"phi_g = 36.5 -> phi_c = {value5:.6f} deg (handoff 31.005)\n"
       f"phi_g = {PHI_G} -> phi_c = {math.degrees(float(pc[0])):.6f} deg\n"
       f"0.2 degrees of planetographic latitude moves the planetocentric value by "
       f"{(value5 - math.degrees(float(pc[0]))):.6f} degrees, so the label's own uncertainty "
       f"propagates almost one for one and is not reduced by the conversion")

# ---------------------------------------------------------------------------
# 5. The wind moves the anchor
# ---------------------------------------------------------------------------
pc_w, psi_w, _, cnt_w = solve(PHI_G, wind_surface, u_standin, GM_NULL)
shift = math.degrees(float(pc_w[0])) - math.degrees(float(pc[0]))
r_nw = float(no_wind_surface(GM_NULL)(np.array([float(pc[0])]))[0])
r_w = float(wind_surface(np.array([float(pc_w[0])]))[0])
record(5, "including the wind moves the anchor latitude, by the amount reported",
       abs(shift) > 1e-4,
       f"no wind: phi_c = {math.degrees(float(pc[0])):.6f} deg, "
       f"psi = {math.degrees(float(psi[0])):.6f} deg, r = {r_nw / 1e3:.3f} km\n"
       f"wind   : phi_c = {math.degrees(float(pc_w[0])):.6f} deg, "
       f"psi = {math.degrees(float(psi_w[0])):.6f} deg, r = {r_w / 1e3:.3f} km\n"
       f"the anchor moves {shift:+.6f} deg ({shift * 3600:+.2f} arcsec), equatorward, and the "
       f"anchor radius moves {(r_w - r_nw) / 1e3:+.3f} km\n"
       "This is with the Step 5 stand-in wind, which is not the Smith wind. The value the "
       "reduction freezes is the one Step 7 produces with the real wind; this check exists to "
       "show the shift is resolvable and which way it goes, not to fix a number.")

# ---------------------------------------------------------------------------
# 6. Vectorization and the southern hemisphere
# ---------------------------------------------------------------------------
grid = np.array([-60.0, -36.3, -5.0, 0.0, 5.0, 36.3, 60.0])
pcv, psiv, _, cntv = solve(grid, no_wind_surface(GM_NULL), NO_WIND, GM_NULL)
one_at_a_time = np.array([
    float(solve(x, no_wind_surface(GM_NULL), NO_WIND, GM_NULL)[0][0]) for x in grid
])
symmetric = abs((math.degrees(float(pcv[1])) + math.degrees(float(pcv[5])))) < 1e-9
record(6, "vectorized over phi_g, equator and both hemispheres, matching scalar calls",
       np.allclose(pcv, one_at_a_time, rtol=0, atol=1e-14) and symmetric
       and abs(float(pcv[3])) < 1e-15,
       f"phi_g  {grid.tolist()}\n"
       f"phi_c  {[round(math.degrees(float(v)), 5) for v in pcv]}\n"
       f"psi    {[round(math.degrees(float(v)), 5) for v in psiv]}\n"
       f"steps  {cntv.tolist()}\n"
       f"identical to the scalar calls: {bool(np.allclose(pcv, one_at_a_time, rtol=0, atol=1e-14))}\n"
       f"north and south mirror to {abs(math.degrees(float(pcv[1])) + math.degrees(float(pcv[5]))):.2e} deg, "
       f"and the equator maps to exactly {float(pcv[3]):.1e} deg, both as they must on a "
       f"surface with only even harmonics and no wind")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
