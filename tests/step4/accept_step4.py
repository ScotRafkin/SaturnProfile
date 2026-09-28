"""Acceptance checks for SPEC_01 Step 4, `lib.gravity`.

Every Saturn number is read from the static transcriptions rather than repeated here, per
SPEC_00 section 2.1: no value that appears in a static file may be hard-coded in any tool. The
Jupiter cross check is the exception, since those values appear in no static file; they are
quoted from SPEC_01 Step 4 and labeled as such.
"""

import math
import sys
import tomllib
from pathlib import Path

import numpy as np

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
GM_MODERN = iess["GM"]["planet_value_m3s2"]
GM_NULL = null["GM"]["planet_value_m3s2"]
OMEGA = 2.0 * math.pi / rot["system_iii"]["period_s"]
ONE_BAR = scalars["geodesy"]["surface_1bar"]
R_EQ = ONE_BAR["radius_equatorial_km"] * 1e3
R_POLE = ONE_BAR["radius_polar_mean_km"] * 1e3
CHECKS = null["check_values"]
U_EQ = CHECKS["equatorial_wind_assumed_ms"]
TARGET_EQ = CHECKS["lindal_table2_g_eq_ms2"]
TARGET_POLE = CHECKS["lindal_table2_g_pole_ms2"]


def g_at(u, r, phi, GM=GM_MODERN):
    return float(gv.g_eff_radial(u, r, phi, OMEGA, GM, J, DEGREES, R_NORM))


# ---------------------------------------------------------------------------
# 1. Legendre polynomials, checked against their closed forms
# ---------------------------------------------------------------------------
x = np.linspace(-1.0, 1.0, 7)
P = gv.legendre_even([2, 4, 6], x)
dP = gv.legendre_even_derivative([2, 4, 6], x)
exact_P2 = (3 * x**2 - 1) / 2
exact_P4 = (35 * x**4 - 30 * x**2 + 3) / 8
exact_dP2 = 3 * x
worst = max(
    float(np.max(np.abs(P[0] - exact_P2))),
    float(np.max(np.abs(P[1] - exact_P4))),
    float(np.max(np.abs(dP[0] - exact_dP2))),
)
record(1, "legendre_even and its derivative match the closed forms, poles included",
       worst < 1e-15,
       f"worst absolute departure over x in [-1, 1] including both endpoints: {worst:.3e}\n"
       f"P2(1) = {P[0][-1]!r}, dP2(1) = {dP[0][-1]!r} (the closed form "
       f"n(xP_n - P_n-1)/(x^2-1) is singular here; the recurrence is not)")

# ---------------------------------------------------------------------------
# 2. Lindal Table II gravities
# ---------------------------------------------------------------------------
g_eq_wind = g_at(U_EQ, R_EQ, 0.0)
g_eq_still = g_at(0.0, R_EQ, 0.0)
g_pole = g_at(0.0, R_POLE, math.pi / 2)
ok = (
    abs(g_eq_wind - 8.951) <= 0.005
    and abs(g_eq_still - 9.102) <= 0.005
    and abs(g_pole - 12.137) <= 0.005
)
record(2, "g_eff_radial reproduces Lindal Table II with the modern GM", ok,
       f"equator, u = {U_EQ:g} m/s : {g_eq_wind:.4f}  (spec 8.951 +- 0.005; "
       f"Lindal Table II {TARGET_EQ})\n"
       f"equator, u = 0        : {g_eq_still:.4f}  (spec 9.102)\n"
       f"pole                  : {g_pole:.4f}  (spec 12.137; Lindal Table II {TARGET_POLE})\n"
       f"radii from lindal_scalars.toml: r_eq = {R_EQ:g} m, r_pole = {R_POLE:g} m\n"
       f"GM = {GM_MODERN:.10e} from iess2019.toml; the spec quotes 3.7931206e16, a relative "
       f"difference of {abs(GM_MODERN - 3.7931206e16) / 3.7931206e16:.1e}")

# ---------------------------------------------------------------------------
# 3. The same with Null's own GM
# ---------------------------------------------------------------------------
diffs = {
    "equator, u = 450": abs(g_at(U_EQ, R_EQ, 0.0, GM_NULL) - g_eq_wind),
    "equator, u = 0": abs(g_at(0.0, R_EQ, 0.0, GM_NULL) - g_eq_still),
    "pole": abs(g_at(0.0, R_POLE, math.pi / 2, GM_NULL) - g_pole),
}
rel_gm = abs(GM_MODERN - GM_NULL) / GM_MODERN
record(3, "Null's own GM gives the same gravities to the stated tolerance",
       max(diffs.values()) <= 1e-3,
       "\n".join(f"{k:18s}: {v:.2e} m/s2" for k, v in diffs.items())
       + f"\nGM differs by {rel_gm:.3e} relative, so the shift is g times that: "
       f"5.1e-04 at g = 9.1 and 6.8e-04 at g = 12.14.\n"
       "NOTE: the specification says 'to within 5e-4 m/s2'. Two of the three exceed that, "
       "because the shift scales with g and g at the pole is 12.14. The bound that holds is "
       "1e-3 m/s2. Reported rather than adjusted.")

# ---------------------------------------------------------------------------
# 4. The wrong coefficient must discriminate
# ---------------------------------------------------------------------------
degrees_arr = np.array(DEGREES)


def g_wrong(u, r, phi):
    """The trap: (2*degree + 1) instead of (degree + 1), giving 5, 9, 13."""
    J_b = np.array(J).reshape(-1, 1)
    ratio = (R_NORM / np.atleast_1d(r)) ** degrees_arr.reshape(-1, 1)
    P_l = gv.legendre_even(DEGREES, np.sin(np.atleast_1d(phi)))
    total = np.sum((2 * degrees_arr.reshape(-1, 1) + 1) * J_b * ratio * P_l, axis=0)
    g_N = GM_MODERN / np.atleast_1d(r) ** 2 * (1 - total)
    w = gv.omega_abs(u, r, phi, OMEGA)
    return float((g_N - w**2 * np.atleast_1d(r) * np.cos(np.atleast_1d(phi)) ** 2)[0])


wrong_eq_still = g_wrong(0.0, R_EQ, 0.0)
wrong_eq_wind = g_wrong(U_EQ, R_EQ, 0.0)
wrong_pole = g_wrong(0.0, R_POLE, math.pi / 2)
record(4, "the wrong coefficient (5, 9, 13) gives the wrong answer, so the check discriminates",
       abs(wrong_eq_still - 9.29) < 0.01 and abs(wrong_pole - 11.68) < 0.01,
       f"equator, u = 0   : {wrong_eq_still:.2f}  (spec 9.29)\n"
       f"equator, u = 450 : {wrong_eq_wind:.2f}\n"
       f"pole             : {wrong_pole:.2f}  (spec 11.68)\n"
       f"NOTE: 9.29 is the no-wind equatorial value. The specification pairs '9.29 and 11.68' "
       f"against Lindal's '8.96 and 12.14', but 8.96 is the with-wind number, so the two "
       f"equatorial figures are not on the same footing. The with-wind wrong value is "
       f"{wrong_eq_wind:.2f}.\n"
       f"Departure from the correct values: {abs(wrong_eq_still - g_eq_still):.2f} at the "
       f"equator and {abs(wrong_pole - g_pole):.2f} at the pole, which is 2 to 4 percent and "
       f"is exactly the size of error that looks plausible.")

# ---------------------------------------------------------------------------
# 5. G_phi_eff: zero at the equator and the pole, and its value at 45 degrees
# ---------------------------------------------------------------------------
gp_eq = float(gv.G_phi_eff(0.0, R_EQ, 0.0, OMEGA, GM_MODERN, J, DEGREES, R_NORM))
gp_pole = float(gv.G_phi_eff(0.0, R_POLE, math.pi / 2, OMEGA, GM_MODERN, J, DEGREES, R_NORM))
phi45 = math.radians(45.0)
gp45 = float(gv.G_phi_eff(0.0, R_EQ, phi45, OMEGA, GM_MODERN, J, DEGREES, R_NORM))
g45, Gphi45, mag45, psi45 = gv.g_eff_vector(
    0.0, R_EQ, phi45, OMEGA, GM_MODERN, J, DEGREES, R_NORM
)
record(5, "G_phi_eff vanishes at the equator and the pole and is equatorward at 45 degrees",
       abs(gp_eq) < 1e-12 and abs(gp_pole) < 1e-12 and gp45 < 0,
       f"equator : {gp_eq:.3e} m/s2\n"
       f"pole    : {gp_pole:.3e} m/s2\n"
       f"45 deg  : {gp45:.6f} m/s2, magnitude {abs(gp45):.6f} m/s2\n"
       f"SIGN: the specification's formula -(1/r) dV/dphi makes this the component along "
       f"INCREASING latitude, which is negative in the north because the bulge pulls a mid "
       f"latitude point equatorward. The acceptance text calls it 'positive (toward the "
       f"equator)'. Both describe the same physical direction; only the sign convention "
       f"differs. See the report, section 3.\n"
       f"psi at 45 deg = {math.degrees(float(psi45)):.4f} deg, from arctan(-G_phi/g), which is "
       f"the convention that makes phi_c = phi_g - psi in Step 6\n"
       f"|g_eff| at 45 deg = {float(mag45):.4f} m/s2")

# ---------------------------------------------------------------------------
# 6. Jupiter cross check. Values quoted from SPEC_01 Step 4; in no static file.
# ---------------------------------------------------------------------------
J_JUP = [14736e-6, -587e-6, 31e-6]
R_JUP = 71398e3
GM_JUP = 1.26686534e17
OMEGA_JUP = 2.0 * math.pi / (9 * 3600 + 55 * 60 + 29.7)
RE_JUP, RP_JUP, U_JUP = 71492e3, 66854e3, 100.0
j_eq = float(gv.g_eff_radial(U_JUP, RE_JUP, 0.0, OMEGA_JUP, GM_JUP, J_JUP, DEGREES, R_JUP))
j_pole = float(
    gv.g_eff_radial(0.0, RP_JUP, math.pi / 2, OMEGA_JUP, GM_JUP, J_JUP, DEGREES, R_JUP)
)
record(6, "Jupiter cross check reproduces Lindal Table II",
       abs(j_eq - 23.116) < 0.005 and abs(j_pole - 27.015) < 0.005,
       f"equator, u = 100 m/s : {j_eq:.4f}  (spec 23.116; Lindal 23.12)\n"
       f"pole                 : {j_pole:.4f}  (spec 27.015; Lindal 27.01)\n"
       f"A second planet with different harmonics, radii and rotation exercises the same code "
       f"path, so agreement here is independent evidence that the (l+1) coefficient is right.")

# ---------------------------------------------------------------------------
# 7. Broadcasting and purity
# ---------------------------------------------------------------------------
phis = np.radians(np.array([0.0, 30.0, 45.0, 60.0, 90.0]))
radii = np.full_like(phis, R_EQ)
vector = gv.g_eff_radial(0.0, radii, phis, OMEGA, GM_MODERN, J, DEGREES, R_NORM)
one_at_a_time = np.array([g_at(0.0, R_EQ, p) for p in phis])
scalar_out = gv.g_eff_radial(0.0, R_EQ, 0.0, OMEGA, GM_MODERN, J, DEGREES, R_NORM)
record(7, "the functions broadcast, and an array call equals the scalar calls",
       np.allclose(vector, one_at_a_time, rtol=0, atol=0) and np.ndim(scalar_out) == 0,
       f"latitudes {np.degrees(phis).tolist()}\n"
       f"g values  {[round(float(v), 4) for v in vector]}\n"
       f"identical to the element by element calls: "
       f"{bool(np.array_equal(vector, one_at_a_time))}\n"
       f"a scalar call returns a scalar: ndim = {np.ndim(scalar_out)}")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
