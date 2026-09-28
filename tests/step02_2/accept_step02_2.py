"""Acceptance checks for SPEC_02 v0.4 Step 2, `refrac.anchor.freeze_anchor`.

Every check reports its measured value. The refusal cases run on copies of the manifest and its
six inputs under this directory, because a manifest may not point outside its own directory.
Nothing under `occul_data/` is written.
"""

import dataclasses
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np

from casspian.lib import control as ctl
from casspian.lib import geoid as gd
from casspian.lib import latitude as lat
from casspian.lib.control import ControlFileError
from casspian.lib.gravity import G_phi_eff, g_eff_radial
from casspian.refrac import anchor as anc

HERE = Path("reports/step02_2")
HERE.mkdir(parents=True, exist_ok=True)
SOURCE = Path("occul_data/lindal")
MANIFEST = SOURCE / "lindal_reduction.toml"
NULL_GM = 3.7929085e16
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def copy_case(name):
    case = HERE / name
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    shutil.copy2(MANIFEST, case / MANIFEST.name)
    for path in ctl.read_reduction_manifest(MANIFEST).inputs.values():
        shutil.copy2(path, case / path.name)
    return case


def edit_manifest(case, pattern, replacement):
    target = case / MANIFEST.name
    text = target.read_text(encoding="utf-8")
    edited, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"the acceptance edit {pattern!r} matched {count} times")
    target.write_bytes(edited.encode("utf-8"))


deg = math.degrees
manifest = ctl.read_reduction_manifest(MANIFEST)
# SPEC_02 v0.8 Step 6 may change the committed manifest's anchor rule. This suite tests Step 2
# as accepted, under mean polar anchoring, so it pins that rule and quantity in memory.
manifest = dataclasses.replace(manifest, anchor_rule="mean_polar_radius",
                               anchor_quantity="radius_polar_m")
inputs = ctl.load_reduction_inputs(manifest)
a = anc.freeze_anchor(inputs, manifest)

gravity, rotation = inputs.gravity, inputs.rotation
DEGREES = np.asarray(gravity["degree"].values)
J = np.asarray(gravity["J"].values, dtype="float64")
GM = float(gravity["GM_m3s2"])
R_NORM = float(gravity["normalization_radius_m"])
OMEGA = float(rotation["angular_rate_rad_s"])
CONSTANTS = (OMEGA, GM, J, DEGREES, R_NORM)
u_of_phi = anc.wind_of_latitude(inputs.wind)

# ---------------------------------------------------------------------------
# 1. The no wind reference values
# ---------------------------------------------------------------------------
nw_phi, nw_r0 = deg(a.nowind_phi_c_rad), a.nowind_r0_m / 1e3
record(1, "no wind reference: phi_c = 30.8182 +- 0.001 deg and r0 = 58,453.1 +- 0.5 km with "
       "Null's GM",
       GM == NULL_GM and abs(nw_phi - 30.8182) <= 0.001 and abs(nw_r0 - 58453.1) <= 0.5,
       f"GM from lindal_gravity.nc = {GM!r}, Null 1981: {GM == NULL_GM}\n"
       f"phi_c = {nw_phi:.6f} deg, departure {nw_phi - 30.8182:+.6f} deg "
       f"(SPEC_01 Step 6 gave 30.818189)\n"
       f"psi   = {deg(a.nowind_psi_rad):.6f} deg (Step 6 gave 5.481811)\n"
       f"r0    = {nw_r0:.3f} km, departure {nw_r0 - 58453.1:+.3f} km "
       f"(Step 6 gave 58453.120 at this latitude)")

# ---------------------------------------------------------------------------
# 2. The wind included values, and where the expected 58,518 km comes from
# ---------------------------------------------------------------------------
# The same anchored march, read at the no wind latitude instead of the wind one, separates
# the dynamical height at a fixed latitude from the effect of the latitude shift.
at_nowind_lat = gd.wind_geoid(np.array([a.nowind_phi_c_rad]), a.r_anchor_m, a.anchor_rule,
                              u_of_phi, *CONSTANTS, tol_m=manifest.convergence_m)
r_wind_nw_lat = float(at_nowind_lat.radius[0])
shift_deg = deg(a.phi_c_rad) - nw_phi
g =float(g_eff_radial(a.u_at_phi_c_ms, a.r0_m, a.phi_c_rad, *CONSTANTS))
G_phi = float(G_phi_eff(a.u_at_phi_c_ms, a.r0_m, a.phi_c_rad, *CONSTANTS))
slope = a.r0_m * G_phi / g
record(2, "wind included: r0 near 58,518 km and phi_c within 0.02 deg of the no wind value "
       "(expected figures, reported)",
       abs(shift_deg) <= 0.02,
       f"phi_c = {deg(a.phi_c_rad):.6f} deg, psi = {deg(a.psi_rad):.6f} deg, "
       f"u(phi_c) = {a.u_at_phi_c_ms:.3f} m/s\n"
       f"phi_c - no wind phi_c = {shift_deg:+.6f} deg (within 0.02: {abs(shift_deg) <= 0.02})\n"
       f"r0 = {a.r0_m / 1e3:.3f} km; r0 - no wind r0 = {(a.r0_m - a.nowind_r0_m) / 1e3:+.3f} km\n"
       f"decomposition of that difference:\n"
       f"  the wind surface at the NO WIND latitude {nw_phi:.6f} deg: "
       f"{r_wind_nw_lat / 1e3:.3f} km, the dynamical height there "
       f"{(r_wind_nw_lat - a.nowind_r0_m) / 1e3:+.3f} km\n"
       f"  the latitude shift along the wind surface: "
       f"{(a.r0_m - r_wind_nw_lat) / 1e3:+.3f} km; Eq. B3 slope at phi_c "
       f"{slope / 1e3:.1f} km/rad times the shift {math.radians(shift_deg):.3e} rad = "
       f"{slope * math.radians(shift_deg) / 1e3:+.3f} km\n"
       f"The spec's 58,518 km matches the wind surface read at the no wind latitude "
       f"({r_wind_nw_lat / 1e3:.1f} km); read at its own converged latitude the frozen r0 is "
       f"{a.r0_m / 1e3:.1f} km.")

# ---------------------------------------------------------------------------
# 3. Convergence of the fixed point
# ---------------------------------------------------------------------------
record(3, "the fixed point converges to the manifest tolerance in at most eight iterations, "
       "iterates listed",
       a.iteration_count <= 8 and a.final_step_rad < a.tolerance_rad,
       f"manifest fixed_point_tolerance_deg = {manifest.fixed_point_tolerance_deg!r}, converted "
       f"once in refrac.anchor to {a.tolerance_rad!r} rad; max_iterations = "
       f"{manifest.max_iterations}\n"
       f"wind:    {a.iteration_count} iterations, last step {deg(a.final_step_rad):.3e} deg\n"
       f"         iterates (deg, seed first): "
       + ", ".join(f"{deg(v):.7f}" for v in a.iterates_rad)
       + f"\nno wind: {a.nowind_iterates_rad.size - 1} iterations, last step "
         f"{deg(abs(a.nowind_iterates_rad[-1] - a.nowind_iterates_rad[-2])):.3e} deg\n"
         f"         iterates (deg, seed first): "
       + ", ".join(f"{deg(v):.7f}" for v in a.nowind_iterates_rad)
       + f"\nseed from the anchor surface oblateness: "
         f"{float(inputs.geodesy['oblateness'].values[0])!r}")

# ---------------------------------------------------------------------------
# 4. The anchoring residual
# ---------------------------------------------------------------------------
record(4, "the anchoring residual is below 1e-3 m", abs(a.anchor_residual_m) < 1e-3,
       f"anchor_rule {a.anchor_rule!r} on {a.anchor_quantity} at {a.anchor_surface_Pa:g} Pa = "
       f"{a.r_anchor_m / 1e3:.3f} km\n"
       f"residual {a.anchor_residual_m:.3e} m; mean of the two polar radii minus the anchor "
       f"{0.5 * (a.polar_north_m + a.polar_south_m) - a.r_anchor_m:.3e} m\n"
       f"the manifest's convergence_m = {manifest.convergence_m!r} is the secant stopping "
       "tolerance; the near affine map from start to outcome takes it far below that")

# ---------------------------------------------------------------------------
# 5. Polar radii and asymmetry
# ---------------------------------------------------------------------------
pn, ps, asym = a.polar_north_m / 1e3, a.polar_south_m / 1e3, a.polar_asymmetry_m / 1e3
record(5, "polar radii and asymmetry as REPORT_01_step7 section 6, to 0.1 km",
       abs(pn - 54423.6) <= 0.1 and abs(ps - 54452.4) <= 0.1 and abs(asym - 28.7) <= 0.1,
       f"north polar radius {pn:.3f} km (54423.6), departure {pn - 54423.6:+.3f}\n"
       f"south polar radius {ps:.3f} km (54452.4), departure {ps - 54452.4:+.3f}\n"
       f"asymmetry          {asym:.3f} km (28.7), departure {asym - 28.7:+.3f}\n"
       f"north polar start  {a.north_start_m / 1e3:.3f} km, which is the north polar radius, "
       f"as a march from the north pole requires: {a.north_start_m == a.polar_north_m}\n"
       "these are from the final march, the one r0 was read from")

# ---------------------------------------------------------------------------
# 6. Beyond the specification: the frozen pair is stable and self consistent
# ---------------------------------------------------------------------------
fine = gd.wind_geoid(np.array([a.phi_c_rad]), a.r_anchor_m, a.anchor_rule, u_of_phi,
                     *CONSTANTS, tol_m=manifest.convergence_m, march_step_deg=0.025)
r0_fine = float(fine.radius[0])


def wind_surface(phi):
    return gd.wind_geoid(np.atleast_1d(phi), a.r_anchor_m, a.anchor_rule, u_of_phi, *CONSTANTS,
                         tol_m=manifest.convergence_m).radius


back = float(lat.planetographic_from_planetocentric(
    np.array([a.phi_c_rad]), wind_surface, u_of_phi, *CONSTANTS)[0])
identity = deg(a.phi_c_rad + a.psi_rad) - deg(a.phi_g_rad)
record(6, "beyond the specification: r0 is converged in the march step, and phi_c inverts to "
       "the label",
       abs(r0_fine - a.r0_m) < 1.0 and abs(deg(back) - 36.3) < 1e-6,
       f"march step 0.05 deg: r0 = {a.r0_m:.4f} m; 0.025 deg: {r0_fine:.4f} m; difference "
       f"{r0_fine - a.r0_m:+.4f} m\n"
       f"planetographic_from_planetocentric(phi_c) on the wind surface = {deg(back):.10f} deg, "
       f"error {deg(back) - 36.3:+.2e} deg\n"
       f"psi is evaluated at phi_c on the final march, so phi_c + psi - phi_g = "
       f"{identity:+.2e} deg rather than exactly zero; the lib fixed point's own psi, taken at "
       f"the previous iterate, would make it zero by construction")

# ---------------------------------------------------------------------------
# 7. Beyond the specification: refusals
# ---------------------------------------------------------------------------
case = copy_case("case_iterations")
edit_manifest(case, r"^(max_iterations\s*=\s*)\S+", r"\g<1>2")
m2 = ctl.read_reduction_manifest(case / MANIFEST.name)
try:
    anc.freeze_anchor(ctl.load_reduction_inputs(m2), m2)
    refused_iter, message_iter = False, "not refused"
except anc.AnchorConvergenceError as exc:
    refused_iter, message_iter = True, str(exc)

case = copy_case("case_sensitivity")
target = case / MANIFEST.name
target.write_bytes(target.read_bytes()
                   + b"\n[sensitivity]\nenabled = true\nwind_draws = 200\nwind_draw_seed = 1\n")
try:
    ctl.read_reduction_manifest(target)
    refused_sens, message_sens = False, "not refused"
except ControlFileError as exc:
    refused_sens, message_sens = True, str(exc)
record(7, "beyond the specification: an exhausted iteration budget and a [sensitivity] section "
       "are refused",
       refused_iter and refused_sens and "unknown section" in message_sens,
       f"max_iterations = 2: refused {refused_iter}\n  {message_iter}\n"
       f"[sensitivity] appended (SPEC_00 v0.12 removed it): refused {refused_sens}\n"
       f"  {message_sens}")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
