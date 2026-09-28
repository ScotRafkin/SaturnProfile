"""Acceptance checks for SPEC_03 v0.6 Step 1, `lib.geopotential`.

Every check prints its measured value. Checks the specification did not ask for are labeled
"beyond the specification". The Lindal values are computed on the clean product of the Step 0
sweep (`occul_data/lindal/lindal_refractivity.nc`, commit 2149b64) with the gravity, rotation and
wind the reduction used, loaded through `lib.control`. Nothing is written.
"""

import math
import sys

import numpy as np

from casspian.lib import control as ctl
from casspian.lib import geopotential as gp
from casspian.lib import io as cio
from casspian.lib.gravity import g_eff_radial
from casspian.refrac.anchor import geoid_setup, wind_of_latitude

D = "occul_data/lindal/"
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def figures(value, spec_text):
    """True when `value` printed in the format of `spec_text` reads as `spec_text`."""
    mantissa = spec_text.lower().split("e")[0]
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    fmt = f"{{:.{decimals}e}}" if "e" in spec_text.lower() else f"{{:.{decimals}f}}"
    return float(fmt.format(value)) == float(spec_text)


# ---------------------------------------------------------------------------
# The Lindal column
# ---------------------------------------------------------------------------
tree = cio.read(D + "lindal_refractivity.nc", "refractivity")
root = tree.to_dataset(inherit=False)
thermo = tree["inputs/thermo"].to_dataset(inherit=False)
commit = str(tree.attrs["casspian_git_commit"])
u_record = float(tree["reduction_record"].attrs["u_at_anchor_ms"])
r = np.asarray(root["radius_m"].values, dtype="float64")
h = np.asarray(root["height_above_anchor_isobar_m"].values, dtype="float64")
p = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
phi_c = math.radians(float(root["latitude_planetocentric_deg"].values))
psi_anchor = math.radians(float(root["psi_deg"].values))
tree.close()

manifest = ctl.read_reduction_manifest(D + "lindal_reduction.toml")
inputs = ctl.load_reduction_inputs(manifest)
constants = geoid_setup(inputs, manifest).constants
u = float(wind_of_latitude(inputs.wind)(np.array([phi_c]))[0])
a = int(np.flatnonzero(p == 1.0e4)[0])
prof = gp.geopotential_along_profile(u, r, h, phi_c, a, *constants)
Phi, dPhi, gm, psi = prof.geopotential_m2s2, prof.layer_increments_m2s2, prof.g_magnitude_ms2, prof.psi_rad
print(f"product commit {commit}; {r.size} levels; gauge level {a} at {p[a]!r} Pa; "
      f"u(phi_c) = {u!r} m/s (record {u_record!r}); phi_c = {math.degrees(phi_c)!r} deg\n")

# ---------------------------------------------------------------------------
# 1. Uniform gravity
# ---------------------------------------------------------------------------
R_BIG = 1.0e18
GM_UNIFORM = 10.0 * R_BIG**2
r_uniform = R_BIG + (h - h[a])
g_levels = GM_UNIFORM / r_uniform**2
variation = float(g_levels.max() / g_levels.min() - 1.0)
uni = gp.geopotential_along_profile(0.0, r_uniform, h, phi_c, a, 0.0, GM_UNIFORM,
                                    np.array([0.0]), np.array([2]), 6.0e7)
g_uniform = GM_UNIFORM / R_BIG**2
expected = g_uniform * (h - h[a])
off = np.arange(h.size) != a
relative = np.abs(uni.geopotential_m2s2[off] - expected[off]) / np.abs(expected[off])
record(1, "uniform gravity returns Phi_k - Phi_a = g (h_k - h_a) to 1e-12 relative",
       variation < 1e-12 and float(relative.max()) < 1e-12 and uni.geopotential_m2s2[a] == 0.0,
       f"r = 1e18 m + (h - h_a), GM = {GM_UNIFORM:.3e}, J = 0, Omega = 0, u = 0; GM / r^2 varies by "
       f"{variation:.2e} across the column\n"
       f"largest relative departure over the {int(off.sum())} levels off the gauge: "
       f"{relative.max():.2e}; Phi_a = {uni.geopotential_m2s2[a]!r}")

# ---------------------------------------------------------------------------
# 2. Phi_a is exactly zero; 3. Phi strictly monotonic in h
# ---------------------------------------------------------------------------
record(2, "Phi at the gauge level is exactly 0.0", Phi[a] == 0.0,
       f"Lindal: Phi[{a}] = {Phi[a]!r}; uniform case: {uni.geopotential_m2s2[a]!r}")
same_sense = np.sign(np.diff(Phi)) == np.sign(np.diff(h))
record(3, "Phi is strictly monotonic in h", bool(np.all(same_sense)) and bool(np.all(np.diff(Phi) != 0)),
       f"sign(diff Phi) equals sign(diff h) at {int(same_sense.sum())} of {same_sense.size} layers; "
       f"Phi decreases with level index (upward positive): {bool(np.all(np.diff(Phi) < 0))}")

# ---------------------------------------------------------------------------
# 4. The expected values
# ---------------------------------------------------------------------------
pairs = [
    ("|g_eff| top (m/s2)", gm[0], "9.892477"),
    ("|g_eff| anchor (m/s2)", gm[a], "10.005625"),
    ("|g_eff| bottom (m/s2)", gm[-1], "10.047103"),
    ("psi top (deg)", math.degrees(psi[0]), "5.546374"),
    ("psi bottom (deg)", math.degrees(psi[-1]), "5.475835"),
    ("Phi top (m2/s2)", Phi[0], "2.852355e6"),
    ("Phi bottom (m2/s2)", Phi[-1], "-1.043743e6"),
    ("first layer dPhi (m2/s2)", dPhi[0], "-1.2666e5"),
]
rows = [(name, value, text, figures(value, text)) for name, value, text in pairs]
record(4, "the expected values reproduced to the figures printed",
       all(ok for *_, ok in rows),
       "\n".join(f"{name:26s} {value!r:24s} expected {text:12s} {'agrees' if ok else 'DIFFERS'}"
                 for name, value, text, ok in rows))

# The radial rule, in this script only (no code path in lib). Two readings: the radial component
# times the tabulated altitude increment (the pre-Step 0 reading of B1 and B2, where h was taken as
# radial), and the radial component times the radial increment r_{k+1} - r_k, which after Step 0 is
# (h_{k+1} - h_k) cos psi at the anchor, the form SPEC_03 Step 4 writes for its negative control.
g_radial = prof.g_radial_ms2
g_layer = 0.5 * (g_radial[:-1] + g_radial[1:])
radial_h = gp.geopotential(g_layer * np.diff(h), a)
radial_r = gp.geopotential(g_layer * np.diff(r), a)
record("4b", "beyond the specification: which radial reading gives the specification's "
             "2.839126e6 and -1.038963e6",
       figures(radial_h[0], "2.839126e6") and figures(radial_h[-1], "-1.038963e6"),
       f"g_k with (h_(k+1) - h_k):          top {radial_h[0]!r}, bottom {radial_h[-1]!r}; "
       f"field line over it {Phi[0] / radial_h[0] - 1:+.4e}, {Phi[-1] / radial_h[-1] - 1:+.4e}\n"
       f"g_k with (h_(k+1) - h_k) cos psi:  top {radial_r[0]!r}, bottom {radial_r[-1]!r}; "
       f"field line over it {Phi[0] / radial_r[0] - 1:+.4e}, {Phi[-1] / radial_r[-1] - 1:+.4e}\n"
       f"expected 2.839126e6 and -1.038963e6; 1 / cos psi_anchor - 1 = "
       f"{1 / math.cos(psi_anchor) - 1:.4e}, 1 / cos^2 psi_anchor - 1 = "
       f"{1 / math.cos(psi_anchor) ** 2 - 1:.4e}")

# ---------------------------------------------------------------------------
# 5. |g_eff| at the anchor against g_eff_radial / cos psi
# ---------------------------------------------------------------------------
g_direct = float(g_eff_radial(u, r[a], phi_c, *constants))
ratio = gm[a] / (g_direct / math.cos(psi[a])) - 1.0
record(5, "|g_eff| at the anchor equals lib.gravity.g_eff_radial there divided by cos psi to round-off",
       abs(ratio) < 1e-14,
       f"|g_eff| {gm[a]!r}; g_eff_radial {g_direct!r}; / cos psi {g_direct / math.cos(psi[a])!r}; "
       f"relative difference {ratio:.2e}\n"
       f"psi at the anchor level {math.degrees(psi[a])!r} deg; the product's frozen psi "
       f"{math.degrees(psi_anchor)!r} deg")

# ---------------------------------------------------------------------------
# 6. geopotential refuses a gauge index off the grid
# ---------------------------------------------------------------------------
lines, ok = [], True
for bad in (-1, r.size, r.size + 5, 2.5, True, np.float64(29.0)):
    try:
        gp.geopotential(dPhi, bad)
        ok = False
        lines.append(f"gauge_index {bad!r}: not refused")
    except ValueError as exc:
        lines.append(f"gauge_index {bad!r}: {exc}")
accepted = gp.geopotential(dPhi, np.int64(a))
lines.append(f"gauge_index np.int64({a}): accepted, equal to the int result "
             f"{np.array_equal(accepted, Phi)}")
record(6, "geopotential refuses a gauge index off the grid", ok and np.array_equal(accepted, Phi),
       "\n".join(lines))

# ---------------------------------------------------------------------------
# 7. Beyond the specification: layer_increments refuses a malformed column
# ---------------------------------------------------------------------------
lines, ok = [], True
h_flat = h.copy(); h_flat[5] = h_flat[4]
h_back = h.copy(); h_back[5], h_back[6] = h_back[6], h_back[5]
for label, args in (("a zero thickness layer", (gm, h_flat)), ("a reversed layer", (gm, h_back)),
                    ("a shape mismatch", (gm[:-1], h)), ("one level", (gm[:1], h[:1])),
                    ("a NaN height", (gm, np.where(np.arange(h.size) == 3, np.nan, h)))):
    try:
        gp.layer_increments(*args)
        ok = False
        lines.append(f"{label}: not refused")
    except ValueError as exc:
        lines.append(f"{label}: {exc}")
record(7, "beyond the specification: layer_increments refuses a malformed column", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 8. Beyond the specification: the trapezoid error, measured
# ---------------------------------------------------------------------------
# Simpson on each layer with |g_eff| at the midpoint radius; the difference from the trapezoid is
# the trapezoid error to leading order.
r_mid, h_mid = 0.5 * (r[:-1] + r[1:]), 0.5 * (h[:-1] + h[1:])
g_mid = gp.effective_gravity_magnitude(u, r_mid, phi_c, *constants)[0]
simpson = (gm[:-1] + 4.0 * g_mid + gm[1:]) / 6.0 * np.diff(h)
layer_error = np.abs(dPhi / simpson - 1.0)
thick = int(np.argmax(np.abs(np.diff(h))))
bound = np.diff(h) ** 2 / (2.0 * r[:-1] ** 2)
Phi_simpson = gp.geopotential(simpson, a)
record(8, "beyond the specification: the trapezoid error per layer, against Simpson, is below "
          "1e-7 of the layer as the staggering table states; compared with its estimate "
          "(dh)^2 / (2 r^2)",
       bool(np.all(layer_error < 1e-7)),
       f"largest per layer |trapezoid / Simpson - 1| = {layer_error.max():.2e} at layer "
       f"{int(np.argmax(layer_error))}; the thickest layer ({abs(np.diff(h)[thick]) / 1e3:.2f} km, "
       f"layer {thick}) {layer_error[thick]:.2e} against the bound {bound[thick]:.2e}\n"
       f"largest ratio of the measured error to the bound: {np.max(layer_error / bound):.3f}\n"
       f"effect on Phi at the top {Phi[0] - Phi_simpson[0]:+.3e} m2/s2 "
       f"({Phi[0] / Phi_simpson[0] - 1:+.2e}), at the bottom {Phi[-1] - Phi_simpson[-1]:+.3e} "
       f"({Phi[-1] / Phi_simpson[-1] - 1:+.2e})\n"
       f"|g_eff| varies by {gm.max() / gm.min() - 1:.4f} over the profile")

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
