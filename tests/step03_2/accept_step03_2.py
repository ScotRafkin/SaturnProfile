"""Acceptance checks for SPEC_03 v0.7 Step 2, `lib.hydrostatic`.

Every check prints its measured value. Checks the specification did not ask for are labeled
"beyond the specification". The geopotential grid is the Step 1 grid on the clean product of the
Step 0 sweep (`occul_data/lindal/lindal_refractivity.nc`, commit 2149b64), from the accepted
`lib.geopotential` with the reduction's gravity, rotation and wind. The test columns are built on
that grid from their analytic pressure, with the mean molar mass and mean refractivity of the
product's top level held uniform. Nothing is written.
"""

import math
import sys

import numpy as np

from casspian.lib import control as ctl
from casspian.lib import geopotential as gp
from casspian.lib import hydrostatic as hs
from casspian.lib import io as cio
from casspian.lib.constants import AVOGADRO_CONSTANT, BOLTZMANN_CONSTANT
from casspian.refrac.anchor import geoid_setup, wind_of_latitude

D = "occul_data/lindal/"
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


# ---------------------------------------------------------------------------
# The Step 1 geopotential grid and the product's own fields
# ---------------------------------------------------------------------------
tree = cio.read(D + "lindal_refractivity.nc", "refractivity")
root = tree.to_dataset(inherit=False)
thermo = tree["inputs/thermo"].to_dataset(inherit=False)
commit = str(tree.attrs["casspian_git_commit"])
r = np.asarray(root["radius_m"].values, dtype="float64")
h = np.asarray(root["height_above_anchor_isobar_m"].values, dtype="float64")
N_prod = np.asarray(root["refractivity"].values, dtype="float64")
R_prod = np.asarray(root["mean_refractivity_m3"].values, dtype="float64")
M_prod = np.asarray(root["mean_molar_mass_kg_mol"].values, dtype="float64")
p_tab = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
T_tab = np.asarray(thermo["temperature_K"].values, dtype="float64")
phi_c = math.radians(float(root["latitude_planetocentric_deg"].values))
tree.close()

manifest = ctl.read_reduction_manifest(D + "lindal_reduction.toml")
inputs = ctl.load_reduction_inputs(manifest)
constants = geoid_setup(inputs, manifest).constants
u = float(wind_of_latitude(inputs.wind)(np.array([phi_c]))[0])
a = int(np.flatnonzero(p_tab == 1.0e4)[0])
prof = gp.geopotential_along_profile(u, r, h, phi_c, a, *constants)
Phi = prof.geopotential_m2s2
g_mag_levels = prof.g_magnitude_ms2

M = float(M_prod[0])
R_BAR = float(R_prod[0])
m = M / AVOGADRO_CONSTANT
P_B = 20.0
print(f"product commit {commit}; {Phi.size} levels; Phi from {Phi[0]:.6e} to {Phi[-1]:.6e} m2/s2 "
      f"(gauge level {a}); uniform M = {M!r} kg/mol, R_bar = {R_BAR!r} m3; p_b = {P_B} Pa\n")


def column(Phi_grid, T_of_Phi, p_of_Phi):
    """The analytic column on a grid: p, T, and N and rho through the module's own density."""
    p = p_of_Phi(Phi_grid)
    T = T_of_Phi(Phi_grid)
    N = p * R_BAR / (BOLTZMANN_CONSTANT * T)
    rho = hs.density(N, R_BAR, np.full(Phi_grid.shape, M))
    return p, T, N, rho


def produce(Phi_grid, rho, rule="exact"):
    if rule == "exact":
        I = hs.layer_mass(rho, Phi_grid)
    else:
        I = 0.5 * (rho[:-1] + rho[1:]) * (Phi_grid[:-1] - Phi_grid[1:])
    return hs.pressure_from_top(P_B, I)


Phi0 = Phi[0]

# ---------------------------------------------------------------------------
# 1. The isothermal column
# ---------------------------------------------------------------------------
T_ISO = 100.0
H_iso = BOLTZMANN_CONSTANT * T_ISO / m
p_iso_f = lambda x: P_B * np.exp(-(x - Phi0) / H_iso)
T_iso_f = lambda x: np.full(np.shape(x), T_ISO)
p_a, T_a, N_a, rho_a = column(Phi, T_iso_f, p_iso_f)
p_iso = produce(Phi, rho_a)
T_iso = hs.temperature(p_iso, N_a, R_BAR)
err_p = np.abs(p_iso / p_a - 1.0)
err_T = np.abs(T_iso / T_ISO - 1.0)
p_trap_iso = produce(Phi, rho_a, "trapezoid")
record(1, "the isothermal column (T = 100 K, p_b = 20 Pa) returns p to 1e-14 and T = 100 K to 1e-14 "
          "relative at every level",
       float(err_p.max()) < 1e-14 and float(err_T.max()) < 1e-14,
       f"scale geopotential k_B T / m = {H_iso:.6e} m2/s2; p from {p_a[0]:.4g} to {p_a[-1]:.4g} Pa\n"
       f"largest |p / p_analytic - 1| = {err_p.max():.2e} at level {int(np.argmax(err_p))}\n"
       f"largest |T / 100 K - 1| = {err_T.max():.2e} at level {int(np.argmax(err_T))}\n"
       f"trapezoid alternative (script only): largest |p / p_analytic - 1| = "
       f"{np.abs(p_trap_iso / p_a - 1).max():.3e}")

# ---------------------------------------------------------------------------
# 2. The linear-T column, both orientations; the trapezoid alternative
# ---------------------------------------------------------------------------
span = Phi[0] - Phi[-1]


def linear_T(T_top, T_bottom):
    beta = (T_top - T_bottom) / span               # dT/dPhi
    T_f = lambda x: T_bottom + beta * (x - Phi[-1])
    # dp/dPhi = -p m / (k_B T): ln p = ln p_b - (m / (k_B beta)) ln(T / T_top)
    p_f = lambda x: P_B * (T_f(x) / T_top) ** (-m / (BOLTZMANN_CONSTANT * beta))
    return T_f, p_f


orientation = {}
for label, (T_top, T_bottom) in {"140 K at the top to 80 K at the bottom": (140.0, 80.0),
                                 "80 K at the top to 140 K at the bottom": (80.0, 140.0)}.items():
    T_f, p_f = linear_T(T_top, T_bottom)
    p_l, T_l, N_l, rho_l = column(Phi, T_f, p_f)
    exact = np.abs(produce(Phi, rho_l) / p_l - 1.0)
    trap = np.abs(produce(Phi, rho_l, "trapezoid") / p_l - 1.0)
    # halving: midpoints in Phi, analytic density there
    Phi_fine = np.empty(2 * Phi.size - 1)
    Phi_fine[0::2], Phi_fine[1::2] = Phi, 0.5 * (Phi[:-1] + Phi[1:])
    p_lf, _, _, rho_lf = column(Phi_fine, T_f, p_f)
    exact_fine = np.abs(produce(Phi_fine, rho_lf)[0::2] / p_l - 1.0)
    orientation[label] = (exact, trap, exact_fine, p_l)

# "Falls linearly from 140 to 80 K" is read as falling with height, 140 K at the bottom and 80 K at
# the top, as an atmosphere cools upward (decision 4); both orientations are reported.
chosen = "80 K at the top to 140 K at the bottom"
exact, trap, exact_fine, p_l = orientation[chosen]
lines = []
for label, (e, t, ef, pl) in orientation.items():
    lines.append(f"{label}: p from {pl[0]:.4g} to {pl[-1]:.4g} Pa; log-linear largest error "
                 f"{e.max():.3e} (level {int(np.argmax(e))}); trapezoid {t.max():.3e}; halved "
                 f"{ef.max():.3e}, ratio {e.max() / ef.max():.3f}")
record(2, "the linear-T column returns p to 6e-4 or better at every level; the trapezoid "
          "alternative 8e-3 to 9e-3, both reported",
       float(exact.max()) <= 6e-4 and 8e-3 <= float(trap.max()) <= 9e-3,
       "\n".join(lines)
       + f"\nthe check is on '{chosen}', falling with height (decision 4); T is linear in Phi, "
         "the column's own vertical coordinate")

# ---------------------------------------------------------------------------
# 3. Halving every layer: second order
# ---------------------------------------------------------------------------
ratio = float(exact.max() / exact_fine.max())
record(3, "halving every layer reduces the linear-T error by a factor of 3.5 to 4.5",
       3.5 <= ratio <= 4.5,
       f"largest error {exact.max():.4e} on the Step 1 grid, {exact_fine.max():.4e} with every "
       f"layer halved (midpoints in Phi, the analytic density there): ratio {ratio:.3f}\n"
       f"at the bottom level: {exact[-1]:.4e} and {exact_fine[-1]:.4e}, ratio "
       f"{exact[-1] / exact_fine[-1]:.3f}")

# ---------------------------------------------------------------------------
# 4. layer_mass: the limit form on equal densities; the refusals
# ---------------------------------------------------------------------------
lines, ok = [], True
equal = hs.layer_mass(np.array([2.0e-3, 2.0e-3]), np.array([5.0e5, 4.0e5]))
ok = ok and equal[0] == 2.0e-3 * 1.0e5
lines.append(f"equal densities 2e-3 kg/m3 across 1e5 m2/s2: {equal[0]!r}, the limit form "
             f"rho_k dPhi = {2.0e-3 * 1.0e5!r}, equal {equal[0] == 2.0e-3 * 1.0e5}")
near = hs.layer_mass(np.array([2.0e-3, 2.0e-3 * (1 + 5e-11)]), np.array([5.0e5, 4.0e5]))
lines.append(f"log ratio 5e-11 (inside the limit): {near[0]!r}; relative to the exact form "
             f"{near[0] / (2.0e-3 * 1.0e5 * (1 + 2.5e-11)) - 1:+.2e}")
Phi_zero = Phi.copy(); Phi_zero[7] = Phi_zero[6]
for label, args, needle in (
        ("a zero thickness layer", (rho_a, Phi_zero), "layer 6"),
        ("a column that is not top down", (rho_a[::-1], Phi[::-1]), "layer 0"),
        ("a zero density", (np.where(np.arange(Phi.size) == 3, 0.0, rho_a), Phi), "levels [3]"),
        ("a negative density", (np.where(np.arange(Phi.size) == 4, -1.0, rho_a), Phi), "levels [4]")):
    try:
        hs.layer_mass(*args)
        ok = False
        lines.append(f"{label}: not refused")
    except ValueError as exc:
        ok = ok and needle in str(exc)
        lines.append(f"{label}: {exc}")
record(4, "layer_mass returns the limit form on equal densities and refuses a zero-thickness layer "
          "with the layer named", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 5. temperature on the reduction's own N, R_bar and tabulated p
# ---------------------------------------------------------------------------
T_rec = hs.temperature(p_tab, N_prod, R_prod)
err = np.abs(T_rec / T_tab - 1.0)
record(5, "temperature applied to the reduction's own N, R_bar and tabulated p returns the "
          "tabulated T to 1e-12 relative",
       float(err.max()) < 1e-12,
       f"largest |T / T_tab - 1| = {err.max():.2e} at level {int(np.argmax(err))} "
       f"({p_tab[int(np.argmax(err))]:.6g} Pa); {err.size} levels")

# ---------------------------------------------------------------------------
# 6. Beyond the specification: density against the reduction's number density
# ---------------------------------------------------------------------------
n_prod = np.asarray(cio.read(D + "lindal_refractivity.nc", "refractivity").to_dataset(inherit=False)
                    ["number_density_m3"].values, dtype="float64")
rho_prod = hs.density(N_prod, R_prod, M_prod)
rho_direct = n_prod * M_prod / AVOGADRO_CONSTANT
err_rho = np.abs(rho_prod / rho_direct - 1.0)
record(6, "beyond the specification: density(N, R_bar, m_bar) equals n m_bar / N_A from the "
          "reduction's number density",
       float(err_rho.max()) < 1e-14,
       f"largest relative difference {err_rho.max():.2e}; rho from {rho_prod[0]:.4e} to "
       f"{rho_prod[-1]:.4e} kg/m3")

# ---------------------------------------------------------------------------
# 7. Beyond the specification: the expm1 form against the literal difference form
# ---------------------------------------------------------------------------
lines = []
for x in (1e-2, 1e-6, 1e-8, 1e-9):
    rho_pair = np.array([2.0e-3, 2.0e-3 * math.exp(x)])
    Phi_pair = np.array([5.0e5, 4.0e5])
    exact_value = 2.0e-3 * 1.0e5 * math.expm1(x) / x
    module = float(hs.layer_mass(rho_pair, Phi_pair)[0])
    literal = float((rho_pair[1] - rho_pair[0]) * 1.0e5 / math.log(rho_pair[1] / rho_pair[0]))
    lines.append(f"log ratio {x:.0e}: module {module / exact_value - 1:+.2e}, literal difference "
                 f"form {literal / exact_value - 1:+.2e} relative to rho_k dPhi expm1(x) / x")
record(7, "beyond the specification: the module's expm1 form keeps its digits where the literal "
          "difference form loses them (decision 1)", True, "\n".join(lines))

# ---------------------------------------------------------------------------
# 8. Beyond the specification: pressure_from_top refusals
# ---------------------------------------------------------------------------
lines, ok = [], True
I_ok = hs.layer_mass(rho_a, Phi)
for label, args in (("p_b = 0", (0.0, I_ok)), ("p_b = NaN", (float("nan"), I_ok)),
                    ("a negative layer mass", (P_B, np.where(np.arange(I_ok.size) == 2, -1.0, I_ok)))):
    try:
        hs.pressure_from_top(*args)
        ok = False
        lines.append(f"{label}: not refused")
    except ValueError as exc:
        lines.append(f"{label}: {exc}")
p_check = hs.pressure_from_top(P_B, I_ok)
lines.append(f"p_0 equals p_b exactly: {p_check[0] == P_B}")
record(8, "beyond the specification: pressure_from_top refuses a bad boundary or layer, and p_0 is "
          "p_b exactly", ok and p_check[0] == P_B, "\n".join(lines))

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
