"""Acceptance for SPEC_13 v0.3 Step 2: the three runs.

Builds `forward/lindal_iris_v_ii_110` (the gradient from the derivative of the fitted temperatures)
and `forward/lindal_iris_v_sl_ii_110` (and the Sanchez-Lavega cloud wind) with `casspian-run-inputs`
and runs them with the transfer driver. `lindal_iris_ii_110` is read as SPEC_11 and SPEC_12 left it
and is rebuilt and rerun only if its product is absent. Run inputs and products are written in the
run directories, as a run's are; this script writes under `reports/step13_2/`.

**The relaxation.** The tools are uncommitted, so the runs' inputs carry `-dirty`; the one refusal
point, `casspian.lib.control._refuse_dirty_commit`, is replaced for this script's runs only.

**The round trip** is SPEC_11 Step 2's: the delivered `T(target) - T(anchor)` at 110, 290 and 730
mbar against the IRIS fit's change in value between the two planetographic latitudes. Its
difference is split as `tests/step11_2/diagnose_round_trip.py` splits it (REPORT_11_step2 section 3),
with the gradient the run's wind was built from in place of the slope:

Each integral runs from the anchor's latitude to the target's exactly, the integrand linear between
the wind grid's nodes. (The SPEC_11 diagnosis integrated over the nodes inside the path only, which
leaves out the partial cells at the two ends.)

* (A) the gradient's change, `T_anchor` times the integral of `d ln T` along the path, against the
  fit's change in value. It is given in two pieces that add to it: (A1) the plain integral of
  `dT/dphi` against the change in value, which the derivative makes zero to the differencing; and
  (A2) `T_anchor` times the integral of `d ln T` against the plain integral, the linearization;
* (B) inside the equatorial band, the wind's implied change against the gradient's;
* (C) outside the band, the same;
* (D) the delivered change against the wind's implied one, the transfer and its linearization.

Run from the repository root.
"""

import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.lib import geoid as gd
from casspian.lib.constants import AVOGADRO_CONSTANT, MOLAR_GAS_CONSTANT
from casspian.lib.gravity import g_eff_radial, omega_abs
from casspian.lib.windfield import WindField
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.run import run_inputs

HERE = Path("reports/step13_2")
HERE.mkdir(parents=True, exist_ok=True)
#: Each run's gradient and fit window (SPEC_13 v0.5 section 2a item 1: the two new runs use the smooth one).
RUNS = {"lindal_iris_ii_110": ("slope", "adaptive"), "lindal_iris_v_ii_110": ("value", "smooth"),
        "lindal_iris_v_sl_ii_110": ("value", "smooth")}
SHORT = {"lindal_iris_ii_110": "slope, Ingersoll and Pollard",
         "lindal_iris_v_ii_110": "value, smooth window, Ingersoll and Pollard",
         "lindal_iris_v_sl_ii_110": "value, smooth window, Sanchez-Lavega"}
#: Lindal's two labels, planetographic, where the IRIS fit's values are drawn (as SPEC_12's figure).
INGRESS_G, EGRESS_G = 36.3, -31.2
TABLE = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)
EGRESS = np.genfromtxt("occul_data/lindal/voyager2_egress.csv", delimiter=",", names=True)
EGRESS_AMU = 2.135
LEVELS = {110: 11000.0, 290: 29000.0, 730: 73000.0}
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def tree(path):
    handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


def at_pressure(p, values, targets):
    order = np.argsort(np.log(p))
    return np.interp(np.log(targets), np.log(p)[order], np.asarray(values, dtype="float64")[order])


warnings.simplefilter("ignore")
_original = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(attrs.get("casspian_git_commit", ""))
print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this script's runs only\n")
products, passes, rerun = {}, {}, {}
for name in RUNS:
    run = Path("forward") / name
    product_path = run / "output" / f"{name}_profile.nc"
    rerun[name] = name != "lindal_iris_ii_110" or not product_path.exists()
    if rerun[name]:
        run_inputs.build(run / f"{name}_build.toml")
        state = forward_transfer.run(run / f"{name}.toml")
        product_path = state.product
        passes[name] = int(state.state.record["passes"])
    products[name] = tree(product_path)
ctl._refuse_dirty_commit = _original

# ---------------------------------------------------------------------------
# 1. The round trip of each run
# ---------------------------------------------------------------------------
lines, summary = [], {}
for name, (gradient_rule, window) in RUNS.items():
    product = products[name]
    root = product.to_dataset(inherit=False)
    anchor_root = product["anchors/lindal"].to_dataset(inherit=False)
    thermo = product["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)
    wind = product["inputs/wind"].to_dataset(inherit=False)
    gravity = product["inputs/gravity"].to_dataset(inherit=False)
    rotation = product["inputs/rotation"].to_dataset(inherit=False)
    composition = product["inputs/composition"]
    p, T = np.asarray(root["pressure_Pa"].values), np.asarray(root["temperature_K"].values)
    p_a, T_a = np.asarray(thermo["pressure_Pa"].values), np.asarray(thermo["temperature_K"].values)
    identity = np.abs(T * np.asarray(root["pressure_identity_residual"].values)).max()
    lat_c = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    lat_g = np.asarray(wind["latitude_planetographic_deg"].values, dtype="float64")
    phi = np.radians(lat_c)
    dg_dc = np.gradient(lat_g, lat_c)
    anchor_c = float(anchor_root["latitude_planetocentric_deg"])
    anchor_g = float(anchor_root["latitude_planetographic_deg"])
    target_c = float(root["latitude_planetocentric_deg"])
    target_g = float(np.interp(target_c, lat_c, lat_g))
    Omega = float(rotation["angular_rate_rad_s"])
    GM, R_norm = float(gravity["GM_m3s2"]), float(gravity["normalization_radius_m"])
    J, degrees = np.asarray(gravity["J"].values), np.asarray(gravity["degree"].values)
    field = WindField(wind)
    edge = lw.EQUATORIAL_BAND_DEG
    on_path = [(target_c, anchor_c)]
    band = [(max(target_c, -edge), min(anchor_c, edge))]
    outside = [(target_c, max(target_c, -edge)), (min(anchor_c, edge), anchor_c)]

    def integral(values, pieces):
        """The change from north to south over latitude intervals (planetocentric deg), by the
        trapezoid on the grid's nodes inside each interval and on its two ends, linear between."""
        total = 0.0
        for lo, hi in pieces:
            if hi <= lo:
                continue
            inside = (lat_c > lo) & (lat_c < hi)
            x = np.concatenate([[lo], lat_c[inside], [hi]])
            y = np.interp(x, lat_c, values)
            total += float(np.trapezoid(y, np.radians(x)))
        return -total          # the path runs from north to south

    lines.append(f"{name} (gradient {gradient_rule!r}, window {window!r}): anchor {anchor_g:.2f} N ({anchor_c:.3f} planetocentric), "
                 f"target {abs(target_g):.2f} S ({abs(target_c):.3f} S planetocentric); "
                 f"{'rebuilt and rerun, ' + str(passes[name]) + ' passes' if rerun[name] else 'read as SPEC_12 left it'}; "
                 f"pressure identity at most {identity:.3f} K")
    summary[name] = {}
    for level, pressure in LEVELS.items():
        rows = TABLE["pressure_mbar"] == level
        f_grid = iris.fit(TABLE["latitude_planetographic_deg"][rows], TABLE["temperature_K"][rows], grid=lat_g,
                          window=window)
        f_ends = iris.fit(TABLE["latitude_planetographic_deg"][rows], TABLE["temperature_K"][rows],
                          grid=np.array([anchor_g, target_g]), window=window)
        iris_change = float(f_ends.value_K[1] - f_ends.value_K[0])
        delivered = float(at_pressure(p, T, [pressure])[0] - at_pressure(p_a, T_a, [pressure])[0])
        per_deg = (f_grid.gradient_K_per_deg if gradient_rule == "slope"
                   else lw.derivative_of_value(f_grid.value_K, lat_g))
        dT = per_deg * dg_dc * (180.0 / np.pi)            # K per planetocentric radian
        T_fit = f_grid.value_K
        u = field.wind_at(phi, np.full(phi.shape, pressure))
        _, s = field.wind_derivatives(phi, np.full(phi.shape, pressure))
        r = np.asarray(gd.wind_geoid(phi, R_norm, "equatorial_radius", lambda x: np.interp(x, phi, u),
                                     Omega, GM, J, degrees, R_norm).radius)
        g = g_eff_radial(u, r, phi, Omega, GM, J, degrees, R_norm)
        m_bar = lw.mean_molar_mass(composition, np.array([pressure]), lat_c)[0]
        d = -g * m_bar / (MOLAR_GAS_CONSTANT * T_fit) * (np.sin(phi) - np.cos(phi) / r * np.gradient(r, phi))
        d_ln_T_wind = -2.0 * omega_abs(u, r, phi, Omega) * r * (s * d + np.cos(phi) / r * np.gradient(u, phi)) / g
        T_anchor = float(np.interp(anchor_c, lat_c, T_fit))
        plain = integral(dT, on_path)
        gradient_change = T_anchor * integral(dT / T_fit, on_path)
        wind_out, grad_out = T_anchor * integral(d_ln_T_wind, outside), T_anchor * integral(dT / T_fit, outside)
        wind_band, grad_band = T_anchor * integral(d_ln_T_wind, band), T_anchor * integral(dT / T_fit, band)
        parts = {"A": gradient_change - iris_change, "A1": plain - iris_change, "A2": gradient_change - plain,
                 "B": wind_band - grad_band, "C": wind_out - grad_out,
                 "D": delivered - (wind_out + wind_band)}
        difference = delivered - iris_change
        summary[name][level] = {"delivered_K": delivered, "iris_change_K": iris_change, "difference_K": difference,
                                "iris_at_anchor_K": float(f_ends.value_K[0]), "iris_at_target_K": float(f_ends.value_K[1]),
                                **{f"part_{k}_K": v for k, v in parts.items()}}
        lines.append(f"  {level} mbar: delivered {delivered:+.2f} K, IRIS fit's change {iris_change:+.2f} K "
                     f"({f_ends.value_K[0]:.2f} to {f_ends.value_K[1]:.2f}), difference {difference:+.2f} K = "
                     f"(A) {parts['A']:+.2f} [(A1) {parts['A1']:+.3f}, (A2) {parts['A2']:+.3f}] + (B) {parts['B']:+.2f} "
                     f"+ (C) {parts['C']:+.2f} + (D) {parts['D']:+.2f}")
record(1, "the round trip of each run, as SPEC_11 Step 2's, with its difference split as REPORT_11_step2 "
          "section 3 splits it (reported)", True, "\n".join(lines))

# ---------------------------------------------------------------------------
# 2. For the author: the profiles against the egress, and the winds
# ---------------------------------------------------------------------------
p_e, T_e, s_e = EGRESS["pressure_mbar"] * 100.0, EGRESS["temperature_K"], EGRESS["temperature_uncertainty_K"]
rho_e = p_e * EGRESS_AMU * 1e-3 / (MOLAR_GAS_CONSTANT * T_e)
colours = {"lindal_iris_ii_110": "0.35", "lindal_iris_v_ii_110": "C0", "lindal_iris_v_sl_ii_110": "C2"}
fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
axes[0].fill_betweenx(p_e / 100, T_e - s_e, T_e + s_e, color="C3", alpha=0.18, linewidth=0, label="egress +- sigma_T")
axes[0].plot(T_e, p_e / 100, color="C3", linewidth=1.0, label="Voyager 2 egress, 31.2 S (26.06 S planetocentric)")
axes[1].fill_betweenx(p_e / 100, -s_e, s_e, color="C3", alpha=0.18, linewidth=0, label="egress sigma_T")
axes[2].plot(rho_e, p_e / 100, color="C3", linewidth=1.0, label="egress, 2.135 amu")
level_lines = []
for name, product in products.items():
    root = product.to_dataset(inherit=False)
    p, T = np.asarray(root["pressure_Pa"].values), np.asarray(root["temperature_K"].values)
    rho = (np.asarray(root["number_density_m3"].values) * np.asarray(root["mean_molar_mass_kg_mol"].values)
           / AVOGADRO_CONSTANT)
    within = (p_e >= p.min()) & (p_e <= p.max())
    axes[0].plot(T, p / 100, color=colours[name], linewidth=1.2, label=f"{name} ({SHORT[name]})")
    axes[1].plot(at_pressure(p, T, p_e[within]) - T_e[within], p_e[within] / 100, color=colours[name],
                 linewidth=1.2, label=f"{name} - egress")
    axes[2].plot(rho, p / 100, color=colours[name], linewidth=1.2, label=name)
    values = []
    for level in LEVELS:
        pa = level * 100.0
        dT = float(at_pressure(p, T, [pa])[0] - at_pressure(p_e, T_e, [pa])[0])
        dr = float(at_pressure(p, rho, [pa])[0] / at_pressure(p_e, rho_e, [pa])[0] - 1.0)
        values.append(f"{level} mbar {dT:+.2f} K ({100 * dr:+.1f} % density)")
        summary[name][level].update({"minus_egress_K": dT, "density_relative_to_egress": dr})
    level_lines.append(f"{name} - egress: " + "; ".join(values))
sigma_lines = "egress sigma_T: " + ", ".join(f"{level} mbar {float(at_pressure(p_e, s_e, [level * 100.0])[0]):.2f} K"
                                             for level in LEVELS)
# The IRIS fit's values at Lindal's two labels, as SPEC_12's figure draws them, with each window.
levels_Pa = np.array(list(LEVELS.values()))
iris_lines = []
egress_at = at_pressure(p_e, T_e, levels_Pa)
for window, face in (("adaptive", "C0"), ("smooth", "C9")):
    at = {}
    for level in LEVELS:
        rows = TABLE["pressure_mbar"] == level
        at[level] = iris.fit(TABLE["latitude_planetographic_deg"][rows], TABLE["temperature_K"][rows],
                             grid=np.array([INGRESS_G, EGRESS_G]), window=window).value_K
    north = np.array([at[v][0] for v in LEVELS])
    south = np.array([at[v][1] for v in LEVELS])
    axes[0].plot(north, levels_Pa / 100, "o", markerfacecolor="none", markeredgecolor=face, markersize=8,
                 label=f"IRIS fit ({window} window), 36.3 N (30.81 N planetocentric)")
    axes[0].plot(south, levels_Pa / 100, "s", color=face, markersize=7,
                 label=f"IRIS fit ({window} window), 31.2 S (26.06 S planetocentric)")
    axes[1].plot(south - egress_at, levels_Pa / 100, "s", color=face, markersize=7,
                 label=f"IRIS fit ({window} window) at 31.2 S - egress")
    iris_lines.append(f"IRIS fit, {window} window: " + "; ".join(
        f"{level} mbar {at[level][0]:.2f} K at 36.3 N and {at[level][1]:.2f} K at 31.2 S, change "
        f"{at[level][1] - at[level][0]:+.2f} K, 31.2 S - egress {at[level][1] - egress_at[k]:+.2f} K"
        for k, level in enumerate(LEVELS)))
axes[1].axvline(0, color="0.6", linewidth=0.8)
axes[0].set_xlabel("temperature (K)")
axes[0].set_ylabel("pressure (mbar)")
axes[1].set_xlabel("temperature minus the egress (K)")
axes[2].set_xscale("log")
axes[2].set_xlabel("mass density (kg/m3)")
for ax in axes:
    ax.set_yscale("log")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
axes[0].invert_yaxis()
fig.suptitle("The three runs at 31.2 S (26.06 S planetocentric) against the Voyager 2 egress", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "profiles.png", dpi=130)
plt.close(fig)

winds = {name: product["inputs/wind"].to_dataset(inherit=False) for name, product in products.items()}
names = list(RUNS)
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
wind_lines = []
for k, name in enumerate(names):
    w = winds[name]
    g = np.asarray(w["latitude_planetographic_deg"].values)
    pw = np.asarray(w["pressure_Pa"].values) / 100.0
    mesh = axes[0, k].pcolormesh(g, pw, np.asarray(w["u_total_ms"].values).T, cmap="RdBu_r", vmin=-500, vmax=500,
                                 shading="nearest")
    axes[0, k].set_title(f"u_total, {name}\n({SHORT[name]})", fontsize=9)
fig.colorbar(mesh, ax=axes[0, :], label="u_total (m/s)", shrink=0.8)
for k, (before, after) in enumerate(zip(names[:-1], names[1:])):
    a, b = winds[before], winds[after]
    g = np.asarray(b["latitude_planetographic_deg"].values)
    pw = np.asarray(b["pressure_Pa"].values) / 100.0
    du = np.asarray(b["u_total_ms"].values) - np.asarray(a["u_total_ms"].values)
    ds = np.asarray(b["u_shear_ms"].values) - np.asarray(a["u_shear_ms"].values)
    dr = np.asarray(b["u_reference_ms"].values) - np.asarray(a["u_reference_ms"].values)
    limit = float(np.nanpercentile(np.abs(du), 99.5))
    m = axes[1, k].pcolormesh(g, pw, du.T, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest")
    fig.colorbar(m, ax=axes[1, k], label="m/s")
    axes[1, k].set_title(f"u_total: {after} - {before}", fontsize=9)
    worst = np.unravel_index(np.nanargmax(np.abs(ds)), ds.shape)
    band = np.abs(np.asarray(b["latitude_planetocentric_deg"].values)) < 15.0
    wind_lines.append(
        f"{after} - {before}: u_total largest {du.flat[np.nanargmax(np.abs(du))]:+.1f} m/s, RMS "
        f"{np.sqrt(np.nanmean(du ** 2)):.1f}; reference RMS {np.sqrt(np.nanmean(dr ** 2)):.1f} m/s; shear "
        f"(u_shear) largest {ds[worst]:+.1f} m/s at {g[worst[0]]:.1f} deg and {pw[worst[1]]:.4g} mbar, RMS within 15 deg "
        f"of the equator {np.sqrt(np.nanmean(ds[band] ** 2)):.1f} and beyond {np.sqrt(np.nanmean(ds[~band] ** 2)):.1f} m/s")
w3, w2 = winds[names[2]], winds[names[1]]
g = np.asarray(w3["latitude_planetographic_deg"].values)
pw = np.asarray(w3["pressure_Pa"].values) / 100.0
ds = np.asarray(w3["u_shear_ms"].values) - np.asarray(w2["u_shear_ms"].values)
limit = float(np.nanpercentile(np.abs(ds), 99.5))
m = axes[1, 2].pcolormesh(g, pw, ds.T, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest")
fig.colorbar(m, ax=axes[1, 2], label="m/s")
axes[1, 2].set_title(f"u_shear: {names[2]} - {names[1]}\n(the shear's change with the cloud wind, section 3 ruling 3)",
                     fontsize=9)
for ax in axes.flat:
    ax.set_yscale("log")
    ax.set_ylim(1e4, 1e-2)
    ax.axhline(398.107, color="0.4", linewidth=0.6, linestyle=":")
for ax in axes[1]:
    ax.set_xlabel("planetographic latitude (deg)")
for ax in axes[:, 0]:
    ax.set_ylabel("pressure (mbar)")
fig.suptitle("The three runs' winds, and the difference of each from the one before (dotted: 398 mbar, the "
             "reference level)", fontsize=10)
fig.savefig(HERE / "winds.png", dpi=110)
plt.close(fig)
record(2, "for the author: the delivered temperature and density of the three runs with the egress and its "
          "sigma_T and the IRIS fit's values at 36.3 N and 31.2 S; u_total of the three winds and the difference of each from the one before", True,
       "\n".join(level_lines) + "\n" + sigma_lines + "\n" + "\n".join(iris_lines) + "\n" + "\n".join(wind_lines)
       + f"\nfigures {HERE / 'profiles.png'}, {HERE / 'winds.png'}")

(HERE / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
