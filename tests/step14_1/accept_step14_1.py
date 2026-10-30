"""Acceptance for SPEC_14 v0.2 Step 1: hypothesis 2, the shear reversed aloft.

Builds `forward/lindal_iris_v_sl_ii_110` (the base, case ii) and `forward/lindal_iris_v_sl_iii_110`
(case iii, `above_top = "reversed"`, tuned to the Voyager 2 egress) with `casspian-run-inputs` and
runs both with the transfer driver, so that the two are compared on the same code. Run inputs and
products are written in the run directories, as a run's are; this script writes under
`reports/step14_1/`.

Latitudes are given planetographic, with planetocentric in parentheses. The target is Lindal's
egress label, 31.2 S (26.06 S planetocentric); the anchor 36.3 N (30.81 N planetocentric).

Run from the repository root.
"""

import json
import sys
import tomllib
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from casspian.forward import transfer as forward_transfer
from casspian.lib.constants import AVOGADRO_CONSTANT, MOLAR_GAS_CONSTANT
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.run import run_inputs
from casspian.tools.wind import shear as shear_tool

HERE = Path("reports/step14_1")
HERE.mkdir(parents=True, exist_ok=True)
BASE, NEW = "lindal_iris_v_sl_ii_110", "lindal_iris_v_sl_iii_110"
LABEL = {BASE: "case (ii), relaxed above 110 mbar", NEW: "case (iii), reversed aloft, tuned to the egress"}
COLOUR = {BASE: "C2", NEW: "C4"}
INGRESS_G, EGRESS_G = 36.3, -31.2
TABLE = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)
EGRESS = np.genfromtxt("occul_data/lindal/voyager2_egress.csv", delimiter=",", names=True)
EGRESS_AMU = 2.135
IRIS_LEVELS = {110: 11000.0, 290: 29000.0, 730: 73000.0}
TABLE_LEVELS_MBAR = (730.0, 290.0, 110.0, 10.0, 4.5, 1.0, 0.3)
EQUATOR_LEVELS_MBAR = (10.0, 1.0, 0.2)
TOP_Pa = 11000.0
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


def column(wind, latitude_c):
    """`u_total(p)` at a planetocentric latitude, linear in latitude between the grid's nodes."""
    lat = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    u = np.asarray(wind["u_total_ms"].values, dtype="float64")
    return np.array([np.interp(latitude_c, lat, u[:, j]) for j in range(u.shape[1])])


warnings.simplefilter("ignore")
products, passes = {}, {}
for name in (BASE, NEW):
    run = Path("forward") / name
    run_inputs.build(run / f"{name}_build.toml")
    result = forward_transfer.run(run / f"{name}.toml")
    passes[name] = int(result.state.record["passes"])
    products[name] = tree(result.product)
roots = {n: p.to_dataset(inherit=False) for n, p in products.items()}
winds = {n: p["inputs/wind"].to_dataset(inherit=False) for n, p in products.items()}
summary = {"passes": passes}

# ---------------------------------------------------------------------------
# 1. Below the top level the wind is the base's, and so is the shear at the three levels
# ---------------------------------------------------------------------------
pw = np.asarray(winds[NEW]["pressure_Pa"].values, dtype="float64")
below = pw >= TOP_Pa
u_base, u_new = (np.asarray(winds[n]["u_total_ms"].values) for n in (BASE, NEW))
same_grid = np.array_equal(pw, np.asarray(winds[BASE]["pressure_Pa"].values))
wind_equal = same_grid and np.array_equal(u_base[:, below], u_new[:, below])
others = {v: np.array_equal(np.asarray(winds[BASE][v].values), np.asarray(winds[NEW][v].values), equal_nan=True)
          for v in ("u_reference_ms", "u_total_uncertainty_ms", "latitude_planetocentric_deg")}
largest_below = float(np.max(np.abs(u_new[:, below] - u_base[:, below]))) if same_grid else float("nan")

# The shear at the three levels, from the construction itself with each run's own inputs.
shears = {}
for name in (BASE, NEW):
    run = Path("forward") / name
    with open(run / f"{name}_build.toml", "rb") as handle:
        section = tomllib.load(handle)["shear"]
    inputs = {key: shear_tool._read(run / section[key], kind) for key, kind in shear_tool.INPUT_KEYS.items()}
    source = shear_tool._read(run / section["source"], "wind")
    built = lw.construct(source, inputs["temperatures"], inputs["composition"], inputs["gravity"],
                         inputs["rotation"], float(section["top_level_Pa"]), section["above_top"],
                         gradient=section.get("gradient", "slope"), window=section.get("window", "adaptive"),
                         hold_Pa=section.get("hold_Pa"), cross_Pa=section.get("cross_Pa"))
    shears[name] = built
shear_equal = np.array_equal(shears[BASE].shear_ms, shears[NEW].shear_ms)
iterations_equal = shears[BASE].iterations == shears[NEW].iterations

p_b, T_b = (np.asarray(roots[BASE][v].values, dtype="float64") for v in ("pressure_Pa", "temperature_K"))
p_n, T_n = (np.asarray(roots[NEW][v].values, dtype="float64") for v in ("pressure_Pa", "temperature_K"))
deep = p_b >= TOP_Pa
dT_deep = at_pressure(p_n, T_n, p_b[deep]) - T_b[deep]
worst = int(np.argmax(np.abs(dT_deep)))
largest_deep = float(dT_deep[worst])
summary["largest_change_at_and_below_110mbar_K"] = largest_deep
# SPEC_14 v0.3 (author, 30 October): equal to within the construction's tolerance, TOLERANCE_MS. The
# stopping rule takes the largest change anywhere on the grid, so case (iii) can stop on a later
# iterate than the base (REPORT_14_step1 section 3, diagnose_below_top.py).
shear_largest = float(np.max(np.abs(shears[NEW].shear_ms - shears[BASE].shear_ms)))
record(1, "below the top level the wind is the base's, and so is the shear at the three levels, to within "
          "the construction's tolerance; the delivered temperature at and below 110 mbar changes by at "
          "most 0.1 K, or the change is explained",
       same_grid and largest_below < lw.TOLERANCE_MS and shear_largest < lw.TOLERANCE_MS
       and all(others.values()) and abs(largest_deep) <= 0.1,
       f"u_total at the {int(below.sum())} pressure nodes at and below 110 mbar, every latitude: array-equal "
       f"{wind_equal}, largest difference {largest_below:.3e} m/s; "
       f"{', '.join(f'{k} equal {v}' for k, v in others.items())}\n"
       f"the shear at 110, 290 and 730 mbar from the construction: array-equal {shear_equal}, largest "
       f"difference {shear_largest:.3e} m/s, against the tolerance {lw.TOLERANCE_MS:g} m/s "
       f"(iterations {shears[BASE].iterations} and {shears[NEW].iterations}, equal {iterations_equal}; "
       f"see diagnose_below_top.py)\n"
       f"the delivered temperature at and below 110 mbar, new minus base: largest {largest_deep:+.4f} K at "
       f"{p_b[deep][worst] / 100:.4g} mbar (on the base's levels); outer loop passes {passes[BASE]} and "
       f"{passes[NEW]}")

# ---------------------------------------------------------------------------
# 2. The pressure identity, in kelvins, level by level
# ---------------------------------------------------------------------------
identity = T_n * np.asarray(roots[NEW]["pressure_identity_residual"].values, dtype="float64")
identity_base = T_b * np.asarray(roots[BASE]["pressure_identity_residual"].values, dtype="float64")
with open(HERE / "pressure_identity.csv", "w", encoding="utf-8", newline="\n") as handle:
    handle.write("pressure_mbar,temperature_K,identity_K\n")
    for p, T, e in zip(p_n, T_n, identity):
        handle.write(f"{p / 100:.6g},{T:.4f},{e:+.5f}\n")
k = int(np.argmax(np.abs(identity)))
summary["pressure_identity_largest_K"] = float(identity[k])
record(2, "the pressure identity of the run, in kelvins, level by level (reported)", True,
       f"{NEW}: largest |T x residual| {abs(identity[k]):.4f} K at {p_n[k] / 100:.4g} mbar, RMS "
       f"{np.sqrt(np.mean(identity ** 2)):.4f} K over {identity.size} levels; base {np.abs(identity_base).max():.4f} K\n"
       f"every level in {HERE / 'pressure_identity.csv'}")

# ---------------------------------------------------------------------------
# 3. For the author
# ---------------------------------------------------------------------------
anchor = products[NEW]["anchors/lindal"].to_dataset(inherit=False)
anchor_c = float(anchor["latitude_planetocentric_deg"])
target_c = float(roots[NEW]["latitude_planetocentric_deg"])
thermo = products[NEW]["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)
p_a, T_a = np.asarray(thermo["pressure_Pa"].values), np.asarray(thermo["temperature_K"].values)
p_e, T_e, s_e = EGRESS["pressure_mbar"] * 100.0, EGRESS["temperature_K"], EGRESS["temperature_uncertainty_K"]
rho_e = p_e * EGRESS_AMU * 1e-3 / (MOLAR_GAS_CONSTANT * T_e)
lines = []

# (a) f(p) and the shear at the anchor's and the target's latitudes
grid = pw[pw <= TOP_Pa]
f_new = lw.reversed_factor(grid, 1500.0, 450.0)
f_base = np.exp(-np.log(TOP_Pa / grid) / lw.ABOVE_SCALE_HEIGHTS)
fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True)
axes[0].plot(f_base, grid / 100, color=COLOUR[BASE], label=f"{BASE}: exp(-x / 2)")
axes[0].plot(f_new, grid / 100, color=COLOUR[NEW], label=f"{NEW}: held, then reversed")
axes[0].axvline(0, color="0.6", linewidth=0.8)
for p_mark, text in ((1500.0, "hold_Pa"), (450.0, "cross_Pa"), (20.0, "0.2 mbar")):
    axes[0].axhline(p_mark / 100, color="0.5", linewidth=0.6, linestyle=":")
    axes[0].text(-2.5, p_mark / 100, f" {text}", fontsize=7, va="bottom")
axes[0].set_xlabel("f(p), the factor on the 110 mbar shear")
for ax, lat_c, where in ((axes[1], anchor_c, "anchor, 36.3 N (30.81 N planetocentric)"),
                         (axes[2], target_c, "target, 31.2 S (26.06 S planetocentric)")):
    for name in (BASE, NEW):
        u = column(winds[name], lat_c)
        s = np.gradient(u, np.log(pw))
        ax.plot(s, pw / 100, color=COLOUR[name], label=f"{name}")
    ax.axvline(0, color="0.6", linewidth=0.8)
    ax.set_title(f"du/dln p at the {where}", fontsize=9)
    ax.set_xlabel("shear du/dln p (m/s)")
for ax in axes:
    ax.set_yscale("log")
    ax.set_ylim(1e3, 1e-2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
axes[0].set_ylabel("pressure (mbar)")
fig.suptitle("Case (iii) against case (ii): f(p) and the shear along the two columns", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "shear.png", dpi=120)
plt.close(fig)

# (b) the profiles at the target, with the anchor, the egress and the IRIS values
fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
axes[0].fill_betweenx(p_e / 100, T_e - s_e, T_e + s_e, color="C3", alpha=0.18, linewidth=0, label="egress +- sigma_T")
axes[0].plot(T_e, p_e / 100, color="C3", linewidth=1.0, label="Voyager 2 egress, 31.2 S (26.06 S planetocentric)")
axes[0].plot(T_a, p_a / 100, color="0.3", linewidth=1.0, linestyle="--",
             label="anchor, Lindal ingress, 36.3 N (30.81 N planetocentric)")
axes[1].fill_betweenx(p_e / 100, -s_e, s_e, color="C3", alpha=0.18, linewidth=0, label="egress sigma_T")
axes[2].plot(rho_e, p_e / 100, color="C3", linewidth=1.0, label="egress, 2.135 amu")
table = {}
for name in (BASE, NEW):
    root = roots[name]
    p, T = np.asarray(root["pressure_Pa"].values), np.asarray(root["temperature_K"].values)
    rho = (np.asarray(root["number_density_m3"].values) * np.asarray(root["mean_molar_mass_kg_mol"].values)
           / AVOGADRO_CONSTANT)
    within = (p_e >= p.min()) & (p_e <= p.max())
    axes[0].plot(T, p / 100, color=COLOUR[name], linewidth=1.2, label=f"{name} ({LABEL[name]})")
    axes[1].plot(at_pressure(p, T, p_e[within]) - T_e[within], p_e[within] / 100, color=COLOUR[name],
                 linewidth=1.2, label=f"{name} - egress")
    axes[2].plot(rho, p / 100, color=COLOUR[name], linewidth=1.2, label=name)
    table[name] = {}
    for level in TABLE_LEVELS_MBAR:
        pa = level * 100.0
        table[name][level] = {"minus_egress_K": float(at_pressure(p, T, [pa])[0] - at_pressure(p_e, T_e, [pa])[0]),
                              "density_relative": float(at_pressure(p, rho, [pa])[0] / at_pressure(p_e, rho_e, [pa])[0] - 1.0)}
levels_Pa = np.array(list(IRIS_LEVELS.values()))
egress_at = at_pressure(p_e, T_e, levels_Pa)
at = {}
for level in IRIS_LEVELS:
    rows = TABLE["pressure_mbar"] == level
    at[level] = iris.fit(TABLE["latitude_planetographic_deg"][rows], TABLE["temperature_K"][rows],
                         grid=np.array([INGRESS_G, EGRESS_G]), window="smooth").value_K
north, south = np.array([at[v][0] for v in IRIS_LEVELS]), np.array([at[v][1] for v in IRIS_LEVELS])
axes[0].plot(north, levels_Pa / 100, "o", markerfacecolor="none", markeredgecolor="C9", markersize=8,
             label="IRIS fit (smooth window), 36.3 N (30.81 N planetocentric)")
axes[0].plot(south, levels_Pa / 100, "s", color="C9", markersize=7,
             label="IRIS fit (smooth window), 31.2 S (26.06 S planetocentric)")
axes[1].plot(south - egress_at, levels_Pa / 100, "s", color="C9", markersize=7, label="IRIS fit at 31.2 S - egress")
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
fig.suptitle("The base and case (iii) at 31.2 S (26.06 S planetocentric) against the Voyager 2 egress; "
             "case (iii) is tuned to the egress", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "profiles.png", dpi=130)
plt.close(fig)
sigma_at = {level: float(at_pressure(p_e, s_e, [level * 100.0])[0]) for level in TABLE_LEVELS_MBAR}
lines.append("delivered minus egress at 31.2 S (26.06 S planetocentric), K (density relative to the egress):")
lines.append("  level (mbar)  " + "  ".join(f"{level:>8g}" for level in TABLE_LEVELS_MBAR))
for name in (BASE, NEW):
    lines.append(f"  {name:24s}" + "  ".join(f"{table[name][l]['minus_egress_K']:+8.2f}" for l in TABLE_LEVELS_MBAR))
    lines.append(f"  {'':24s}" + "  ".join(f"{100 * table[name][l]['density_relative']:+7.1f}%" for l in TABLE_LEVELS_MBAR))
lines.append(f"  {'egress sigma_T':24s}" + "  ".join(f"{sigma_at[l]:8.2f}" for l in TABLE_LEVELS_MBAR))
lines.append("IRIS fit (smooth window) at 31.2 S minus the egress: " + "; ".join(
    f"{level} mbar {at[level][1] - egress_at[k]:+.2f} K" for k, level in enumerate(IRIS_LEVELS)))
summary["minus_egress"] = table

# (c) u_total of the new wind and its difference from the base
g = np.asarray(winds[NEW]["latitude_planetographic_deg"].values)
lat_c = np.asarray(winds[NEW]["latitude_planetocentric_deg"].values)
du = u_new - u_base
fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
m0 = axes[0].pcolormesh(g, pw / 100, u_new.T, cmap="RdBu_r", vmin=-1000, vmax=1000, shading="nearest")
fig.colorbar(m0, ax=axes[0], label="u_total (m/s)")
axes[0].set_title(f"u_total, {NEW}", fontsize=9)
limit = float(np.nanpercentile(np.abs(du), 99.5))
m1 = axes[1].pcolormesh(g, pw / 100, du.T, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest")
fig.colorbar(m1, ax=axes[1], label="m/s")
axes[1].set_title(f"u_total: {NEW} - {BASE}", fontsize=9)
for ax in axes:
    ax.set_yscale("log")
    ax.set_ylim(1e4, 1e-2)
    ax.axhline(TOP_Pa / 100, color="0.4", linewidth=0.6, linestyle=":")
    ax.axhline(0.2, color="0.4", linewidth=0.6, linestyle=":")
    ax.set_xlabel("planetographic latitude (deg)")
axes[0].set_ylabel("pressure (mbar)")
fig.suptitle("The case (iii) wind and its difference from the base (dotted: 110 and 0.2 mbar)", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "winds.png", dpi=110)
plt.close(fig)
aloft = pw < TOP_Pa
with open(HERE / "largest_wind_aloft.csv", "w", encoding="utf-8", newline="\n") as handle:
    handle.write("pressure_mbar,largest_abs_u_base_ms,largest_abs_u_new_ms,latitude_g_new_deg\n")
    for j in np.flatnonzero(aloft)[::-1]:
        i = int(np.argmax(np.abs(u_new[:, j])))
        handle.write(f"{pw[j] / 100:.6g},{np.abs(u_base[:, j]).max():.1f},{abs(u_new[i, j]):.1f},{g[i]:.1f}\n")
picks = [11000.0, 1000.0, 450.0, 100.0, 20.0, 1.0]
for name, u in ((BASE, u_base), (NEW, u_new)):
    parts = []
    for p_pick in picks:
        j = int(np.argmin(np.abs(np.log(pw / p_pick))))
        i = int(np.argmax(np.abs(u[:, j])))
        parts.append(f"{pw[j] / 100:.4g} mbar {u[i, j]:+.0f} m/s at {g[i]:.1f}")
    lines.append(f"largest |u| aloft, {name}: " + "; ".join(parts))
equator = {}
for name, u in ((BASE, u_base), (NEW, u_new)):
    u0 = np.array([np.interp(0.0, lat_c, u[:, j]) for j in range(pw.size)])
    equator[name] = {level: float(at_pressure(pw, u0, [level * 100.0])[0]) for level in EQUATOR_LEVELS_MBAR}
    lines.append(f"equatorial wind, {name}: " + "; ".join(f"{l:g} mbar {v:+.0f} m/s" for l, v in equator[name].items()))
lines.append(f"every wind grid level above 110 mbar in {HERE / 'largest_wind_aloft.csv'}")
summary["equatorial_wind_ms"] = equator

# (d) the delivered lapse rate above 110 mbar
fig, ax = plt.subplots(figsize=(7, 7))
for name in (BASE, NEW):
    root = roots[name]
    p = np.asarray(root["pressure_Pa"].values)
    T = np.asarray(root["temperature_K"].values)
    z = np.asarray(root["altitude_m"].values)
    gamma = -np.gradient(T, z) * 1000.0
    up = p <= TOP_Pa
    ax.plot(gamma[up], p[up] / 100, color=COLOUR[name], label=name)
    picks_l = [at_pressure(p, gamma, [level * 100.0])[0] for level in (110.0, 10.0, 4.5, 1.0, 0.3)]
    lines.append(f"lapse rate -dT/dz, {name}: " + "; ".join(
        f"{level:g} mbar {v:+.2f} K/km" for level, v in zip((110.0, 10.0, 4.5, 1.0, 0.3), picks_l))
        + f"; range above 110 mbar {gamma[up].min():+.2f} to {gamma[up].max():+.2f} K/km")
ax.axvline(0, color="0.6", linewidth=0.8)
ax.set_yscale("log")
ax.invert_yaxis()
ax.set_xlabel("lapse rate -dT/dz (K/km)")
ax.set_ylabel("pressure (mbar)")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
ax.set_title("The delivered lapse rate above 110 mbar at 31.2 S (26.06 S planetocentric)", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "lapse_rate.png", dpi=120)
plt.close(fig)

record(3, "for the author: f(p) and the shear at the two latitudes; the profiles at 31.2 S against the "
          "egress; the table of delivered minus egress; the new wind and its difference; the lapse rate "
          "above 110 mbar (reported)", True,
       "\n".join(lines) + f"\nfigures {HERE / 'shear.png'}, {HERE / 'profiles.png'}, {HERE / 'winds.png'}, "
       f"{HERE / 'lapse_rate.png'}")

(HERE / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8", newline="\n")
sys.exit(1 if failed else 0)
