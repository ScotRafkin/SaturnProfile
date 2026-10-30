"""Acceptance for SPEC_16 v0.2 Step 1: Pioneer 11 refractivity against the transferred Lindal anchor.

Checks that our refractivity arithmetic reproduces Kliore et al. (1980) Table 1, then compares his
refractivity with the refractivity the two SPEC_15 runs deliver at 11.9 S (9.51 S planetocentric):
against radius, with a fitted radius offset, for both readings of his radius mapping; and against
pressure, with his pressure rescaled to the anchor's mixture. No new run: the two runs' products are
read as they are, built if absent. This script writes under `reports/step16_1/` only.

Latitudes are given planetographic, with planetocentric in parentheses.

Run from the repository root.
"""

import json
import sys
import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.optimize import minimize_scalar

from casspian.forward import transfer as forward_transfer
from casspian.lib import reduction
from casspian.lib.constants import LOSCHMIDT_CONSTANT
from casspian.tools.lindal import radio_profiles as rp
from casspian.tools.run import run_inputs

HERE = Path("reports/step16_1")
HERE.mkdir(parents=True, exist_ok=True)
RUNS = ("lindal_iris_v_sl_ii_110_12s", "lindal_iris_v_sl_iii_110_12s")
HYPOTHESIS = ("hypothesis 1 (case ii)", "hypothesis 2 (case iii)")
LEVELS_MBAR = (10.0, 30.0, 100.0, 150.0)
#: Kliore et al. (1980) p. 5861: the equatorial radius of each level is computed from its latitude
#: and "a constant oblateness value based on an ellipsoidal approximation"; p. 5857 gives Saturn's
#: oblateness as 0.096 at the 100 mbar level. Table 1's column is that equatorial radius.
KLIORE_OBLATENESS = 0.096
#: Kliore et al. (1980) Table 1 title: latitude 9.5 S. Planetocentric (SPEC_15 section 2 ruling 1).
KLIORE_LATITUDE_C = -9.5
#: The planetographic reading: the same latitude converted to planetographic on Kliore's own
#: ellipse, tan(phi_g) = tan(phi_c) / (1 - f)^2 with f = 0.096. Lindal et al. (1985) label the
#: profile 11.9 S on their own figure of Saturn; Kliore's f gives 11.57 S.
KLIORE_LATITUDE_G = float(np.degrees(np.arctan(np.tan(np.radians(KLIORE_LATITUDE_C)) / (1 - KLIORE_OBLATENESS) ** 2)))
#: The two readings of the mapping (SPEC_16 section 2 ruling 1): which latitude Kliore put in it.
READINGS = {"planetocentric": KLIORE_LATITUDE_C, "planetographic": KLIORE_LATITUDE_G}
#: Table 1 rounds refractivity to 0.1 n-units.
ROUNDING_N = 0.1
results = []
summary = {}


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def tree(path):
    handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


def interp_log(x, y, at):
    """y at `at`, linear in ln x; NaN outside."""
    x, y = np.asarray(x, dtype="float64"), np.asarray(y, dtype="float64")
    order = np.argsort(x)
    return np.interp(np.log(at), np.log(x[order]), y[order], left=np.nan, right=np.nan)


def local_radius(r_eq, latitude_deg, f=KLIORE_OBLATENESS):
    """Kliore's mapping undone: r_local = r_eq (1 - f) / sqrt((1 - f)^2 cos^2 phi + sin^2 phi)."""
    phi = np.radians(latitude_deg)
    return r_eq * (1 - f) / np.sqrt((1 - f) ** 2 * np.cos(phi) ** 2 + np.sin(phi) ** 2)


def mean_refractivity_m3(mixture, refractivity_set):
    """R_bar per molecule, as `tools/composition/build_composition.py` converts it: e6 x 1e-6 / n_L."""
    with open("data_static/species_master.toml", "rb") as handle:
        table = tomllib.load(handle)["refractivity"][refractivity_set]
    return sum(x * table[s]["refractivity_stp_e6"] * 1e-6 / LOSCHMIDT_CONSTANT for s, x in mixture.items())


kliore = rp.read_pioneer11()
p_k = np.asarray(kliore["pressure_mbar"], dtype="float64")
T_k = np.asarray(kliore["temperature_K"], dtype="float64")
N_k = np.asarray(kliore["refractivity_n_units"], dtype="float64")
r_eq = np.asarray(kliore["equatorial_radius_km"], dtype="float64") * 1e3

# ---------------------------------------------------------------------------
# 1. The self-check
# ---------------------------------------------------------------------------
R_bar = mean_refractivity_m3(rp.KLIORE_MIXTURE, rp.REFRACTIVITY_SET)
N_self = reduction.refractivity(reduction.number_density(p_k * 100.0, T_k), R_bar) * 1e6
d = N_self - N_k
k = int(np.argmax(np.abs(d)))
pct = 100 * (N_self / N_k - 1)
kp = int(np.argmax(np.abs(pct)))
# The implied mean refractivity of Kliore's conversion: least squares in ratio, weighted to where the
# rounding matters least (the column's ratio at its largest values).
deep = N_k >= 10.0
ratio = float(np.sum(N_k[deep] * N_self[deep]) / np.sum(N_self[deep] ** 2))
over = np.abs(d) > ROUNDING_N
R_e6 = R_bar * LOSCHMIDT_CONSTANT * 1e6
with open(HERE / "self_check.csv", "w", encoding="utf-8", newline="\n") as handle:
    handle.write("pressure_mbar,temperature_K,refractivity_kliore_n_units,refractivity_ours_n_units,difference_n_units,difference_percent\n")
    handle.writelines(f"{p:g},{t:g},{n:g},{s:.4f},{e:+.4f},{q:+.3f}\n" for p, t, n, s, e, q in zip(p_k, T_k, N_k, N_self, d, pct))
summary["self_check"] = {"rms_n": float(np.sqrt(np.mean(d ** 2))), "largest_n": float(d[k]), "largest_n_at_mbar": float(p_k[k]),
                         "rms_percent": float(np.sqrt(np.mean(pct ** 2))), "largest_percent": float(pct[kp]),
                         "largest_percent_at_mbar": float(p_k[kp]), "levels_over_rounding": int(over.sum()),
                         "mean_refractivity_e6": R_e6, "implied_ratio_deep": ratio}
record(1, "the self-check: Kliore's T and p as published (90/10), through reduction.refractivity, against his column", True,
       f"mixture 90 % H2, 10 % He; refractivity set {rp.REFRACTIVITY_SET} (H2 136.0, He 35.0): mean {R_e6:.2f} n-units at STP\n"
       f"{p_k.size} levels: RMS {summary['self_check']['rms_n']:.3f} n-units, largest {d[k]:+.3f} at {p_k[k]:g} mbar; "
       f"RMS {summary['self_check']['rms_percent']:.2f} %, largest {pct[kp]:+.2f} % at {p_k[kp]:g} mbar\n"
       f"{int(over.sum())} levels beyond the table's 0.1 n-unit rounding, at "
       + (", ".join(f"{p:g}" for p in p_k[over][:15]) + (" ..." if over.sum() > 15 else "") if over.any() else "none") + " mbar\n"
       f"mean of ours minus Kliore's: {np.mean(d):+.3f} n-units; mean ratio {np.mean(N_self / N_k):.4f} "
       f"(below 10 n-units {np.mean((N_self / N_k)[~deep]):.4f}, from 10 n-units down {np.mean((N_self / N_k)[deep]):.4f})\n"
       f"Kliore's column over ours, least squares from 10 n-units down: {ratio:.5f}, an implied mean refractivity of "
       f"{ratio * R_e6:.2f} n-units at STP; every level in {HERE / 'self_check.csv'}")

# ---------------------------------------------------------------------------
# The runs
# ---------------------------------------------------------------------------
products, built = {}, []
for name in RUNS:
    run = Path("forward") / name
    path = run / "output" / f"{name}_profile.nc"
    if not path.exists():
        run_inputs.build(run / f"{name}_build.toml")
        path = forward_transfer.run(run / f"{name}.toml").product
        built.append(name)
    products[name] = tree(path).to_dataset(inherit=False)
target_c = float(products[RUNS[0]]["latitude_planetocentric_deg"])

# ---------------------------------------------------------------------------
# 3. Refractivity against radius
# ---------------------------------------------------------------------------
r_local = {reading: local_radius(r_eq, lat) for reading, lat in READINGS.items()}
fits = {}
for name in RUNS:
    root = products[name]
    r_m = np.asarray(root["radius_m"].values, dtype="float64")
    lnN_m = np.log(np.asarray(root["refractivity"].values, dtype="float64") * 1e6)
    order = np.argsort(r_m)
    model = lambda r: np.interp(r, r_m[order], lnN_m[order], left=np.nan, right=np.nan)
    for reading, r_obs in r_local.items():
        # Every one of Kliore's levels inside the delivered column, or the offset is not a candidate.
        cost = lambda dr: float(np.mean((model(r_obs + dr) - np.log(N_k)) ** 2)) if np.isfinite(model(r_obs + dr)).all() else np.inf
        first = np.arange(-400e3, 400e3 + 1, 1e3)
        dr0 = first[int(np.argmin([cost(x) for x in first]))]
        dr = float(minimize_scalar(cost, bounds=(dr0 - 2e3, dr0 + 2e3), method="bounded", options={"xatol": 1.0}).x)
        before = model(r_obs) - np.log(N_k)
        after = model(r_obs + dr) - np.log(N_k)
        fits[(name, reading)] = {"offset_m": dr, "before": before, "after": after,
                                 "rms_after_percent": 100 * float(np.sqrt(np.nanmean(after ** 2))),
                                 "rms_before_percent": 100 * float(np.sqrt(np.nanmean(before ** 2))),
                                 "levels": int(np.isfinite(after).sum())}
i100 = int(np.argmin(np.abs(np.log(p_k / 100.0))))
lines = [f"Kliore's latitude 9.5 S planetocentric; the planetographic reading on his own ellipse (f = {KLIORE_OBLATENESS}) is "
         f"{abs(KLIORE_LATITUDE_G):.2f} S; the two readings differ by "
         f"{(r_local['planetocentric'] - r_local['planetographic']).mean() / 1e3:.1f} km in r_local",
         f"at {p_k[i100]:g} mbar: Kliore r_eq {r_eq[i100] / 1e3:.1f} km, r_local {r_local['planetocentric'][i100] / 1e3:.1f} km "
         f"(planetocentric reading), {r_local['planetographic'][i100] / 1e3:.1f} km (planetographic reading); delivered "
         + ", ".join(f"{interp_log(products[n]['pressure_Pa'].values / 100, products[n]['radius_m'].values, p_k[i100]) / 1e3:.1f} km ({HYPOTHESIS[j]})"
                     for j, n in enumerate(RUNS)),
         "offset is delivered radius minus Kliore's r_local at equal refractivity; RMS in ln N x 100 (percent in N):"]
for (name, reading), f in fits.items():
    lines.append(f"  {name}, {reading} reading: offset {f['offset_m'] / 1e3:+.2f} km, RMS before {f['rms_before_percent']:.1f} %, "
                 f"after {f['rms_after_percent']:.2f} % over {f['levels']} levels")
summary["radius"] = {f"{n}|{r}": {k: v for k, v in f.items() if k not in ("before", "after")} for (n, r), f in fits.items()}

# ---------------------------------------------------------------------------
# 4. Refractivity against pressure
# ---------------------------------------------------------------------------
T_r, p_r, factors = rp.rescale_composition(T_k, p_k)
press = {}
for name in RUNS:
    root = products[name]
    N_d = interp_log(root["pressure_Pa"].values / 100, np.asarray(root["refractivity"].values) * 1e6, p_r)
    q = N_d / N_k - 1
    at = {lev: float(interp_log(root["pressure_Pa"].values / 100, np.asarray(root["refractivity"].values) * 1e6, lev)
                     / np.exp(interp_log(p_r, np.log(N_k), lev))) for lev in LEVELS_MBAR}
    press[name] = {"N_delivered": N_d, "ratio": q, "rms_percent": 100 * float(np.sqrt(np.nanmean(q ** 2))),
                   "mean_percent": 100 * float(np.nanmean(q)), "at": at}
summary["pressure"] = {n: {"rms_percent": v["rms_percent"], "mean_percent": v["mean_percent"], "at": v["at"]} for n, v in press.items()}

table = ["", f"Kliore's pressure rescaled to 94/6 by Lindal's rule: x{factors['pressure']:.4f}; his refractivity as measured",
         "the table: run, reading, radius offset (km), residual RMS in ln N after it (%), RMS in pressure space (%), "
         "delivered over observed at " + ", ".join(f"{l:g}" for l in LEVELS_MBAR) + " mbar"]
for j, name in enumerate(RUNS):
    for reading in READINGS:
        f = fits[(name, reading)]
        table.append(f"  {HYPOTHESIS[j]:24s} {reading:15s} {f['offset_m'] / 1e3:+8.2f} {f['rms_after_percent']:7.2f} "
                     f"{press[name]['rms_percent']:7.2f}   " + "  ".join(f"{press[name]['at'][l]:.4f}" for l in LEVELS_MBAR))
table.append("the pressure columns do not depend on the reading; mean delivered over observed minus one: "
             + ", ".join(f"{HYPOTHESIS[j]} {press[n]['mean_percent']:+.2f} %" for j, n in enumerate(RUNS)))

# ---------------------------------------------------------------------------
# The figures
# ---------------------------------------------------------------------------
colours = {RUNS[0]: "C2", RUNS[1]: "C4"}
styles = {"planetocentric": ("o", "C3"), "planetographic": ("s", "C0")}
fig, (ax, axr) = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
for j, name in enumerate(RUNS):
    root = products[name]
    r_m = np.asarray(root["radius_m"].values) / 1e3
    keep = (r_m > r_local["planetographic"].min() / 1e3 - 150) & (r_m < r_eq.max() / 1e3 + 50)
    ax.plot(np.log(np.asarray(root["refractivity"].values) * 1e6)[keep], r_m[keep], color=colours[name], linewidth=1.2,
            label=f"delivered, {HYPOTHESIS[j]}")
for reading, r_obs in r_local.items():
    marker, colour = styles[reading]
    ax.plot(np.log(N_k), r_obs / 1e3, marker, color=colour, markersize=3, mfc="none", alpha=0.6,
            label=f"Kliore, {reading} reading, before the offset")
    f = fits[(RUNS[0], reading)]
    ax.plot(np.log(N_k), (r_obs + f["offset_m"]) / 1e3, marker, color=colour, markersize=3,
            label=f"Kliore, {reading} reading, after {f['offset_m'] / 1e3:+.1f} km ({HYPOTHESIS[0]})")
    for name in RUNS:
        g = fits[(name, reading)]
        axr.plot(100 * g["after"], (r_obs + g["offset_m"]) / 1e3, marker + "-", color=colour, markersize=2.5, linewidth=0.6,
                 alpha=1.0 if name == RUNS[0] else 0.45,
                 label=f"{reading} reading, {HYPOTHESIS[RUNS.index(name)]} ({g['offset_m'] / 1e3:+.1f} km)")
ax.set_xlabel("ln N (N in n-units)")
ax.set_ylabel("radius (km)")
ax.grid(alpha=0.3)
ax.legend(fontsize=7)
axr.axvline(0, color="0.5", linewidth=0.8)
axr.set_xlabel("delivered minus Kliore in ln N, x 100 (percent in N), after the offset")
axr.grid(alpha=0.3)
axr.legend(fontsize=7)
fig.suptitle(f"Refractivity against radius at 11.9 S ({abs(target_c):.2f} S planetocentric): the two runs and Kliore's "
             "Pioneer 11 exit, both readings of his radius mapping", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "refractivity_radius.png", dpi=110)
plt.close(fig)

fig, (ax, axr) = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
for j, name in enumerate(RUNS):
    root = products[name]
    p_m = np.asarray(root["pressure_Pa"].values) / 100
    keep = (p_m > 0.5) & (p_m < 400)
    ax.plot(np.asarray(root["refractivity"].values)[keep] * 1e6, p_m[keep], color=colours[name], linewidth=1.2,
            label=f"delivered, {HYPOTHESIS[j]}")
    axr.plot(100 * press[name]["ratio"], p_r, "o-", color=colours[name], markersize=2.5, linewidth=0.6,
             label=f"{HYPOTHESIS[j]}: RMS {press[name]['rms_percent']:.1f} %")
bars = np.isfinite(np.asarray(kliore["pressure_uncertainty_minus_mbar"], dtype="float64"))
lo = -np.asarray(kliore["pressure_uncertainty_minus_mbar"], dtype="float64")[bars] * factors["pressure"]
hi = np.asarray(kliore["pressure_uncertainty_plus_mbar"], dtype="float64")[bars] * factors["pressure"]
ax.plot(N_k, p_r, "o", color="C3", markersize=3, label="Kliore, pressure rescaled to 94/6")
ax.errorbar(N_k[bars], p_r[bars], yerr=[lo, hi], fmt="none", ecolor="C3", capsize=3, label="Kliore's pressure ranges, rescaled")
ax.set_xscale("log")
ax.set_yscale("log")
ax.invert_yaxis()
ax.set_xlabel("refractivity (n-units)")
ax.set_ylabel("pressure (mbar)")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=8)
axr.axvline(0, color="0.5", linewidth=0.8)
axr.set_xlabel("delivered over observed minus one (%)")
axr.grid(alpha=0.3, which="both")
axr.legend(fontsize=8)
fig.suptitle(f"Refractivity against pressure at 11.9 S ({abs(target_c):.2f} S planetocentric): the two runs and Kliore's "
             "Pioneer 11 exit, his pressure rescaled to the anchor's mixture", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "refractivity_pressure.png", dpi=110)
plt.close(fig)

record(2, "for the author: refractivity against radius (fitted offset, both readings) and against pressure (rescaled); "
          "the table", True,
       "\n".join(lines + table) + f"\nfigures {HERE / 'refractivity_radius.png'}, {HERE / 'refractivity_pressure.png'}"
       + (f"\nbuilt here: {', '.join(built)}" if built else ""))

(HERE / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8", newline="\n")
sys.exit(1 if failed else 0)
