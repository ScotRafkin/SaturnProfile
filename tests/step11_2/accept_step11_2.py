"""Acceptance for SPEC_11 v0.4 Step 2: the Lindal anchor transferred under the Lindal-only wind.

Builds `forward/lindal_iris_ii_110/` with `casspian-run-inputs` (its `[shear]` the `lindal_iris`
case, case (ii) at 110 mbar) and runs it with the transfer driver, to 26.4 S planetocentric. The
run's inputs and product are written in the run directory, as a run's are, and are not committed;
this script writes its figures and numbers under `reports/step11_2/`.

**The relaxation.** The tool and the shear case are uncommitted, so the run's inputs carry `-dirty`,
and the one refusal point, `casspian.lib.control._refuse_dirty_commit`, is replaced for this script
only (SPEC_03 section 0, the in-memory candidate rule). Nothing else is patched.

Two conversions, stated here and in the report:

* **The pressure identity in kelvins.** The residual is `p_produced / p_label - 1`. At fixed
  refractivity `T = p R_bar / (k_B N)` is proportional to `p`, so a residual `d` is a temperature
  error `T d`.
* **The dry adiabat.** `g / c_p`, with `g = dPhi/dz` along the target's local vertical from the
  product, and `c_p` per unit mass of the run's composition with frozen molecular rotation: 7/2 R
  per mole of H2, 5/2 R of He and 4 R of NH3, divided by the mean molar mass.

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
from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.run import run_inputs

HERE = Path("reports/step11_2")
HERE.mkdir(parents=True, exist_ok=True)
RUN = Path("forward/lindal_iris_ii_110")
TABLE = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)
LEVELS = {110: 11000.0, 290: 29000.0, 730: 73000.0}
STATED_TARGET_GRAPHIC_DEG = -31.35
EXPLAIN_K = 0.5
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def tree(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


def at_pressure(p, values, targets):
    order = np.argsort(np.log(p))
    return np.interp(np.log(targets), np.log(p)[order], np.asarray(values)[order])


_original_refuse = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(attrs.get("casspian_git_commit", ""))
print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this script only\n")
warnings.simplefilter("ignore")

written = run_inputs.build(RUN / f"{RUN.name}_build.toml")
state = forward_transfer.run(RUN / f"{RUN.name}.toml")
product = tree(state.product)
root = product.to_dataset(inherit=False)
anchor_root = product["anchors/lindal"].to_dataset(inherit=False)
thermo = product["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)
wind = product["inputs/wind"].to_dataset(inherit=False)
ctl._refuse_dirty_commit = _original_refuse

p = np.asarray(root["pressure_Pa"].values, dtype="float64")
T = np.asarray(root["temperature_K"].values, dtype="float64")
p_anchor = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
T_anchor = np.asarray(thermo["temperature_K"].values, dtype="float64")
target_c = float(root["latitude_planetocentric_deg"])
anchor_c = float(anchor_root["latitude_planetocentric_deg"])
anchor_g = float(anchor_root["latitude_planetographic_deg"])
lat_c = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
lat_g = np.asarray(wind["latitude_planetographic_deg"].values, dtype="float64")
target_g = float(np.interp(target_c, lat_c, lat_g))
identity = np.asarray(root["pressure_identity_residual"].values, dtype="float64")
passes = int(state.state.record["passes"])

# ---------------------------------------------------------------------------
# 1. The round trip
# ---------------------------------------------------------------------------
lines, rows = [], {}
for level, pressure in LEVELS.items():
    delivered = float(at_pressure(p, T, [pressure])[0] - at_pressure(p_anchor, T_anchor, [pressure])[0])
    use = TABLE["pressure_mbar"] == level
    f = iris.fit(TABLE["latitude_planetographic_deg"][use], TABLE["temperature_K"][use],
                 grid=np.array([anchor_g, target_g]))
    iris_change = float(f.value_K[1] - f.value_K[0])
    difference = delivered - iris_change
    rows[level] = {"delivered_K": delivered, "iris_K": iris_change, "difference_K": difference,
                   "iris_at_anchor_K": float(f.value_K[0]), "iris_at_target_K": float(f.value_K[1])}
    lines.append(f"{level} mbar: delivered T(target) - T(anchor) {delivered:+.2f} K; the IRIS fit's change "
                 f"{iris_change:+.2f} K ({f.value_K[0]:.2f} to {f.value_K[1]:.2f} K); difference "
                 f"{difference:+.2f} K{'  (over 0.5 K: explained in the report)' if abs(difference) > EXPLAIN_K else ''}")
record(1, "the round trip: the delivered temperature change from the anchor to the target against the "
          "IRIS fit's change between the same two latitudes, at 110, 290 and 730 mbar (reported)", True,
       f"anchor {anchor_c:.3f} planetocentric, {anchor_g:.2f} planetographic; target {target_c:.2f} "
       f"planetocentric, {target_g:.2f} planetographic by the wind tool's rule (stated "
       f"{STATED_TARGET_GRAPHIC_DEG}); outer loop passes {passes}\n" + "\n".join(lines))

# ---------------------------------------------------------------------------
# 2. For the author: the profile, its lapse rate, the identity in kelvins
# ---------------------------------------------------------------------------
z = np.asarray(root["altitude_m"].values, dtype="float64")
Phi = np.asarray(root["geopotential_m2s2"].values, dtype="float64")
m_bar = np.asarray(root["mean_molar_mass_kg_mol"].values, dtype="float64")
composition = product["inputs/composition"]
comp_root = composition.to_dataset(inherit=False)
column = int(np.argmin(np.abs(np.asarray(comp_root["latitude_planetocentric_deg"].values) - target_c)))
p_comp = np.asarray(comp_root["pressure_Pa"].values, dtype="float64")
x = {name: at_pressure(p_comp, np.asarray(comp_root[f"x_{name}"].values)[:, column], p)
     for name in ("H2", "He", "NH3")}
c_p = MOLAR_GAS_CONSTANT * (3.5 * x["H2"] + 2.5 * x["He"] + 4.0 * x["NH3"]) / m_bar
order = np.argsort(z)
g = np.gradient(Phi[order], z[order])
lapse = -np.gradient(T[order], z[order])
adiabat = g / c_p[order]
identity_K = T * identity

fig, axes = plt.subplots(1, 3, figsize=(15, 6.5), sharey=True)
axes[0].plot(T, p / 100.0, color="C3", label=f"delivered at {target_c:g} deg")
axes[0].plot(T_anchor, p_anchor / 100.0, color="k", linewidth=0.9, label=f"anchor at {anchor_c:.2f} deg")
for level, pressure in LEVELS.items():
    axes[0].plot([rows[level]["iris_at_anchor_K"]], [pressure / 100.0], "ko", markersize=4)
    axes[0].plot([rows[level]["iris_at_target_K"]], [pressure / 100.0], "o", color="C3", markersize=4)
axes[0].set_xlabel("temperature (K); dots: the IRIS fit at the two latitudes")
axes[0].set_ylabel("pressure (mbar)")
axes[0].legend(fontsize=8)
axes[1].plot(lapse * 1e3, p[order] / 100.0, color="C3", label="delivered -dT/dz")
axes[1].plot(adiabat * 1e3, p[order] / 100.0, color="0.4", linestyle="--", label="dry adiabat g/c_p")
axes[1].set_xlabel("lapse rate (K/km)")
axes[1].legend(fontsize=8)
axes[2].plot(identity_K, p / 100.0, ".-", color="C0", markersize=3)
axes[2].axvline(0, color="0.6", linewidth=0.8)
axes[2].set_xlabel("pressure identity, T (p_produced / p_label - 1), K")
for ax in axes:
    ax.set_yscale("log")
    ax.grid(alpha=0.3)
axes[0].invert_yaxis()
fig.suptitle(f"lindal_iris_ii_110: the Lindal anchor transferred to {target_c:g} deg planetocentric "
             "under the Lindal-only wind, case (ii)", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "target_profile.png", dpi=130)
plt.close(fig)
worst = int(np.argmax(np.abs(identity_K)))
super_adiabatic = p[order][(lapse > adiabat) & (p[order] > 1e3)]
record(2, "for the author: the delivered profile with the anchor's, its lapse rate with the dry "
          "adiabat, and the pressure identity in kelvins level by level", True,
       f"{p.size} levels, {p.min():g} to {p.max():g} Pa; the identity at most {abs(identity_K[worst]):.3f} K "
       f"(residual {identity[worst]:.2e}) at {p[worst]:.4g} Pa; its RMS {np.sqrt(np.mean(identity_K ** 2)):.3f} K\n"
       f"lapse rate above the dry adiabat at {super_adiabatic.size} levels deeper than 10 mbar"
       + (f", {super_adiabatic.min():.4g} to {super_adiabatic.max():.4g} Pa" if super_adiabatic.size else "")
       + f"\nfigure {HERE / 'target_profile.png'}")

(HERE / "summary.json").write_text(json.dumps({"round_trip": rows, "passes": passes,
                                               "target_planetographic_deg": target_g,
                                               "identity_K_max": float(abs(identity_K[worst]))},
                                              indent=2), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
