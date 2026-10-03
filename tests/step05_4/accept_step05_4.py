"""Acceptance checks for SPEC_05 v0.7 Step 4: the named experiments.

Runs `tests/step05_4/run_experiments.py` first (all eleven runs, both spacings) unless `--no-run`
is given, then reads `reports/step05_4/results.json` and the products. The four checks of v0.7:
run 2 reproduces the registered transfer product; runs 3a and 3b, the exact null, agree in N, p
and T to round-off; run 3a's reference surface at 10 N is the reviewing agent's no-wind value; the
largest pressure identity of every run agrees between the two spacings within 10 percent, its size
reported (REVIEW_05_step4 bounded it at 1e-2; restated by SPEC_06 v0.4).
Every other outcome is reported, not checked (v0.7): the table of the report and one comparison
figure, `reports/figures/step05_4_temperature_minus_run2.png`.

Delivered quantities are compared level by level, by the anchor's level index: every run delivers
the anchor's 66 levels, and a run's isobar labels are its own anchor production under its own wind,
so the same level has slightly different pressures in different runs. Run from the repository root.
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path("reports/step05_4")
FIGURE = Path("reports/figures/step05_4_temperature_minus_run2.png")
REGISTERED = Path("forward/lindal_transfer/output/lindal_transfer_profile.nc")
NULL_TARGET_RELATIVE = 1e-12
R0_NO_WIND_M, R0_BOUND_M = 60_092_307.69, 0.1
# SPEC_06 v0.4: the identity at a break in the wind's shear is the production's integration over
# the anchor's own levels, which no mesh setting reaches; REVIEW_05_step4's bound of 1e-2 is
# restated as mesh independence, the size reported.
MESH_AGREEMENT = 0.10
results = []

if "--no-run" not in sys.argv:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_experiments  # noqa: E402
    run_experiments.main()

RUNS = json.loads((HERE / "results.json").read_text(encoding="utf-8"))


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def root(run, spacing="5e4"):
    entry = RUNS.get(run, {}).get(spacing, {})
    if entry.get("status") != "completed":
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(entry["product"], engine="netcdf4").load()
    tree.close()
    return tree.to_dataset(inherit=False)


def values(ds, name):
    return np.asarray(ds[name].values, dtype="float64")


def largest_relative(a, b):
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    return float(np.max(np.abs(a - b) / np.abs(b)))


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    _tree = xr.open_datatree(REGISTERED, engine="netcdf4").load()
_tree.close()
registered = _tree.to_dataset(inherit=False)
run2 = root("shear_r2_uniform")

# ---------------------------------------------------------------------------
# 1. Run 2 reproduces the registered transfer product
# ---------------------------------------------------------------------------
NAMES = ("refractivity", "pressure_Pa", "temperature_K", "altitude_m", "reference_surface_radius_m")
same = ({name: bool(np.array_equal(values(run2, name), values(registered, name))) for name in NAMES}
        if run2 is not None else {})
record(1, "run 2 (uniform, c = 1) is array-equal to the registered transfer product in N, p, T, "
          "altitude and r0",
       run2 is not None and all(same.values()),
       "; ".join(f"{name} {value}" for name, value in same.items()) or "run 2 did not complete")

# ---------------------------------------------------------------------------
# 2. Runs 3a and 3b, the exact null, agree in N, p, T
# ---------------------------------------------------------------------------
a, b = root("shear_r3a_nowind"), root("shear_r3b_nowind_anchor")
if a is not None and b is not None:
    null = {name: largest_relative(values(a, name), values(b, name))
            for name in ("refractivity", "pressure_Pa", "temperature_K")}
    record(2, f"runs 3a (10 N) and 3b (the anchor's phi_c) agree in N, p and T to round-off, targeted at "
              f"{NULL_TARGET_RELATIVE:g} relative",
           all(v <= NULL_TARGET_RELATIVE for v in null.values()),
           "largest relative difference, level by level: "
           + "; ".join(f"{name} {value:.3e}" for name, value in null.items())
           + f"\n3a at {float(a['latitude_planetocentric_deg'].values):g} deg, 3b at "
             f"{float(b['latitude_planetocentric_deg'].values):.14g} deg")
else:
    record(2, "runs 3a and 3b agree in N, p and T to round-off", False, "a run did not complete")

# ---------------------------------------------------------------------------
# 3. Run 3a's reference surface at 10 N
# ---------------------------------------------------------------------------
r0 = float(a["reference_surface_radius_m"].values) if a is not None else float("nan")
record(3, f"run 3a's r0(10 N), reference_surface_radius_m, is {R0_NO_WIND_M:,.2f} m within "
          f"{R0_BOUND_M} m (the no-wind surface)",
       abs(r0 - R0_NO_WIND_M) <= R0_BOUND_M,
       f"r0(10 N) {r0:,.3f} m, {r0 - R0_NO_WIND_M:+.3f} m from {R0_NO_WIND_M:,.2f}; under the closure "
       f"wind (run 2) {float(run2['reference_surface_radius_m'].values):,.3f} m")

# ---------------------------------------------------------------------------
# 4. The pressure identity of every completed run at both spacings
# ---------------------------------------------------------------------------
lines, ok = [], True
for run, entry in RUNS.items():
    cells, largest = [], {}
    for spacing in ("5e4", "2.5e4"):
        if spacing not in entry:
            continue
        e = entry[spacing]
        if e["status"] == "completed":
            largest[spacing] = e["pressure_identity_largest_abs"]
            cells.append(f"{spacing}: largest {e['pressure_identity_largest_abs']:.3e} at level "
                         f"{e['pressure_identity_largest_level']}, rms {e['pressure_identity_rms']:.3e}, "
                         f"{e['outer_loop_passes']} passes, {e['wall_s']} s")
        else:
            ok = False
            cells.append(f"{spacing}: {e['status']}, {e['failure_type']}")
    expected = ("5e4",) if run == "shear_r7f_decay_linp_fine" else ("5e4", "2.5e4")
    present = all(s in entry for s in expected)
    ok = ok and present
    if len(largest) == 2:
        apart = abs(largest["5e4"] / largest["2.5e4"] - 1.0)
        ok = ok and apart <= MESH_AGREEMENT
        cells.append(f"apart {100 * apart:.1f} percent")
    lines.append(f"{run}: " + "; ".join(cells))
record(4, f"every run completes, and the largest pressure identity agrees between 5e4 and 2.5e4 within "
          f"{100 * MESH_AGREEMENT:g} percent (7f at 5e4 only); its size is reported, not bounded "
          f"(SPEC_06 v0.4)",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# The results table and the comparison figure (reported, not checked)
# ---------------------------------------------------------------------------
def level_of(ds, pressure_Pa):
    return int(np.argmin(np.abs(values(ds, "pressure_label_Pa") - pressure_Pa)))


rows = ["| Run | Spacing | Largest dT from run 2 (K) | at p (Pa) | dT at 1 bar | dT at gauge | dT at top "
        "| dz at 1 bar (m) | dz at gauge (m) | dz at top (m) | r0 (m) | identity | passes | wall (s) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
T2, z2 = values(run2, "temperature_K"), values(run2, "altitude_m")
bar, gauge = level_of(run2, 1.0e5), int(run2["gauge_level_index"].values)
for run, entry in RUNS.items():
    for spacing in ("5e4", "2.5e4"):
        e = entry.get(spacing)
        if e is None:
            continue
        if e["status"] != "completed":
            rows.append(f"| {run} | {spacing} | stopped: {e['failure_type']} | | | | | | | | | | | {e['wall_s']} |")
            continue
        ds = root(run, spacing)
        dT, dz = values(ds, "temperature_K") - T2, values(ds, "altitude_m") - z2
        k = int(np.argmax(np.abs(dT)))
        rows.append(f"| {run} | {spacing} | {dT[k]:+.3f} | {values(ds, 'pressure_Pa')[k]:.4g} | "
                    f"{dT[bar]:+.3f} | {dT[gauge]:+.3f} | {dT[0]:+.3f} | {dz[bar]:+.1f} | {dz[gauge]:+.1f} | "
                    f"{dz[0]:+.1f} | {float(ds['reference_surface_radius_m'].values):,.2f} | "
                    f"{e['pressure_identity_largest_abs']:.2e} | {e['outer_loop_passes']} | {e['wall_s']} |")
(HERE / "table.md").write_text("\n".join(rows) + "\n", encoding="utf-8")
print("\n".join(rows))

from matplotlib.figure import Figure  # noqa: E402

fig = Figure(figsize=(7.5, 6.0))
ax = fig.subplots()
for run in ("shear_r3a_nowind", "shear_r3c_half", "shear_r4_decay12", "shear_r5_decay20",
            "shear_r6_decay40", "shear_r7_decay_linp", "shear_r7f_decay_linp_fine",
            "shear_r8_increase25", "shear_r9_increase50"):
    ds = root(run)
    if ds is None:
        continue
    ax.plot(values(ds, "temperature_K") - T2, values(ds, "pressure_Pa") / 100.0, linewidth=1.1,
            label=run.replace("shear_", ""))
ax.axvline(0.0, color="0.5", linewidth=0.8)
ax.set_yscale("log")
ax.invert_yaxis()
ax.set_xlabel("delivered temperature at 10 N minus run 2's, level by level (K)")
ax.set_ylabel("pressure (mbar)")
ax.set_title("SPEC_05 Step 4: the shear experiments against run 2 (run 3b, at the anchor, is not on "
             "this axis)", fontsize=8)
ax.legend(fontsize=7)
FIGURE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(FIGURE, dpi=130)
print(f"figure {FIGURE}")

failed = [r for r in results if not r[2]]
summary = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(summary)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{summary}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
