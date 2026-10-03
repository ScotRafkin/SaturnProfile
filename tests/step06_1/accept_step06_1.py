"""Acceptance checks for SPEC_06 v0.2 Step 1: the transfer kernel's shear term read on the isobar.

Every run is made in a copy under `reports/step06_1/<name>/`, its inputs built there from its
build file, from the working tree with the `-dirty` refusal relaxed in this process only, as in
SPEC_05. The SPEC_05 values compared with are the products set aside before this step under
`reports/step06_1/spec05/` (run 7f's for check 3; SPEC_05's measured numbers are quoted for checks
4 and 5). The five checks of v0.2:

1. the transfer run, no vertical shear, against the registered transfer product: delivered T within
   0.1 K at every level, altitude within 10 m, `r0(10 N)` within 1 m;
2. the closure product does not move (it has no transfer): rebuilt in a copy, as `step04_5` check
   11 does, and compared with the registered one;
3. run 7f's case at 5e4 along the 999 mbar isobar: the largest step of `S/g` between consecutive
   curve nodes that do not straddle a wind file latitude node is at least ten times smaller than
   SPEC_05's;
4. run 8's case at 5e4: the delivered temperature change from the no-shear run at the 998.7 mbar
   level within 0.1 K of zero;
5. run 5's case at 5e4 and 2.5e4: the delivered temperature at the 998.7 mbar level agrees between
   the spacings within 0.5 K;
6. (v0.4) the pressure identity does not depend on the mesh: for every run of the SPEC_05 Step 4
   set the largest identity at 5e4 and at 2.5e4 agree within 10 percent; size and level reported,
   SPEC_05's beside them. The set is rerun first by `tests/step05_4/run_experiments.py` (unless
   `--no-run`), into `reports/step05_4/results.json`; SPEC_05's values are
   `reports/step06_1/spec05/results.json`;
7. (v0.4) the identity at the stop pressure explained by a measurement:
   `tests/step06_1/stop_position.py` (unless `--no-run`) places run 5's and run 6's stop pressure,
   on a wind file node, at five positions `s` across one anchor layer. The offset `p - p_label` must
   arise in that layer alone (under 1 Pa in its neighbours), be linear in `s` (largest departure
   from the straight line under 5 percent of its range), with a slope within 5 percent of the
   layer's mass times the jump in `ln T` across it relative to run 2; and the same law, applied to
   the original runs' breaks at their wind file nodes, must give their measured offsets within 5
   percent of the measured range of the offset in `s` (an absolute bound: run 6's break lies near
   the zero crossing, where a relative bound measures nothing; the measured points interpolated
   at the break are reported beside the law).

Levels are compared by the anchor's level index, as in SPEC_05 Step 4. Writes
`reports/step06_1/output.txt` and `reports/step06_1/F8_sawtooth.png`. Run from the repository
root.
"""

import json
import shutil
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import production as fp
from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.tools.run import run_inputs

HERE = Path("reports/step06_1")
SPEC05 = HERE / "spec05"
TRANSFER = Path("forward/lindal_transfer")
CLOSURE = Path("forward/lindal_closure")
REGISTERED_TRANSFER = TRANSFER / "output" / "lindal_transfer_profile.nc"
REGISTERED_CLOSURE = CLOSURE / "output" / "lindal_closure_profile.nc"
LEVEL_PRESSURE_PA = 1.0e5  # the level nearest 1 bar, 998.7 mbar, level 59
results = []
_started = time.time()
NEWLINE = chr(10)


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def tree_of(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(path, engine="netcdf4").load()
    tree.close()
    return tree


def root_of(path):
    return tree_of(path).to_dataset(inherit=False)


def values(ds, name):
    return np.asarray(ds[name].values, dtype="float64")


def copy_run(source: Path, name: str, spacing: str = "5.0e4", figures: bool = False) -> Path:
    """The run's control files under HERE/name, paths for the copy's depth, inputs built there.

    Returns the namelist. `spacing` replaces the namelist's geopotential spacing; `figures`
    keeps or turns off F5 to F9.
    """
    copy = HERE / name
    if copy.exists():
        shutil.rmtree(copy)
    copy.mkdir(parents=True)
    run = source.name
    for suffix in ("_build.toml", ".toml"):
        text = (source / f"{run}{suffix}").read_text(encoding="utf-8")
        text = text.replace('"../../', '"../../../')
        if suffix == ".toml":
            for old, new in (("geopotential_spacing_m2s2 = 5.0e4",
                              f"geopotential_spacing_m2s2 = {spacing}"),
                             ("figures = true", f"figures = {'true' if figures else 'false'}")):
                if text.count(old) != 1:
                    raise RuntimeError(f"{run}{suffix}: {old!r} matched {text.count(old)} times")
                text = text.replace(old, new)
        (copy / f"{run}{suffix}").write_bytes(text.encode("utf-8"))
    run_inputs.build(copy / f"{run}_build.toml")
    return copy / f"{run}.toml"


def transfer_run(source, name, **options):
    started = time.time()
    product = Path(forward_transfer.run(copy_run(Path(source), name, **options)).product)
    print(f"{name}: {time.time() - started:.0f} s, {product}", flush=True)
    return product


original_refusal = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
    attrs.get("casspian_git_commit", ""))
print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this process only\n")
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        products = {
            "transfer": transfer_run(TRANSFER, "lindal_transfer"),
            "r7f": transfer_run("forward/shear_r7f_decay_linp_fine", "r7f_5e4", figures=True),
            "r8": transfer_run("forward/shear_r8_increase25", "r8_5e4"),
            "r5": transfer_run("forward/shear_r5_decay20", "r5_5e4"),
            "r5_half": transfer_run("forward/shear_r5_decay20", "r5_2p5e4", spacing="2.5e4"),
        }
        # The closure, as step04_5 check 11 rebuilds it: its namelist copied, its registered
        # inputs, figures off.
        closure_copy = HERE / "lindal_closure"
        if closure_copy.exists():
            shutil.rmtree(closure_copy)
        (closure_copy / "inputs").mkdir(parents=True)
        for source in ctl.read_run_namelist(CLOSURE / "lindal_closure.toml").inputs.values():
            shutil.copy(source, closure_copy / "inputs" / source.name)
        text = (CLOSURE / "lindal_closure.toml").read_text(encoding="utf-8").replace(
            '"../../occul_data/', '"../../../occul_data/').replace("figures = true", "figures = false")
        (closure_copy / "lindal_closure.toml").write_bytes(text.encode("utf-8"))
        products["closure"] = Path(fp.run(closure_copy / "lindal_closure.toml").product)
finally:
    ctl._refuse_dirty_commit = original_refusal

registered = root_of(REGISTERED_TRANSFER)
level = int(np.argmin(np.abs(values(registered, "pressure_label_Pa") - LEVEL_PRESSURE_PA)))
label_hPa = values(registered, "pressure_label_Pa")[level] / 100.0

# ---------------------------------------------------------------------------------------------
# 1. The transfer run, no vertical shear
# ---------------------------------------------------------------------------------------------
new = root_of(products["transfer"])
dT = values(new, "temperature_K") - values(registered, "temperature_K")
dz = values(new, "altitude_m") - values(registered, "altitude_m")
dr0 = float(new["reference_surface_radius_m"].values - registered["reference_surface_radius_m"].values)
dN = np.max(np.abs(values(new, "refractivity") / values(registered, "refractivity") - 1.0))
record(1, "the transfer run (no vertical shear) against the registered transfer product: T within "
          "0.1 K at every level, altitude within 10 m, r0(10 N) within 1 m",
       np.max(np.abs(dT)) <= 0.1 and np.max(np.abs(dz)) <= 10.0 and abs(dr0) <= 1.0,
       f"largest |dT| {np.max(np.abs(dT)):.3e} K at level {int(np.argmax(np.abs(dT)))}; largest "
       f"|d altitude| {np.max(np.abs(dz)):.3e} m; d r0 {dr0:+.3e} m; largest relative dN {dN:.3e}\n"
       f"pressure identity {float(tree_of(products['transfer'])['transfer_record'].attrs['pressure_identity_largest_abs']):.3e} "
       f"(registered "
       f"{float(tree_of(REGISTERED_TRANSFER)['transfer_record'].attrs['pressure_identity_largest_abs']):.3e})")

# ---------------------------------------------------------------------------------------------
# 2. The closure product does not move
# ---------------------------------------------------------------------------------------------
rebuilt, closure = root_of(products["closure"]), root_of(REGISTERED_CLOSURE)
COMPARED = ("pressure_Pa", "temperature_K", "number_density_m3", "mean_refractivity_m3",
            "mean_molar_mass_kg_mol", "refractivity", "radius_m", "height_above_anchor_isobar_m",
            "geopotential_m2s2")
equal = {name: bool(np.array_equal(rebuilt[name].values, closure[name].values)) for name in COMPARED}
record(2, "the closure product, rebuilt in a copy as step04_5 check 11 does, is array-equal to the "
          "registered one in its nine variables (it has no transfer)",
       all(equal.values()), "; ".join(f"{name} {value}" for name, value in equal.items()))

# ---------------------------------------------------------------------------------------------
# 3. The sawtooth along the 999 mbar isobar, run 7f's case at 5e4
# ---------------------------------------------------------------------------------------------
def along_isobar(path):
    """`(latitude_deg, S/g)` along the 999 mbar isobar of the lindal anchor, and the wind's
    latitude nodes. The lindal levels are matched to the union levels at the gauge column."""
    tree = tree_of(path)
    iso = tree["isobars"].to_dataset(inherit=False)
    est = tree["estimate"].to_dataset(inherit=False)
    root = tree.to_dataset(inherit=False)
    latitude = values(iso, "latitude_planetocentric_deg")
    gauge = int(np.argmin(np.abs(latitude - float(root["gauge_latitude_planetocentric_deg"].values))))
    Phi_lindal = values(iso, "geopotential_lindal_m2s2")[:, gauge]
    union = int(np.argmin(np.abs(values(est, "label_pressure_Pa") - LEVEL_PRESSURE_PA)))
    k = int(np.argmin(np.abs(Phi_lindal - values(est, "geopotential_m2s2")[union])))
    wind = tree["inputs/wind"].to_dataset(inherit=False)
    return (latitude, values(iso, "shear_kernel_lindal_per_rad")[k],
            values(wind, "latitude_planetocentric_deg"), values(est, "label_pressure_Pa")[union])


def largest_off_node_step(latitude, s, wind_nodes):
    lo, hi = latitude[:-1], latitude[1:]
    straddles = np.array([np.any((wind_nodes > a - 1e-9) & (wind_nodes <= b + 1e-9))
                          for a, b in zip(lo, hi)])
    steps = np.abs(np.diff(s))
    off = steps[~straddles]
    return float(off.max()), int(straddles.sum()), int((~straddles).sum()), steps


lat_old, s_old, wind_nodes, label_old = along_isobar(SPEC05 / "shear_r7f_decay_linp_fine_profile.nc")
lat_new, s_new, _, label_new = along_isobar(products["r7f"])
old_step, n_on, n_off, steps_old = largest_off_node_step(lat_old, s_old, wind_nodes)
new_step, _, _, steps_new = largest_off_node_step(lat_new, s_new, wind_nodes)
largest_s = float(np.max(np.abs(s_new)))
record(3, "run 7f's case at 5e4, along the 999 mbar isobar: the largest step of S/g between "
          "consecutive curve nodes not straddling a wind file latitude node is at least ten times "
          "smaller than SPEC_05's",
       new_step * 10.0 <= old_step,
       f"isobar label {label_new:.1f} Pa (SPEC_05 {label_old:.1f} Pa); {n_off} off-node steps, "
       f"{n_on} straddling a wind node\n"
       f"largest off-node step: SPEC_05 {old_step:.3e} per rad, now {new_step:.3e} per rad, ratio "
       f"{old_step / new_step if new_step else float('inf'):.1f}\n"
       f"largest |S/g| on the isobar: SPEC_05 {np.max(np.abs(s_old)):.3e}, now {largest_s:.3e}; "
       f"largest step at a wind node now {float(np.max(steps_new)):.3e}")

from matplotlib.figure import Figure  # noqa: E402

fig = Figure(figsize=(8.0, 5.0))
ax = fig.subplots()
ax.plot(lat_old, s_old, linewidth=0.9, label="SPEC_05 (S/g bilinear on the mesh)")
ax.plot(lat_new, s_new, linewidth=0.9, label="SPEC_06 (the wind read on the isobar)")
ax.set_xlabel("planetocentric latitude (deg)")
ax.set_ylabel("S/g along the isobar (per rad)")
ax.set_title(f"Run 7f at 5e4, the {label_new / 100:.1f} mbar isobar", fontsize=9)
ax.legend(fontsize=8)
fig.savefig(HERE / "F8_sawtooth.png", dpi=130)

# ---------------------------------------------------------------------------------------------
# 4. The kink does not leak across p_s, run 8's case at 5e4
# ---------------------------------------------------------------------------------------------
r8 = root_of(products["r8"])
dT8 = values(r8, "temperature_K") - values(new, "temperature_K")
record(4, f"run 8's case at 5e4: the delivered temperature change from the no-shear run at the "
          f"{label_hPa:.1f} mbar level (above p_s) within 0.1 K of zero",
       abs(dT8[level]) <= 0.1,
       f"dT at level {level}: {dT8[level]:+.4f} K (SPEC_05 -5.872 K); the levels below it: "
       + ", ".join(f"{v:+.3f}" for v in dT8[level + 1:])
       + f"\nthe no-shear run is check 1's transfer run, which is run 2 (uniform, c = 1)")

# ---------------------------------------------------------------------------------------------
# 5. The kink level does not depend on the mesh, run 5's case
# ---------------------------------------------------------------------------------------------
r5, r5_half = root_of(products["r5"]), root_of(products["r5_half"])
T5, T5h = values(r5, "temperature_K"), values(r5_half, "temperature_K")
dT5, dT5h = T5 - values(new, "temperature_K"), T5h - values(new, "temperature_K")
zone = int(np.argmax(np.abs(dT5)))
record(5, f"run 5's case: the delivered temperature at the {label_hPa:.1f} mbar level agrees between "
          f"5e4 and 2.5e4 within 0.5 K",
       abs(T5[level] - T5h[level]) <= 0.5,
       f"T at level {level}: {T5[level]:.4f} K at 5e4, {T5h[level]:.4f} K at 2.5e4, difference "
       f"{T5[level] - T5h[level]:+.4f} K\n"
       f"dT from the no-shear run at level {level}: {dT5[level]:+.3f} K at 5e4, {dT5h[level]:+.3f} K "
       f"at 2.5e4 (SPEC_05 -8.746 and -9.058); in the zone, largest {dT5[zone]:+.3f} K at level {zone}")

# ---------------------------------------------------------------------------------------------
# 6. The identity does not depend on the mesh: the SPEC_05 Step 4 set
# ---------------------------------------------------------------------------------------------
if "--no-run" not in sys.argv:
    sys.path.insert(0, str(Path("tests/step05_4").resolve()))
    import run_experiments  # noqa: E402
    run_experiments.main()
now = json.loads(Path("reports/step05_4/results.json").read_text(encoding="utf-8"))
before = json.loads((SPEC05 / "results.json").read_text(encoding="utf-8"))
lines, ok = [], True
for run, entry in now.items():
    cells, largest = [], {}
    for spacing in ("5e4", "2.5e4"):
        e, b = entry.get(spacing), before.get(run, {}).get(spacing, {})
        if e is None:
            continue
        if e["status"] != "completed":
            ok = False
            cells.append(f"{spacing} {e['status']}, {e.get('failure_type')}")
            continue
        largest[spacing] = e["pressure_identity_largest_abs"]
        old = b.get("pressure_identity_largest_abs") if b.get("status") == "completed" else None
        cells.append(f"{spacing} {largest[spacing]:.3e} at level {e['pressure_identity_largest_level']}"
                     f" (SPEC_05 {'n/a' if old is None else f'{old:.3e}'}), "
                     f"{e['outer_loop_passes']} passes")
    if len(largest) == 2:
        apart = abs(largest["5e4"] / largest["2.5e4"] - 1.0)
        ok = ok and apart <= 0.10
        cells.append(f"apart {100 * apart:.1f} percent")
    lines.append(f"{run}: " + "; ".join(cells))
record(6, "the largest pressure identity of every run of the SPEC_05 Step 4 set agrees between 5e4 and "
          "2.5e4 within 10 percent (7f at 5e4 only); size and level reported, not bounded",
       ok, NEWLINE.join(lines))

# ---------------------------------------------------------------------------------------------
# 7. The identity at the stop pressure, explained by a measurement
# ---------------------------------------------------------------------------------------------
if "--no-run" not in sys.argv:
    import runpy  # noqa: E402
    runpy.run_path("tests/step06_1/stop_position.py", run_name="__main__")
measured = json.loads((HERE / "stop_position.json").read_text(encoding="utf-8"))
T_run2 = values(root_of(now["shear_r2_uniform"]["5e4"]["product"]), "temperature_K")
labels_run2 = values(root_of(now["shear_r2_uniform"]["5e4"]["product"]), "pressure_label_Pa")


def ln_T_jump(T, top):
    """The jump in `ln T` from level `top` to `top + 1`, relative to run 2."""
    return float(np.log(T[top + 1] / T_run2[top + 1]) - np.log(T[top] / T_run2[top]))


lines, ok = [], True
fits = {}
for run, top in (("shear_r5_decay20", 15), ("shear_r6_decay40", 27)):
    rows = {float(k.split()[1]): v for k, v in measured.items() if k.split()[0] == run}
    s_values = np.array(sorted(rows))
    inside = np.array([rows[x]["offset_Pa"][str(top + 1)] - rows[x]["offset_Pa"][str(top)]
                       for x in s_values])
    outside = max(abs(rows[x]["offset_Pa"][str(k + 1)] - rows[x]["offset_Pa"][str(k)])
                  for x in s_values for k in range(top - 2, top + 3) if k != top)
    slope, intercept = np.polyfit(s_values, inside, 1)
    departure = float(np.max(np.abs(inside - (slope * s_values + intercept))) / np.ptp(inside))
    T_mid = values(root_of(f"reports/step06_1/stop_position/{run}_s0.50/output/{run}_profile.nc"),
                   "temperature_K")
    mass = float(labels_run2[top + 1] - labels_run2[top])
    expected = mass * abs(ln_T_jump(T_mid, top))
    zero = -intercept / slope
    fits[run] = (zero, slope / expected, float(np.ptp(inside)), s_values, inside)
    passed = outside < 1.0 and departure < 0.05 and abs(slope / expected - 1.0) <= 0.05
    ok = ok and passed
    lines.append(f"{run}, layer {top} to {top + 1}: offset gained there at s = "
                 + ", ".join(f"{x:.2f}: {v:+.2f}" for x, v in zip(s_values, inside))
                 + f" Pa; in the neighbouring layers at most {outside:.2f} Pa\n"
                 f"    linear in s: slope {slope:.1f} Pa, zero at s = {zero:.3f}, largest departure "
                 f"{100 * departure:.1f} percent of the range; layer mass {mass:.1f} Pa times "
                 f"|jump in ln T| {abs(ln_T_jump(T_mid, top)):.4f} = {expected:.1f} Pa, ratio "
                 f"{slope / expected:.4f}")


def law(run, top, node_Pa):
    """The offset the law gives for a break at a wind file node inside layer `top`: the layer's
    mass times its jump in ln T, times (s - s0), with s0 and the ratio from the measurement."""
    product = now[run]["5e4"]["product"]
    d = root_of(product)
    T = values(d, "temperature_K")
    s = float(np.log(node_Pa / labels_run2[top]) / np.log(labels_run2[top + 1] / labels_run2[top]))
    zero, ratio = fits[run][:2]
    mass = float(labels_run2[top + 1] - labels_run2[top])
    return ratio * mass * abs(ln_T_jump(T, top)) * (s - zero), s


for run, breaks, level in (("shear_r5_decay20", ((14, 630.957), (15, 794.328)), 16),
                           ("shear_r6_decay40", ((27, 7943.28),), 28)):
    d = root_of(now[run]["5e4"]["product"])
    offset = float(values(d, "pressure_Pa")[level] - values(d, "pressure_label_Pa")[level])
    parts = [law(run, top, node) for top, node in breaks]
    predicted = sum(v for v, _ in parts)
    span = fits[run][2]
    passed = abs(predicted - offset) <= 0.05 * span
    ok = ok and passed
    s_points, measured_points = fits[run][3], fits[run][4]
    measured_layer = {"shear_r5_decay20": 15, "shear_r6_decay40": 27}[run]
    interpolated = [f"{float(np.interp(s, s_points, measured_points)):+.2f} Pa (layer {top})"
                    for (top, _), (_, s) in zip(breaks, parts) if top == measured_layer]
    lines.append(f"{run} as run, offset at level {level}: measured {offset:+.2f} Pa, the law "
                 f"{predicted:+.2f} Pa from its breaks at "
                 + ", ".join(f"{node:g} Pa (layer {top}, s = {s:.3f}, {v:+.2f} Pa)"
                             for (top, node), (v, s) in zip(breaks, parts))
                 + f"; apart {abs(predicted - offset):.2f} Pa, {100 * abs(predicted - offset) / span:.1f}"
                 f" percent of the measured range {span:.1f} Pa; the measured points interpolated at "
                 f"the break's s in the measured layer: " + ", ".join(interpolated))
record(7, "the identity at the stop pressure: the offset arises in the layer holding the wind's break, "
          "linear in the break's position s, slope the layer mass times the jump in ln T, and the "
          "same law gives the original runs' offsets",
       ok, NEWLINE.join(lines))

failed = [r for r in results if not r[2]]
summary = f"{len(results) - len(failed)} of {len(results)} checks pass, {time.time() - _started:.0f} s"
print(summary)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{summary}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
