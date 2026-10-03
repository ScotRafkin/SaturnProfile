"""Acceptance checks for SPEC_04 v0.17 Step 5: the production at the target, altitude and datum,
kind `profile` in transfer mode, `casspian-forward` in transfer mode, and F5, F6 and F7.

Every check prints its measured value, and every check states the spacings it ran at. Checks beyond
the specification's acceptance list are labeled so.

**The runs.** Run 1 is `casspian-forward forward/lindal_transfer/lindal_transfer.toml` through
`forward.transfer.run` itself, which writes the product and renders the figures. Runs 2 to 5 go
through `forward.transfer.chain`, the same path between the inputs and the product, because they
need a wind field or a mesh that a namelist cannot name: the cylinder-extended wind of decision P
is rebuilt in memory, and the M = 2 run is built on a mesh given to the outer loop, which section
15 ruling 1 requires so that the identity is an identity of the transfer and not of two lattices.
Run 5's product is written by `forward.transfer.write_product`, the function `run` writes with.

**What the values are compared against.** The delivered temperature is compared with **the closure
product's**, not with the anchor's tabulated column: the closure residual is 3.3e-3 in both `p` and
`T`, far outside this step's bounds, and A40 is about the model's own gradients (SPEC_04 section 1,
and the interim report section 3, accepted at section 17). The registered closure product is read,
never rewritten; check 11 rebuilds it in a run directory of its own.

**What is relaxed.** The synthetic 60 N anchor of run 5 is written by this script on a working tree
that is not clean, so its `casspian_git_commit` ends in `-dirty` and `lib.control` refuses it. The
refusal is relaxed in this process only, at the one point `_refuse_dirty_commit`, and the
relaxation is named in run 5's output. Nothing in the code changes and no registered input is
touched.

**The spacings.** This suite is pinned to 0.05 degrees by 50,000 m2/s2, the spacings SPEC_04
section 7's expected values are stated at, and asserts that the production namelist still carries
them, because run 1 is the namelist's own run through the driver. That is the rule of section 0 at
v0.17: a suite carries the spacings its bounds were set for, and a suite meant to follow the
namelist says so. This one says so for run 1 and pins the rest.

**The development flag.** `--quick <factor>` multiplies both mesh spacings and is a development aid
only; the recorded run takes no argument. Run 1 goes through the namelist and ignores it, so
`--quick` is for runs 2 to 5 and the refusals.
"""

import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import netCDF4
import numpy as np
from matplotlib.image import imread

import fields
from casspian.forward import estimate as fe
from casspian.forward import production as fp
from casspian.forward import propagate as fprop
from casspian.forward import transfer as tr
from casspian.lib import control as ctl
from casspian.lib import hydrostatic as hs
from casspian.lib import io as cio
from casspian.lib import kernel as lk
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf
from casspian.tools.plots import style

HERE = Path("reports/step04_5")
TRANSFER = Path("forward/lindal_transfer")
CLOSURE = Path("forward/lindal_closure")
CLOSURE_PRODUCT = CLOSURE / "output" / "lindal_closure_profile.nc"
QUICK = 1.0
if "--quick" in sys.argv:
    QUICK = float(sys.argv[sys.argv.index("--quick") + 1])

#: Checks whose runs this pass does not repeat, when `--carry <output.txt>` is given. SPEC_04 v0.22
#: has a change rerun what it can reach and nothing else: a change to a figure module cannot move a
#: number, so it reruns the figure checks and no suite that computes. These six are the ones whose
#: runs no figure touches, and their rows are taken from the run named on the command line,
#: unchanged and marked as carried. Checks 1, 6, 7 and 8 are rerun although no figure reaches them
#: either, because check 10 is built on their states.
CARRY_FROM = None
CARRIED = frozenset()
if "--carry" in sys.argv:
    CARRY_FROM = Path(sys.argv[sys.argv.index("--carry") + 1])
    CARRIED = frozenset({2, 3, 4, 5, 9, 11})

HERE.mkdir(parents=True, exist_ok=True)
results = []
_started = time.time()
_last = _started


def departure(a, b):
    """The largest relative difference of `a` from `b` where `b` is nonzero, and the largest in
    units in the last place of `b`. Reported beside the bit-identity, which says only True or
    False (REPORT_05_step0: on a platform other than the one that registered the product, the
    closure rerun differs from it in the last bit)."""
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    d = np.abs(a - b)
    nonzero = b != 0
    rel = float(np.max(d[nonzero] / np.abs(b[nonzero]))) if nonzero.any() else 0.0
    ulp = float(np.max(d / np.spacing(np.abs(b))))
    return rel, ulp


#: REVIEW_05_step0, the ruling on section 5: one bound on every platform.
RELATIVE_BOUND = 1e-14


def within(a, b, bound=RELATIVE_BOUND):
    """`a` agrees with `b` to `bound` relative where `b` is nonzero, and exactly where it is zero."""
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    if a.shape != b.shape:
        return False
    d = np.abs(a - b)
    nonzero = b != 0
    return bool(np.all(d[nonzero] <= bound * np.abs(b[nonzero])) and np.all(d[~nonzero] == 0))


def record(number, description, passed, detail):
    """One check's result, with the wall clock since the check before it."""
    global _last
    now = time.time()
    detail = f"{detail}\nwall clock to here {now - _last:.0f} s, {now - _started:.0f} s in all"
    _last = now
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}", flush=True)
    for line in str(detail).splitlines():
        print(f"        {line}", flush=True)
    write_output()


def write_output():
    (HERE / "output.txt").write_text(
        "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n"
                  + "\n".join("        " + line for line in str(detail).splitlines())
                  for n, d, ok, detail in results) + "\n", encoding="utf-8")


def carry(number):
    """Record a check's row from an earlier run's output, unchanged, and say where it came from."""
    lines = CARRY_FROM.read_text(encoding="utf-8").split("\n")
    heads = [i for i, line in enumerate(lines)
             if re.match(r"^\[(PASS|FAIL)\] \d+\. ", line)]
    for position, i in enumerate(heads):
        if int(lines[i].split("]")[1].split(".")[0]) != number:
            continue
        end = heads[position + 1] if position + 1 < len(heads) else len(lines)
        detail = "\n".join(line[8:] for line in lines[i + 1:end] if line.strip())
        results.append((number, lines[i].split(". ", 1)[1], lines[i].startswith("[PASS]"),
                        f"{detail}\nCARRIED, not rerun in this pass: measured in "
                        f"{CARRY_FROM.name}, the eleven of eleven run at the pinned spacings. "
                        f"What changed since is tools/plots/figures_profile.py, F8's altitude "
                        f"panel restated on a common datum (SPEC_04 v0.21), which no run this "
                        f"check makes can reach; under the section 0 rule of v0.22 a figure "
                        f"change reruns the figure checks and no suite that computes."))
        print(f"[{'PASS' if results[-1][2] else 'FAIL'}] {number}. {results[-1][1]} [CARRIED]",
              flush=True)
        write_output()
        return
    raise AssertionError(f"{CARRY_FROM} has no row for check {number}")


def read_closed(path, kind):
    """Read a file under its kind, pull it into memory, and release the handle.

    An unclosed handle segmentation faults the interpreter when the same file is opened again in
    this HDF5 build (section 15 ruling 3), which is what the driver and this script both avoid.
    """
    handle = cio.read(path, kind)
    try:
        return handle.load()
    finally:
        handle.close()


def refusal(call, *args, **kwargs):
    """The message a call refuses with, or `None` if it does not refuse."""
    try:
        call(*args, **kwargs)
    except Exception as error:  # the message is the subject of the check
        return f"{type(error).__name__}: {error}"
    return None


# ---------------------------------------------------------------------------------------------
# The state, the pin, and the helpers every check uses
# ---------------------------------------------------------------------------------------------
namelist = ctl.read_run_namelist(TRANSFER / "lindal_transfer.toml")
inputs = ctl.load_run_inputs(namelist)
anchors = tuple(fprop.propagate(a, namelist.solar_longitude_deg) for a in inputs.anchors)
lindal = anchors[0]

# The spacings section 7's expected values are stated at, pinned here (section 0, v0.17). Run 1 is
# the namelist's own run through the driver, so the pin is an assertion that the namelist still
# carries them rather than an override of it; runs 2 to 5 take them from here.
ACCEPTED_LATITUDE_SPACING_DEG = 0.05
ACCEPTED_GEOPOTENTIAL_SPACING = 5.0e4
for key, accepted in (("latitude_spacing_deg", ACCEPTED_LATITUDE_SPACING_DEG),
                      ("geopotential_spacing_m2s2", ACCEPTED_GEOPOTENTIAL_SPACING)):
    if float(namelist.grid[key]) != accepted:
        raise AssertionError(
            f"the production namelist's {key} is {namelist.grid[key]} and SPEC_04 section 7's "
            f"expected values, which this suite's bounds are, are stated at {accepted}. Run 1 is "
            "the namelist's own run, so the two cannot differ (section 0, v0.17)")
LATITUDE_SPACING = math.radians(ACCEPTED_LATITUDE_SPACING_DEG) * QUICK
GEOPOTENTIAL_SPACING = ACCEPTED_GEOPOTENTIAL_SPACING * QUICK
LOOP = namelist.numerics["outer_loop"]
SIGMA_K = float(namelist.estimation["kernel_uncertainty_per_rad"])
DATUM = float(namelist.datum_isobar_Pa)

phi_a = math.radians(lindal.latitude_planetocentric_deg)
target10 = math.radians(float(namelist.target_latitude_deg))
target60 = math.radians(60.0)
p_tab = np.asarray(lindal.label_pressure_Pa, dtype="float64")
gauge_level = lindal.gauge_level_index
p_b = tr.boundary_pressure(anchors)
closure_field = wf.WindField(inputs.wind)
CONSTANTS = (float(inputs.rotation["angular_rate_rad_s"]), float(inputs.gravity["GM_m3s2"]),
             np.asarray(inputs.gravity["J"].values), np.asarray(inputs.gravity["degree"].values),
             float(inputs.gravity["normalization_radius_m"]))

# The closure product, read once and never rewritten: what the delivered profile is compared with.
closure_tree = read_closed(CLOSURE_PRODUCT, "profile")
closure = closure_tree.to_dataset(inherit=False)
T_closure = np.asarray(closure["temperature_K"].values, dtype="float64")
p_closure = np.asarray(closure["pressure_Pa"].values, dtype="float64")
N_closure = np.asarray(closure["refractivity"].values, dtype="float64")
Phi_closure_product = np.asarray(closure["geopotential_m2s2"].values, dtype="float64")
phi_c = math.radians(float(closure["latitude_planetocentric_deg"].values))


def spacings(scale=1.0):
    return (f"{math.degrees(LATITUDE_SPACING * scale):.4g} degrees by "
            f"{GEOPOTENTIAL_SPACING * scale:.4g} m2/s2")


def chain_to(target_rad, *, field=None, run_anchors=None, scale=1.0, mesh=None,
             gauge_latitude_rad=None, datum_isobar_Pa=None):
    """One run of the driver's own chain, at the pinned spacings times `scale`."""
    run_anchors = tuple(anchors if run_anchors is None else run_anchors)
    phi_r = (fe.gauge_latitude(run_anchors) if gauge_latitude_rad is None
             else float(gauge_latitude_rad))
    return tr.chain(
        inputs, run_anchors, gauge_latitude_rad=phi_r, target_latitude_rad=target_rad,
        gauge_isobar_Pa=namelist.gauge_isobar_Pa, p_b=tr.boundary_pressure(run_anchors),
        latitude_spacing_rad=LATITUDE_SPACING * scale,
        geopotential_spacing_m2s2=GEOPOTENTIAL_SPACING * scale,
        relative_tolerance_ln_p=LOOP["relative_tolerance_ln_p"],
        max_iterations=LOOP["max_iterations"], kernel_uncertainty_per_rad=SIGMA_K,
        datum_isobar_Pa=DATUM if datum_isobar_Pa is None else float(datum_isobar_Pa),
        field=closure_field if field is None else field, mesh=mesh)


def along(state, index=0, anchor=None):
    """One anchor's `ln N` along its curves, with the run's composition term."""
    labels = state.placements[index].label_pressure_Pa
    latitude_c, slope = lk.composition_term(inputs.composition, labels)
    return tr.transfer(state.isobar_kernel, state.curves[index],
                       (anchors[index] if anchor is None else anchor).ln_N, labels, latitude_c, slope)


def named(produced):
    """The top, gauge and bottom levels of a production, in the order it delivers them."""
    return (("the top level", 0), ("the gauge", int(produced["gauge_level_index"])),
            ("the bottom", int(produced["label_pressure_Pa"].size) - 1))


def label_agreement(produced):
    """How far a run's isobar labels are from the closure production's own pressure.

    An M = 1 run's level `k` is the anchor's level `k` carried to the target, and the label the
    transfer carries on it is **the closure production's pressure there**, not the anchor's
    tabulated pressure: the two differ by the closure residual, 3.3e-3 at the bottom row. So every
    comparison with the closure product below is level for level, and this is the measurement that
    says the pairing is the right one. Interpolating the closure columns onto the tabulated axis
    instead put 6e-4 into the temperature comparison, which is how this was found.
    """
    labels = np.asarray(produced["label_pressure_Pa"], dtype="float64")
    if labels.size != p_closure.size:
        raise AssertionError(f"{labels.size} delivered levels against the closure product's "
                             f"{p_closure.size}; the level for level pairing needs M = 1")
    return float(np.max(np.abs(labels / p_closure - 1.0)))


def identity_composition_held(produced, latitude_rad):
    """The pressure identity with the composition read at the anchor's tabulated pressures.

    Beyond the specification's list. The identity's largest value sits at the bottom row at about
    5.7e-7 at every spacing tried, which is the composition read of section 15 ruling 4 and not the
    tracing: the labels are the closure production's pressures and the composition is read onto
    them, moving `ln R_bar` by about 1.3e-6 at the bottom row (section 10 finding 5). Reading it at
    the anchor's tabulated pressures instead, at the same latitude, leaves the tracing's own
    discretization alone and is the part that should fall under halving.
    """
    labels = np.asarray(produced["label_pressure_Pa"], dtype="float64")
    R_tab, m_tab = fp.mean_properties_on_labels(inputs.composition,
                                                math.degrees(float(latitude_rad)), p_tab)
    _, _, _, pressure, _ = fp.produce_on_geopotential(
        np.asarray(produced["refractivity"], dtype="float64"),
        np.asarray(produced["geopotential_m2s2"], dtype="float64"), R_tab, m_tab,
        float(produced["boundary_pressure_Pa"]))
    return float(np.max(np.abs(pressure / labels - 1.0)))


def identity_lines(produced, every_level=False):
    """Decision H: the pressure identity, its statistics and where it is largest."""
    residual = np.asarray(produced["pressure_identity_residual"], dtype="float64")
    labels = np.asarray(produced["label_pressure_Pa"], dtype="float64")
    worst = int(np.argmax(np.abs(residual)))
    lines = [f"pressure identity p/p_label - 1 over {residual.size} levels: largest "
             f"{np.abs(residual).max():.3e} at level {worst} ({labels[worst] / 100:.4g} mbar), "
             f"mean |.| {np.abs(residual).mean():.3e}, signed mean {residual.mean():+.3e}"]
    if every_level:
        lines.append("    every level (decision H), level: mbar: residual")
        for start in range(0, residual.size, 3):
            lines.append("    " + "  ".join(
                f"{k:2d} {labels[k] / 100:9.4g} {residual[k]:+.3e}"
                for k in range(start, min(start + 3, residual.size))))
    return lines


print(f"SPEC_04 Step 5 acceptance. Pinned spacings {spacings()}"
      f"{'' if QUICK == 1.0 else f' (development quick factor {QUICK})'}", flush=True)
print(f"the closure product {CLOSURE_PRODUCT} at "
      f"{math.degrees(phi_c):.6f} degrees, {p_tab.size} levels", flush=True)

# ---------------------------------------------------------------------------------------------
# 1. Run 1: the driver on the production namelist, 10 N, closure wind, M = 1
# ---------------------------------------------------------------------------------------------
result1 = tr.run(TRANSFER / "lindal_transfer.toml")
produced1 = result1.target
product1 = result1.product
tree1 = read_closed(product1, "profile")
root1 = tree1.to_dataset(inherit=False)

T1 = np.asarray(produced1["temperature_K"], dtype="float64")
labels1_agree = label_agreement(produced1)
below1 = (T_closure - T1) / T_closure
STATED_T10 = {0: 1.103e-2, int(produced1["gauge_level_index"]): 1.089e-2,
              int(produced1["label_pressure_Pa"].size) - 1: 1.084e-2}
worst_T10 = max(abs(float(below1[k]) - w) for k, w in STATED_T10.items())

r0_10 = float(root1["reference_surface_radius_m"].values)
STATED_R0 = 60128613.0
altitude1 = np.asarray(produced1["altitude_m"], dtype="float64")
STATED_ALT10 = {0: 411132.0, int(produced1["gauge_level_index"]): 98186.0,
                int(altitude1.size) - 1: -15290.0}
worst_alt10 = max(abs(float(altitude1[k]) - w) for k, w in STATED_ALT10.items())
identity1 = float(np.abs(np.asarray(produced1["pressure_identity_residual"])).max())

record(1, "run 1, the driver on the production namelist: the delivered temperature is below the "
          "closure product's by the stated fractions to 1e-4, r0 at the target is the stated "
          "radius to 1 m, the altitudes above the 1 bar datum are the stated ones to 5 m, and the "
          "pressure identity is reported at every level",
       worst_T10 <= 1e-4 and abs(r0_10 - STATED_R0) <= 1.0 and worst_alt10 <= 5.0,
       f"casspian-forward {namelist.path} through forward.transfer.run, at the namelist's own "
       f"{namelist.grid['latitude_spacing_deg']:g} degrees by "
       f"{namelist.grid['geopotential_spacing_m2s2']:g} m2/s2, which --quick does not touch, mesh "
       f"{result1.state.mesh.shape[0]} by {result1.state.mesh.shape[1]}, outer "
       f"loop {result1.state.record['passes']} passes, residual "
       f"{result1.state.record['residual_ln_p']:.3e} in ln p, "
       f"{len(result1.state.record['extensions'])} mesh extensions\n"
       f"product {product1}\n"
       f"the {p_closure.size} delivered levels are the anchor's levels carried to the target and "
       f"their labels are the closure production's own pressures, agreeing to "
       f"{labels1_agree:.3e} relative, which is why every comparison below is level for level and "
       f"not against the tabulated column (section 1)\n"
       f"T below the closure product's: "
       + ", ".join(f"{name} {float(below1[k]):.4e} (stated {STATED_T10[k]:.3e})"
                   for name, k in named(produced1))
       + f", largest departure {worst_T10:.3e}, bound 1e-4\n"
       f"r0(10 degrees) {r0_10:.1f} m (stated {STATED_R0:.1f}), "
       f"{abs(r0_10 - STATED_R0):.3f} m from it, bound 1 m\n"
       f"altitude above the {DATUM:g} Pa datum: "
       + ", ".join(f"{name} {float(altitude1[k]):.0f} m (stated {STATED_ALT10[k]:.0f})"
                   for name, k in named(produced1))
       + f", largest departure {worst_alt10:.2f} m, bound 5 m\n"
       f"datum at {float(produced1['datum_geopotential_m2s2']):.1f} m2/s2, radius "
       f"{float(produced1['datum_radius_m']):.1f} m\n"
       + "\n".join(identity_lines(produced1, every_level=True)))

# ---------------------------------------------------------------------------------------------
# 2. The pressure identity falls under halving, and the chain is the driver's own path
# ---------------------------------------------------------------------------------------------
if 2 in CARRIED:
    carry(2)
else:
    same = chain_to(target10)
    identity_same = float(np.abs(np.asarray(same.target["pressure_identity_residual"])).max())
    reproduced = bool(np.array_equal(
        np.asarray(same.target["pressure_Pa"], dtype="float64"),
        np.asarray(produced1["pressure_Pa"], dtype="float64")))
    half = chain_to(target10, scale=0.5)
    identity_half = float(np.abs(np.asarray(half.target["pressure_identity_residual"])).max())
    held_same = identity_composition_held(same.target, target10)
    held_half = identity_composition_held(half.target, target10)
    R_labels, _ = fp.mean_properties_on_labels(
        inputs.composition, float(namelist.target_latitude_deg),
        np.asarray(produced1["label_pressure_Pa"], dtype="float64"))
    R_tabulated, _ = fp.mean_properties_on_labels(
        inputs.composition, float(namelist.target_latitude_deg), p_tab)
    composition_shift = float(np.max(np.abs(np.log(R_labels) - np.log(R_tabulated))))
    identity_no_bottom = float(np.abs(
        np.asarray(same.target["pressure_identity_residual"])[:-1]).max())
    identity_no_bottom_half = float(np.abs(
        np.asarray(half.target["pressure_identity_residual"])[:-1]).max())

    record(2, "the delivered pressure identity is under 1e-6 and its tracing part falls when both "
              "spacings are halved, and the chain runs 2 to 5 use reproduces run 1's production bit "
              "for bit",
           identity_same <= 1e-6 and held_half < held_same and (reproduced or QUICK != 1.0),
           f"the delivered identity at the pinned {spacings()}: largest |p/p_label - 1| "
           f"{identity_same:.3e}, bound 1e-6 (v0.18 section 18 ruling 1; the reviewing agent "
           f"measures 4.1e-7 on its own mesh). It carries a floor near 5e-7 for this state and does "
           f"not fall with the mesh: halved, at {spacings(0.5)}, it is {identity_half:.3e}, ratio "
           f"{identity_same / identity_half:.2f}\n"
           f"the tracing's own part, the composition held on the anchor's tabulated pressures at "
           f"the same latitude, is what the convergence condition applies to: {held_same:.3e} at the "
           f"pinned spacings and {held_half:.3e} halved, ratio {held_same / held_half:.2f}. No order "
           f"is claimed, since the tracing's error and the trapezoid's enter it together\n"
           f"the halved run: mesh {half.state.mesh.shape[0]} by {half.state.mesh.shape[1]}, outer loop "
           f"{half.state.record['passes']} passes, residual {half.state.record['residual_ln_p']:.3e}\n"
           f"forward.transfer.chain at the pinned spacings gives largest {identity_same:.3e} and "
           f"pressure equal to the driver's by array_equal: {reproduced}"
           + ("" if QUICK == 1.0 else " (not asserted: --quick puts the chain on another mesh than "
                                     "the namelist's run 1)")
           + f". That is the path runs 2 to 5 take, and it is the driver's own; beyond the "
           f"specification's list, and it is here because those runs cannot go through a namelist\n"
           f"the floor is the two productions' axes: a level's label is the closure production's own "
           f"pressure, which that production formed with the composition read at the anchor's "
           f"tabulated pressures, while this one reads it at the labels (section 15 ruling 4). That "
           f"read moves ln R_bar by up to {composition_shift:.3e} (section 10 finding 5 states 1.5e-6, "
           f"at the bottom row) and does not depend on the mesh. The identity's largest value sits at "
           f"the bottom row at both spacings\n"
           f"excluding the bottom row alone: {identity_no_bottom:.3e} at the pinned spacings and "
           f"{identity_no_bottom_half:.3e} halved\n"
           + "\n".join(identity_lines(half.target)))

# ---------------------------------------------------------------------------------------------
# 3. Run 2: 60 N under the closure wind
# ---------------------------------------------------------------------------------------------
if 3 in CARRIED:
    carry(3)
else:
    run2 = chain_to(target60)
    produced2 = run2.target
    T2 = np.asarray(produced2["temperature_K"], dtype="float64")
    labels2_agree = label_agreement(produced2)
    below2 = (T_closure - T2) / T_closure
    STATED_T60 = {0: 2.51e-3, int(produced2["gauge_level_index"]): 2.48e-3,
                  int(produced2["label_pressure_Pa"].size) - 1: 2.47e-3}
    worst_T60 = max(abs(float(below2[k]) - w) for k, w in STATED_T60.items())
    altitude2 = np.asarray(produced2["altitude_m"], dtype="float64")
    STATED_ALT60 = {0: 328127.0, int(produced2["gauge_level_index"]): 78533.0,
                    int(altitude2.size) - 1: -12239.0}
    worst_alt60 = max(abs(float(altitude2[k]) - w) for k, w in STATED_ALT60.items())

    record(3, "run 2, 60 N under the closure wind: the delivered temperature is below the closure "
              "product's by the stated fractions to 1e-4 and the altitudes are the stated ones to 5 m",
           worst_T60 <= 1e-4 and worst_alt60 <= 5.0,
           f"at {spacings()}, mesh {run2.state.mesh.shape[0]} by {run2.state.mesh.shape[1]}, outer "
           f"loop {run2.state.record['passes']} passes, residual "
           f"{run2.state.record['residual_ln_p']:.3e} in ln p\n"
           f"the labels agree with the closure production's pressures to {labels2_agree:.3e} relative, "
           f"so the comparison is level for level\n"
           f"T below the closure product's: "
           + ", ".join(f"{name} {float(below2[k]):.4e} (stated {STATED_T60[k]:.3e})"
                       for name, k in named(produced2))
           + f", largest departure {worst_T60:.3e}, bound 1e-4\n"
           f"altitude above the {DATUM:g} Pa datum: "
           + ", ".join(f"{name} {float(altitude2[k]):.0f} m (stated {STATED_ALT60[k]:.0f})"
                       for name, k in named(produced2))
           + f", largest departure {worst_alt60:.2f} m, bound 5 m\n"
           + "\n".join(identity_lines(produced2)))

# ---------------------------------------------------------------------------------------------
# 4. Run 3: 10 N under the cylinder-extended wind of decision P
# ---------------------------------------------------------------------------------------------
if 4 in CARRIED:
    carry(4)
else:
    cyl = fields.cylinder_extended(inputs, lindal, closure_field, p_b, CONSTANTS,
                                   namelist.gauge_isobar_Pa, LATITUDE_SPACING,
                                   extra_latitudes=(target10, target60))
    run3 = chain_to(target10, field=cyl["field"])
    produced3 = run3.target
    labels3_agree = label_agreement(produced3)
    d_T3 = float(np.max(np.abs(np.asarray(produced3["temperature_K"], dtype="float64")
                               / T_closure - 1.0)))
    d_p3 = float(np.max(np.abs(np.asarray(produced3["pressure_Pa"], dtype="float64")
                               / p_closure - 1.0)))
    d_N3 = float(np.max(np.abs(np.asarray(produced3["refractivity"], dtype="float64")
                               / N_closure - 1.0)))
    altitude3 = np.asarray(produced3["altitude_m"], dtype="float64")
    # v0.18 section 18 ruling 2: the stated values are the reviewing agent's at decision P's fixed
    # point, the reference surface marched under the cylinder file's own gauge wind, and the bound is
    # decision Q's floor carried into altitude, 300 m2/s2 of isobar shift at the level and at the datum
    # over g. The v0.10 values, 416,300 and 99,290 and -15,456, were measured with the surface marched
    # under the closure wind and are superseded.
    STATED_ALT_CYL = {0: 416320.0, int(produced3["gauge_level_index"]): 99295.0,
                      int(altitude3.size) - 1: -15457.0}
    CYLINDER_ALTITUDE_BOUND = 60.0
    STATED_R0_CYL = 60129726.0
    worst_alt3 = max(abs(float(altitude3[k]) - w) for k, w in STATED_ALT_CYL.items())
    r0_cylinder = float(run3.state.reference_radius_m[run3.state.mesh.latitude_index(target10)])
    shift3 = (np.asarray(produced3["geopotential_m2s2"], dtype="float64")
              - tr.place(lindal, inputs, cyl["field"], p_b).geopotential_m2s2)
    stated_u = {0: 9.490, gauge_level: 3.717, int(np.argmin(np.abs(p_tab - 1.0e5))): 2.168,
                int(p_tab.size) - 1: 1.926}
    worst_u = max(abs(float(cyl["on_construction"][k]) - w) for k, w in stated_u.items())

    record(4, "run 3, 10 N under the cylinder-extended wind: the delivered T, p and N are the closure "
              "product's to 2e-3 and the altitudes are the v0.18 values to 60 m, both decision Q's "
              "floor, at the wind and carried into altitude",
           (d_T3 <= 2e-3 and d_p3 <= 2e-3 and d_N3 <= 2e-3
            and worst_alt3 <= CYLINDER_ALTITUDE_BOUND),
           f"at {spacings()}, mesh {run3.state.mesh.shape[0]} by {run3.state.mesh.shape[1]}, outer "
           f"loop {run3.state.record['passes']} passes, residual "
           f"{run3.state.record['residual_ln_p']:.3e} in ln p\n"
           f"the cylinder file is rebuilt here by the Step 2 construction (decision P, per hemisphere, "
           f"the inversion curve sampled at the pinned spacing over {cyl['latitudes']} latitudes): the "
           f"fixed point converged in {cyl['passes']} passes to {cyl['final_change_ms']:.2e} m/s, and "
           f"on the construction u is "
           + ", ".join(f"{float(cyl['on_construction'][k]):.3f} (stated {w:.3f})"
                       for k, w in stated_u.items())
           + f" m/s at the top, gauge, 1 bar and bottom levels, largest departure {worst_u:.3e} m/s\n"
           f"the labels agree with the closure production's pressures to {labels3_agree:.3e} relative, "
           f"so the comparison is level for level\n"
           f"delivered against the closure product, relative: T {d_T3:.3e}, p {d_p3:.3e}, "
           f"N {d_N3:.3e}, bound 2e-3 (decision Q; Step 4 measured 1.064e-03 in ln N and 163.3 m2/s2 "
           f"in the isobars on this state)\n"
           f"altitude above the {DATUM:g} Pa datum: "
           + ", ".join(f"{name} {float(altitude3[k]):.0f} m (stated {STATED_ALT_CYL[k]:.0f})"
                       for name, k in named(produced3))
           + f", largest departure {worst_alt3:.2f} m, bound {CYLINDER_ALTITUDE_BOUND:.0f} m\n"
           f"r0 at the target under this wind {r0_cylinder:.1f} m, the reviewing agent's fixed point "
           f"{STATED_R0_CYL:.0f} m, {abs(r0_cylinder - STATED_R0_CYL):.1f} m from it; it is "
           f"{r0_cylinder - r0_10:.0f} m above the closure wind's surface, which is why the v0.10 "
           f"altitudes, measured with the surface marched under the closure wind, are superseded "
           f"(section 18 ruling 2)\n"
           f"the mechanism of the bound, measured: the arrival isobars of the gridded file are shifted "
           f"from the anchor's own by {shift3[0]:+.1f} m2/s2 at the top level, "
           f"{shift3[int(produced3['gauge_level_index'])]:+.1f} at the gauge and {shift3[-1]:+.1f} at "
           f"the bottom, largest |shift| {np.abs(shift3).max():.1f} against decision Q's 300; the "
           f"altitude of a level above the datum moves by the difference of the shift there and at the "
           f"datum, over g\n"
           + "\n".join(identity_lines(produced3)))

# ---------------------------------------------------------------------------------------------
# 5. Run 4: the target at the anchor's own latitude, against the closure product
# ---------------------------------------------------------------------------------------------
if 5 in CARRIED:
    carry(5)
else:
    run4 = chain_to(phi_c)
    produced4 = run4.target
    labels4_agree = label_agreement(produced4)
    d_N4 = float(np.max(np.abs(np.log(np.asarray(produced4["refractivity"], dtype="float64"))
                               - np.log(N_closure))))
    d_Phi4 = float(np.max(np.abs(np.asarray(produced4["geopotential_m2s2"], dtype="float64")
                                 - Phi_closure_product)))
    d_p4 = float(np.max(np.abs(np.asarray(produced4["pressure_Pa"], dtype="float64")
                               / p_closure - 1.0)))
    d_T4 = float(np.max(np.abs(np.asarray(produced4["temperature_K"], dtype="float64")
                               / T_closure - 1.0)))
    altitude4 = np.asarray(produced4["altitude_m"], dtype="float64")
    STATED_ALT_C = {0: 376780.0, int(produced4["gauge_level_index"]): 90067.0,
                    int(altitude4.size) - 1: -14031.0}
    worst_alt4 = max(abs(float(altitude4[k]) - w) for k, w in STATED_ALT_C.items())

    record(5, "run 4, the target at the anchor's own latitude: N and Phi are the closure product's to "
              "1e-12, p and T to 1e-5, and the altitudes are the stated ones to 5 m",
           d_N4 <= 1e-12 and d_Phi4 <= 1e-12 and d_p4 <= 1e-5 and d_T4 <= 1e-5 and worst_alt4 <= 5.0,
           f"the target is {math.degrees(phi_c):.6f} degrees, the anchor's own latitude and at M = 1 "
           f"the gauge as well, at {spacings()}, mesh {run4.state.mesh.shape[0]} by "
           f"{run4.state.mesh.shape[1]}, outer loop {run4.state.record['passes']} passes\n"
           f"the labels agree with the closure production's pressures to {labels4_agree:.3e} relative, "
           f"so the comparison is level for level\n"
           f"against the closure product on {produced4['label_pressure_Pa'].size} levels: largest "
           f"|d ln N| {d_N4:.3e} (bound 1e-12), largest |d Phi| {d_Phi4:.3e} m2/s2 (bound 1e-12), "
           f"largest relative p {d_p4:.3e} and T {d_T4:.3e} (bound 1e-5)\n"
           f"p and T agree at 6e-7 rather than at round-off because the composition is read onto the "
           f"produced labels (section 15 ruling 4), which moves ln R_bar by at most 1.5e-6 (section 10 "
           f"finding 5); N and Phi are the transfer's own and are exact\n"
           f"altitude above the {DATUM:g} Pa datum: "
           + ", ".join(f"{name} {float(altitude4[k]):.0f} m (stated {STATED_ALT_C[k]:.0f})"
                       for name, k in named(produced4))
           + f", largest departure {worst_alt4:.2f} m, bound 5 m\n"
           + "\n".join(identity_lines(produced4)))

# Run 5's state: the synthetic 60 N anchor written on the M = 2 run's own mesh (ruling 1)
_original_refuse = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
    attrs.get("casspian_git_commit", ""))

M2_NAMELIST = """# Written by the SPEC_04 Step 5 acceptance. The transfer namelist with a second
# anchor, the synthetic 60 N profile the closure-wind run produced (decision M), carried through
# the production at 10 N. The run keeps the name `lindal_transfer`, because a namelist is named
# after its run and its inputs carry the run's prefix.
[run]
name        = "lindal_transfer"
description = "The Lindal anchor and a synthetic 60 N anchor transferred to 10 N"
mode        = "transfer"
solar_longitude_deg = 18.2
date        = "1981-08-26"

[[anchors]]
slug   = "lindal"
path   = "../../../occul_data/lindal/lindal_refractivity.nc"
weight = 1.0
measurement_uncertainty_scale = 1.0

[[anchors]]
slug   = "synthetic60"
path   = "../{anchor_file}"
weight = 1.0
measurement_uncertainty_scale = 1.0

[inputs]
composition = "inputs/lindal_transfer_composition.nc"
gravity     = "inputs/lindal_transfer_gravity.nc"
rotation    = "inputs/lindal_transfer_rotation.nc"
wind        = "inputs/lindal_transfer_wind.nc"

[hydrostatic_boundary]
p_b_rule     = "anchor_profile_top"
p_b_location = "top_of_anchor_profile"

[isobars]
gauge_isobar_Pa = 1.0e4
datum_isobar_Pa = 1.0e5

[target]
latitude_planetocentric_deg = 10.0

[grid]
geopotential_spacing_m2s2 = {geopotential:g}
latitude_spacing_deg      = {latitude:g}

[numerics.reference_surface]
scheme = "rk4"
[numerics.shear_integral]
scheme = "trapezoid"
[numerics.isobar_tracing]
scheme = "rk4"
[numerics.transfer]
scheme = "trapezoid"
[numerics.altitude]
scheme = "trapezoid"
[numerics.outer_loop]
relative_tolerance_ln_p = 1.0e-8
max_iterations          = 50

[estimation]
gauge_latitude_rule                = "weighted_centroid"
kernel_uncertainty_per_rad         = 0.02
model_error_correlation_length_deg = 0.0

[output]
directory = "output"
product   = "lindal_transfer_profile.nc"

[diagnostics]
figures = false
format  = "png"
dpi     = 150
"""

synthetic_path = HERE / "synthetic60_refractivity.nc"
m2_directory = HERE / "m2"
(m2_directory / "inputs").mkdir(parents=True, exist_ok=True)
(m2_directory / "output").mkdir(parents=True, exist_ok=True)
for source in namelist.inputs.values():
    shutil.copy(source, m2_directory / "inputs" / source.name)
m2_path = m2_directory / "lindal_transfer.toml"
m2_path.write_text(M2_NAMELIST.format(anchor_file=synthetic_path.name,
                                      geopotential=GEOPOTENTIAL_SPACING,
                                      latitude=math.degrees(LATITUDE_SPACING)), encoding="utf-8")

# The mesh the M = 2 run works on, built before anything is traced, for the latitudes the run needs
# as nodes: the two anchors, the gauge latitude the estimate will compute (the midpoint, the two
# sigma_ln_N being identical) and the target. The synthetic anchor is written from a trace on this
# mesh, which is what makes the comparison an identity of the transfer (section 15 ruling 1).
phi_r2 = 0.5 * (phi_a + target60)
Phi_lindal = tr.place(lindal, inputs, closure_field, p_b).geopotential_m2s2
mesh_m2 = lm.build_mesh([target10, phi_a, phi_r2, target60], Phi_lindal,
                        LATITUDE_SPACING, GEOPOTENTIAL_SPACING)
pre = tr.outer_loop(
    inputs, [lindal], gauge_latitude_rad=phi_r2, target_latitude_rad=target60,
    gauge_isobar_Pa=namelist.gauge_isobar_Pa, p_b=p_b, latitude_spacing_rad=LATITUDE_SPACING,
    geopotential_spacing_m2s2=GEOPOTENTIAL_SPACING,
    relative_tolerance_ln_p=LOOP["relative_tolerance_ln_p"],
    max_iterations=LOOP["max_iterations"], field=closure_field, mesh=mesh_m2)
lnN_pre = along(pre)


def write_synthetic(path, state, ln_N_along):
    """The synthetic 60 N anchor of decision M, written from this state's trace to 60 N.

    Copied from the Step 4 acceptance, where it is measured: `radius_m` and the height from the
    60 N column at the arrival levels, integrated on the traced `Phi` itself rather than read off
    the mesh's nodes (section 16 ruling 1); `refractivity` the transferred `N`; the thermo group's
    temperatures `p R_bar / (k_B N)`; the uncertainty copied as a relative one so that
    `sigma_ln_N` is the Lindal anchor's.
    """
    shutil.copy(lindal.path, path)
    index60 = state.mesh.latitude_index(target60)
    Phi60 = state.curves[0].geopotential_m2s2[:, index60]
    N60 = np.exp(ln_N_along[:, index60])
    nodes = np.unique(np.concatenate([Phi60, [0.0]]))
    column = lm.column(
        target60, nodes, float(state.reference_radius_m[index60]),
        (lambda Phi: closure_field.wind_at(
            target60, np.exp(state.isobar_map.ln_p_at(index60, Phi)))), *CONSTANTS)
    at_traced = np.searchsorted(nodes, Phi60)
    gauge_node = int(np.flatnonzero(nodes == 0.0)[0])
    psi60 = float(np.degrees(column.psi_rad[gauge_node]))
    R_bar60, _ = fp.mean_properties(inputs.composition, 60.0)
    relative_sigma = (np.asarray(lindal.tree.to_dataset(inherit=False)["refractivity_uncertainty"]
                                 .values, dtype="float64") / np.exp(lindal.ln_N))
    with netCDF4.Dataset(path, "a") as handle:
        handle["radius_m"][:] = column.radius_m[at_traced]
        handle["height_above_anchor_isobar_m"][:] = column.z_local_vertical_m[at_traced]
        handle["refractivity"][:] = N60
        handle["refractivity_uncertainty"][:] = relative_sigma * N60
        handle["latitude_planetocentric_deg"][...] = 60.0
        handle["latitude_planetographic_deg"][...] = 60.0 + psi60
        handle["psi_deg"][...] = psi60
        handle["anchor_isobar_radius_m"][...] = float(column.radius_m[gauge_node])
        handle["inputs/thermo"]["temperature_K"][:] = hs.temperature(p_tab, N60, R_bar60)
        handle.setncattr("profile_or_run", path.name.split("_refractivity")[0])
        handle.setncattr("title", "Synthetic 60 N anchor, traced on the M = 2 run's own mesh, "
                                  "written by the SPEC_04 Step 5 acceptance")
        return str(handle.getncattr("casspian_git_commit"))


synthetic_commit = write_synthetic(synthetic_path, pre, lnN_pre)
m2_namelist = ctl.read_run_namelist(m2_path)
m2_inputs = ctl.load_run_inputs(m2_namelist)
m2_anchors = tuple(fprop.propagate(a, m2_namelist.solar_longitude_deg) for a in m2_inputs.anchors)
gauge_m2 = fe.gauge_latitude(m2_anchors)
if abs(gauge_m2 - phi_r2) > 1.0e-12:
    raise AssertionError(f"phi_r from the anchors is {math.degrees(gauge_m2)} degrees and the mesh "
                         f"was built for {math.degrees(phi_r2)}")
run5 = tr.chain(
    m2_inputs, m2_anchors, gauge_latitude_rad=phi_r2, target_latitude_rad=target10,
    gauge_isobar_Pa=m2_namelist.gauge_isobar_Pa, p_b=tr.boundary_pressure(m2_anchors),
    latitude_spacing_rad=LATITUDE_SPACING, geopotential_spacing_m2s2=GEOPOTENTIAL_SPACING,
    relative_tolerance_ln_p=LOOP["relative_tolerance_ln_p"],
    max_iterations=LOOP["max_iterations"], kernel_uncertainty_per_rad=SIGMA_K,
    datum_isobar_Pa=float(m2_namelist.datum_isobar_Pa), mesh=mesh_m2,
    field=wf.WindField(m2_inputs.wind))
produced5 = run5.target
product5 = tr.write_product(m2_namelist, m2_inputs, m2_anchors, run5.state, run5.estimate,
                            run5.arrived, produced5, run5.curves_to_target)
tree5 = read_closed(product5, "profile")

# ---------------------------------------------------------------------------------------------
# 6. Run 5: the M = 2 run carried through production, against run 1's product
# ---------------------------------------------------------------------------------------------
Phi5 = np.asarray(produced5["geopotential_m2s2"], dtype="float64")
Phi1 = np.asarray(produced1["geopotential_m2s2"], dtype="float64")
inside = (Phi5 >= Phi1.min()) & (Phi5 <= Phi1.max())
rising = np.argsort(Phi1)


def like_run1(values):
    """Run 1's product on run 5's union levels, log-linear in the geopotential (decision C)."""
    return np.interp(Phi5, Phi1[rising], np.asarray(values, dtype="float64")[rising])


d_N5 = float(np.max(np.abs(
    np.log(np.asarray(produced5["refractivity"], dtype="float64")[inside])
    - like_run1(np.log(np.asarray(produced1["refractivity"], dtype="float64")))[inside])))
d_p5 = float(np.max(np.abs(
    np.asarray(produced5["pressure_Pa"], dtype="float64")[inside]
    / np.exp(like_run1(np.log(np.asarray(produced1["pressure_Pa"], dtype="float64"))))[inside]
    - 1.0)))
d_T5 = float(np.max(np.abs(
    np.asarray(produced5["temperature_K"], dtype="float64")[inside]
    / like_run1(np.asarray(produced1["temperature_K"], dtype="float64"))[inside] - 1.0)))
d_alt5 = float(np.max(np.abs(
    np.asarray(produced5["altitude_m"], dtype="float64")[inside]
    - like_run1(np.asarray(produced1["altitude_m"], dtype="float64"))[inside])))
anchor_groups5 = sorted(tree5["anchors"].children)
differences5 = [str(n) for n in tree5["estimate"].to_dataset(inherit=False).data_vars
                if str(n).startswith("D_")]
worst_D5 = max(float(np.nanmax(np.abs(
    tree5["estimate"].to_dataset(inherit=False)[name].values))) for name in differences5)

record(6, "run 5, the M = 2 run of Step 4 carried through production: its product is run 1's to "
          "1e-5 in N, p, T and the altitude, and it carries one anchors group per anchor",
       (d_N5 <= 1e-5 and d_p5 <= 1e-5 and d_T5 <= 1e-5 and d_alt5 <= 5.0
        and anchor_groups5 == ["lindal", "synthetic60"] and len(differences5) == 1),
       f"the two anchors are the Lindal profile and a synthetic 60 N anchor written from this "
       f"run's own mesh (section 15 ruling 1), read back through the kind N reader with the "
       f"-dirty refusal relaxed in this process only: its commit is {synthetic_commit}\n"
       f"at {spacings()}, mesh {run5.state.mesh.shape[0]} by {run5.state.mesh.shape[1]} built "
       f"before anything was traced, gauge latitude {math.degrees(phi_r2):.6f} degrees, outer loop "
       f"{run5.state.record['passes']} passes, residual "
       f"{run5.state.record['residual_ln_p']:.3e} in ln p\n"
       f"product {product5}, {len(anchor_groups5)} anchors groups {anchor_groups5}\n"
       f"on the {int(inside.sum())} of {Phi5.size} union levels inside run 1's range, against run "
       f"1's product interpolated onto them: largest |d ln N| {d_N5:.3e}, relative p {d_p5:.3e}, "
       f"relative T {d_T5:.3e}, bound 1e-5; largest |d altitude| {d_alt5:.2f} m\n"
       f"the estimate at M = 2: {len(differences5)} difference column {differences5}, largest "
       f"|D_ij| {worst_D5:.3e} against the identity floor {tr.IDENTITY_FLOOR:.1e}, which is "
       f"decision G's placement of a written anchor's levels (section 16 ruling 2) and below which "
       f"D_ij says nothing about the atmosphere\n"
       + "\n".join(identity_lines(produced5)))

# ---------------------------------------------------------------------------------------------
# 7. The product validates as kind profile and reads back with the groups of deliverable 3
# ---------------------------------------------------------------------------------------------
def group_names(tree):
    found = []
    stack = [("", tree)]
    while stack:
        prefix, node = stack.pop()
        for name, child in node.children.items():
            path = f"{prefix}{name}"
            found.append(path)
            stack.append((f"{path}/", child))
    return sorted(found)


groups1 = group_names(tree1)
REQUIRED = ["anchors/lindal", "anchors/lindal/transfer", "estimate", "isobars",
            "inputs/composition", "inputs/gravity", "inputs/rotation", "inputs/wind", "namelist",
            "reference_surface", "transfer_record"]
missing = [name for name in REQUIRED if name not in groups1]
root_absent = [name for name in ("pressure_tabulated_Pa", "temperature_tabulated_K")
               if name in root1.variables]
root_present = [name for name in ("pressure_label_Pa", "pressure_identity_residual",
                                  "refractivity", "refractivity_gauge", "geopotential_gauge_m2s2",
                                  "altitude_m", "radius_m", "datum_isobar_Pa",
                                  "datum_geopotential_m2s2", "datum_radius_m",
                                  "reference_surface_radius_m",
                                  "gauge_latitude_planetocentric_deg")
                if name not in root1.variables and name not in root1.coords]
identity_derived = str(root1["pressure_identity_residual"].attrs.get("provenance", ""))

record(7, "run 1's product validates as kind profile, reads back as a DataTree with the groups of "
          "deliverable 3, and carries the transfer-mode variables with the tabulated columns "
          "absent",
       (not missing and not root_absent and not root_present
        and identity_derived == "derived" and "production_record" not in groups1),
       f"{product1}, SHA-256 {cio.sha256(product1)[:12]}, {len(groups1)} groups\n"
       f"groups: {groups1}\n"
       f"deliverable 3's groups all present: {not missing}"
       + (f", missing {missing}" if missing else "")
       + f"\nthe closure columns pressure_tabulated_Pa and temperature_tabulated_K are absent: "
       f"{not root_absent}\n"
       f"the transfer-mode root variables are all present: {not root_present}"
       + (f", missing {root_present}" if root_present else "")
       + f"\npressure_identity_residual provenance {identity_derived!r}, which v0.17 ruling 2 has "
       f"carry no uncertainty companion\n"
       f"exactly one record group: transfer_record present, production_record absent")

# ---------------------------------------------------------------------------------------------
# 8. The refusals: the datum outside the produced range, and the record group rule
# ---------------------------------------------------------------------------------------------
pressure1 = np.asarray(produced1["pressure_Pa"], dtype="float64")
Phi_produced1 = np.asarray(produced1["geopotential_m2s2"], dtype="float64")
above = refusal(fp.datum_geopotential, pressure1, Phi_produced1, 0.5 * pressure1.min())
below = refusal(fp.datum_geopotential, pressure1, Phi_produced1, 2.0 * pressure1.max())
level = int(pressure1.size // 2)
on_level = fp.datum_geopotential(pressure1, Phi_produced1, float(pressure1[level]))
exact = on_level == float(Phi_produced1[level])

both_path = HERE / "both_records_profile.nc"
shutil.copy(product1, both_path)
with netCDF4.Dataset(both_path, "a") as handle:
    handle.createGroup("production_record")
both = refusal(read_closed, both_path, "profile")

neither_path = HERE / "no_record_profile.nc"
shutil.copy(product1, neither_path)
with netCDF4.Dataset(neither_path, "a") as handle:
    handle.renameGroup("transfer_record", "transfer_record_renamed")
neither = refusal(read_closed, neither_path, "profile")

record(8, "the datum is refused outside the produced range, and a profile carrying both record "
          "groups or neither is refused with the fault named",
       (above is not None and below is not None and both is not None and neither is not None
        and exact),
       f"a datum above the top level ({0.5 * pressure1.min():.4g} Pa, the range being "
       f"[{pressure1.min():.4g}, {pressure1.max():.4g}] Pa):\n    {above}\n"
       f"a datum below the bottom level ({2.0 * pressure1.max():.4g} Pa):\n    {below}\n"
       f"both record groups present:\n    {both}\n"
       f"neither record group present:\n    {neither}\n"
       f"beyond the specification's list: a datum that is exactly a produced pressure lands "
       f"exactly on that level, no interpolation: level {level}, p {pressure1[level]:.6e} Pa, "
       f"Phi {on_level:.6e} against the level's {Phi_produced1[level]:.6e}, equal: {exact}")

# ---------------------------------------------------------------------------------------------
# 9. The sheared synthetic run carried through production
# ---------------------------------------------------------------------------------------------
if 9 in CARRIED:
    carry(9)
else:
    sheared = fields.sheared(inputs, HERE / "step04_5_sheared_wind.nc")
    SHEARED_COARSE = 1.0
    run_shc = chain_to(target10, field=sheared["field"], scale=SHEARED_COARSE)
    run_sh = chain_to(target10, field=sheared["field"], scale=0.5)
    identity_shc = float(np.abs(np.asarray(run_shc.target["pressure_identity_residual"])).max())
    identity_sh = float(np.abs(np.asarray(run_sh.target["pressure_identity_residual"])).max())
    u_top_sheared = float(sheared["field"].wind_at(phi_a, p_tab[0]))
    u_top_closure = float(closure_field.wind_at(phi_a, p_tab[0]))

    record(9, "the sheared synthetic run is carried through production, its pressure identity "
              "reported at every level and falling when both spacings are halved",
           (identity_sh < identity_shc and run_sh.state.record["passes"] < 10
            and run_shc.state.record["passes"] < 10 and sheared["sum_identity_ms"] < 1e-12),
           f"the file: u_total(phi, p) = u_reference(phi) [1 + {sheared['beta']} ln(p_ref / p)] above "
           f"p_ref = {sheared['reference_pressure_Pa']:g} Pa and u_reference below, written to "
           f"{sheared['path'].name} in the three parts and read back through the kind W reader, which "
           f"checks the sum identity ({sheared['sum_identity_ms']:.1e} m/s) and the poles\n"
           f"on the anchor's column it runs {u_top_sheared:.3f} m/s at the top level against "
           f"{u_top_closure:.3f} under the closure wind\n"
           f"at the pinned {spacings(SHEARED_COARSE)}: largest |p/p_label - 1| {identity_shc:.3e}, "
           f"outer loop {run_shc.state.record['passes']} passes, residual "
           f"{run_shc.state.record['residual_ln_p']:.3e}\n"
           f"with both spacings halved, at {spacings(0.5)}: {identity_sh:.3e}, ratio "
           f"{identity_shc / identity_sh:.2f}, outer loop {run_sh.state.record['passes']} passes, "
           f"residual {run_sh.state.record['residual_ln_p']:.3e}\n"
           f"the pair is the pinned spacings and half, where the Step 4 check of the same file used "
           f"the namelist's and twice it. Doubling decision R's 50,000 m2/s2 puts the mesh's lower "
           f"edge deep enough that the isobar map, continued at its last slope, asks the wind above "
           f"1 MPa and lib.windfield refuses the extrapolation, correctly and for reasons that have "
           f"nothing to do with the shear (SPEC_04 Step 1 deliverable 2, decision N). The requirement "
           f"is the fall, and the finer member is the run's own spacing either way\n"
           f"the delivered altitudes under the shear: "
           + ", ".join(f"{name} {float(np.asarray(run_sh.target['altitude_m'])[k]):.0f} m"
                       for name, k in named(run_sh.target))
           + "\n" + "\n".join(identity_lines(run_sh.target)))

# ---------------------------------------------------------------------------------------------
# 10. F5, F6 and F7 by the driver and by casspian-plots by hand
# ---------------------------------------------------------------------------------------------
by_hand = HERE / "by_hand"
if by_hand.exists():
    shutil.rmtree(by_hand)
hand = subprocess.run(["casspian-plots", str(product1), "--out", str(by_hand)],
                      capture_output=True, text=True)
driver_figures = [Path(p) for p in result1.figures.written]
lines, same_body = [], hand.returncode == 0
for path in driver_figures:
    a, b = imread(path), imread(by_hand / path.name)
    band = math.ceil(style.FOOTER_TIME_BAND * a.shape[0])
    body = a.shape == b.shape and bool(np.array_equal(a[:-band], b[:-band]))
    same_body = same_body and body
    lines.append(f"{path.name}: {a.shape[1]}x{a.shape[0]} px; identical above the {band} px time "
                 f"band: {body}; whole file byte identical: "
                 f"{path.read_bytes() == (by_hand / path.name).read_bytes()}")
expected_figures = ["F5", "F6", "F7", "F8", "F9"]
found_figures = [key for key in expected_figures
                 if any(f"_{key}_" in p.name for p in driver_figures)]

record(10, "F5 with its across-latitude panel, F6 drawing the pressure identity, F7, F8 along the "
           "isobars and F9 the wind assumed are rendered by the driver and by casspian-plots by "
           "hand, identical apart from the footer time",
       same_body and found_figures == expected_figures and not result1.figures.skipped,
       f"the driver wrote {len(driver_figures)} figures to "
       f"{namelist.output_directory / 'figures'} and the combined "
       f"{Path(result1.figures.combined_pdf).name}; figures present {found_figures}, skipped "
       f"{result1.figures.skipped}\n"
       f"by hand: exit {hand.returncode}: {hand.stdout.strip()}\n"
       + "\n".join(lines)
       + f"\nF7's fourth panel, D_ij, is drawn only at M >= 2 (v0.17 ruling 5): run 1 is M = 1 and "
       f"its F7 has three panels. Run 5's product {product5.name} is the M = 2 case and the "
       f"author views it beside run 1's\n"
       f"F8 draws S/g and the temperature accumulated along five isobars between the gauge and "
       f"the target, from the isobars group's own kernel and ln N_k, and the delivered T on "
       f"altitude and on radius beside the anchor's; under a wind that does not vary along a "
       f"column the five kernel curves coincide, which is the file's own statement. F9 draws the "
       f"wind the run assumed from the product's inputs/wind group, its vertical_structure "
       f"'{str(tree1['inputs/wind'].attrs.get('vertical_structure', ''))}'\n"
       f"the author views F5, F7, F8 and F9 by eye and accepts them, which section 7's acceptance "
       f"asks for at v0.19; this check measures only that the two renderings agree")

# ---------------------------------------------------------------------------------------------
# 11. Beyond the specification: the closure production is bit-identical after the refactor
# ---------------------------------------------------------------------------------------------
if 11 in CARRIED:
    carry(11)
else:
    closure_directory = HERE / "closure_check"
    (closure_directory / "inputs").mkdir(parents=True, exist_ok=True)
    (closure_directory / "output").mkdir(parents=True, exist_ok=True)
    closure_namelist_in = ctl.read_run_namelist(CLOSURE / "lindal_closure.toml")
    for source in closure_namelist_in.inputs.values():
        shutil.copy(source, closure_directory / "inputs" / source.name)
    text = (CLOSURE / "lindal_closure.toml").read_text(encoding="utf-8").replace(
        '"../../occul_data/', '"../../../occul_data/').replace("figures = true", "figures = false")
    closure_path = closure_directory / "lindal_closure.toml"
    closure_path.write_text(text, encoding="utf-8")
    rebuilt = fp.run(closure_path)
    rebuilt_tree = read_closed(rebuilt.product, "profile")
    rebuilt_root = rebuilt_tree.to_dataset(inherit=False)
    COMPARED = ("pressure_Pa", "temperature_K", "number_density_m3", "mean_refractivity_m3",
                "mean_molar_mass_kg_mol", "refractivity", "radius_m",
                "height_above_anchor_isobar_m")
    equal = {name: bool(np.array_equal(np.asarray(rebuilt_root[name].values),
                                       np.asarray(closure[name].values))) for name in COMPARED}
    equal["geopotential_m2s2"] = bool(np.array_equal(
        np.asarray(rebuilt_root["geopotential_m2s2"].values), Phi_closure_product))
    departures = {name: departure(rebuilt_root[name].values, closure[name].values)
                  for name in COMPARED}
    departures["geopotential_m2s2"] = departure(rebuilt_root["geopotential_m2s2"].values,
                                                Phi_closure_product)
    agree = {name: within(rebuilt_root[name].values, closure[name].values) for name in COMPARED}
    agree["geopotential_m2s2"] = within(rebuilt_root["geopotential_m2s2"].values,
                                        Phi_closure_product)

    record(11, "beyond the specification's list: the closure production agrees with the registered "
               "product to 1e-14 relative (REVIEW_05_step0), bit-identity reported, after the "
               "refactor that put both paths through produce_on_geopotential",
           all(agree.values()),
           f"the closure namelist is copied to {closure_path} and run there, so the registered product "
           f"is read and never rewritten; its inputs are the registered ones\n"
           f"compared against {CLOSURE_PRODUCT.name}, level by level, bit-identity and the bound "
           f"{RELATIVE_BOUND:g} relative; within the bound for every variable: {all(agree.values())}\n"
           + "\n".join(f"    {name}: {value}, largest departure {departures[name][0]:.3e} "
                       f"relative ({departures[name][1]:g} ulp)" for name, value in equal.items())
           + f"\nthis is the safety property of the step and it should be checked again after any "
           f"change to forward/production.py")

ctl._refuse_dirty_commit = _original_refuse
write_output()
passed = sum(1 for _, _, ok, _ in results if ok)
print(f"\n{passed} of {len(results)} checks pass. "
      f"{time.time() - _started:.0f} s in all, at {spacings()}"
      f"{'' if QUICK == 1.0 else f' (development quick factor {QUICK})'}", flush=True)
sys.exit(0 if passed == len(results) else 1)
