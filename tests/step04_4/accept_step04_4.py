"""Acceptance checks for SPEC_04 v0.13 Step 4: `forward.transfer`, `forward.estimate`, the outer
loop, and the M = 2 identity test.

Every check prints its measured value, and every check states the spacings it ran at. Checks
beyond the specification's acceptance list are labeled so.

**The runs.** The transfer namelist `forward/lindal_transfer/lindal_transfer.toml` and its four
registered inputs are the state; the anchor is the swept Lindal product. The closure wind is the
file; the cylinder-extended wind of decision P is rebuilt here, by the construction of Step 2,
so that this suite stands alone; the sheared wind of the Step 4 expected values is written here
as a kind W file in its three parts and read back through the reader; two constructed fields, a
kernel that reverses sign in the vertical and the closure wind with its sign reversed, exercise
the two things the tracing must refuse or report.

**What is relaxed.** The synthetic anchor of the M = 2 test is written by this script on a
working tree that is not clean, so its `casspian_git_commit` ends in `-dirty` and
`lib.control` refuses it. The refusal is relaxed in this process only, at the one point
`_refuse_dirty_commit`, and the relaxation is named in check 8's output. Nothing in the code
changes, and no registered input is touched.

**The development flag.** `--quick <factor>` multiplies both mesh spacings and is a development
aid only: the recorded run takes no argument and every check prints the spacings it used, which
for the record are the namelist's 0.05 degrees and 5000 m2/s2.
"""

import dataclasses
import math
import re
import shutil
import sys
import time
from pathlib import Path

import netCDF4
import numpy as np

from casspian.forward import estimate as fe
from casspian.forward import production as fp
from casspian.forward import propagate as fprop
from casspian.forward import transfer as tr
from casspian.lib import control as ctl
from casspian.lib import geoid
from casspian.lib import hydrostatic as hs
from casspian.lib import io as cio
from casspian.lib import kernel as lk
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf

HERE = Path("reports/step04_4")
TRANSFER = Path("forward/lindal_transfer")
QUICK = 1.0
if "--quick" in sys.argv:
    QUICK = float(sys.argv[sys.argv.index("--quick") + 1])

HERE.mkdir(parents=True, exist_ok=True)
results = []
_started = time.time()
_last = _started


def record(number, description, passed, detail):
    """One check's result, with the wall clock since the check before it.

    The cost is recorded because the author asked what a transfer run costs (REPORT_04_step1
    section 5b): every run here is one outer loop with its columns rebuilt on every pass.
    """
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


@dataclasses.dataclass(frozen=True)
class TargetProfile:
    """What the production at a target latitude gives, top down. Step 5 delivers it as a product."""

    geopotential_m2s2: np.ndarray
    radius_m: np.ndarray
    z_local_vertical_m: np.ndarray
    refractivity: np.ndarray
    pressure_Pa: np.ndarray
    temperature_K: np.ndarray


#: Checks whose runs this pass does not repeat, when `--carry <output.txt>` is given. SPEC_04
#: v0.16 changes the M = 2 construction and nothing else, and REVIEW_04_step4 asks for the M = 2
#: runs only; these five touch neither the synthetic anchor nor the estimate, so their rows are
#: taken from the run named on the command line, unchanged and marked as carried. Checks 1 and 2
#: are rerun although they are unchanged, because checks 8 to 10 are built on their states.
CARRY_FROM = None
CARRIED = frozenset()
if "--carry" in sys.argv:
    CARRY_FROM = Path(sys.argv[sys.argv.index("--carry") + 1])
    CARRIED = frozenset({3, 4, 5, 11})


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
                        f"{CARRY_FROM.name}, the run of 25 September under SPEC_04 v0.15, which "
                        f"REVIEW_04_step4 accepted for every check but 8 and 10. v0.16 changes "
                        f"the M = 2 construction and nothing this check reaches, and the "
                        f"regression's porcelain is the evidence that neither run moved a "
                        f"product (REVIEW_04_step4 order of work: the M = 2 runs only)."))
        print(f"[{'PASS' if results[-1][2] else 'FAIL'}] {number}. {results[-1][1]} [CARRIED]",
              flush=True)
        write_output()
        return
    raise AssertionError(f"{CARRY_FROM} has no row for check {number}")


def read_closed(path, kind):
    """Read a file under its kind, pull it into memory, and release the handle.

    `lib.io.read` returns a handle on the open file; leaving those handles open while the same
    file is opened again segmentation faults in this HDF5 build, which is what the first run of
    this script did. `lib.control` does exactly this wherever it loads a file.
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
# The run, its anchor, and the helpers every check uses
# ---------------------------------------------------------------------------------------------
namelist = ctl.read_run_namelist(TRANSFER / "lindal_transfer.toml")
inputs = ctl.load_run_inputs(namelist)
anchors = tuple(fprop.propagate(a, namelist.solar_longitude_deg) for a in inputs.anchors)
lindal = anchors[0]

# The spacings Steps 0 to 4 were run and accepted at, pinned here rather than read from the
# namelist. Decision R has the production namelist carry 50,000 m2/s2 from Step 5, because the mesh
# is oversampled in geopotential by an order of magnitude for the transfer; SPEC_04 v0.16 records
# Steps 0 to 4 as run at 5,000, and this suite's bounds and every value in REPORT_04_step4 are at
# that. The latitude spacing is the namelist's and unchanged.
ACCEPTED_LATITUDE_SPACING_DEG = 0.05
ACCEPTED_GEOPOTENTIAL_SPACING = 5.0e3
if float(namelist.grid["latitude_spacing_deg"]) != ACCEPTED_LATITUDE_SPACING_DEG:
    raise AssertionError(
        f"the namelist's latitude spacing is {namelist.grid['latitude_spacing_deg']} degrees and "
        f"this suite was accepted at {ACCEPTED_LATITUDE_SPACING_DEG}")
LATITUDE_SPACING = math.radians(ACCEPTED_LATITUDE_SPACING_DEG) * QUICK
GEOPOTENTIAL_SPACING = ACCEPTED_GEOPOTENTIAL_SPACING * QUICK
LOOP = namelist.numerics["outer_loop"]
# The outer loop's tolerance this suite's bounds were set at, pinned here like the spacings rather
# than read from the namelist (SPEC_07 v0.5 closure). The iteration cap is the namelist's.
ACCEPTED_RELATIVE_TOLERANCE_LN_P = 1.0e-8
SIGMA_K = float(namelist.estimation["kernel_uncertainty_per_rad"])

phi_a = math.radians(lindal.latitude_planetocentric_deg)
target10 = math.radians(float(namelist.target_latitude_deg))
target60 = math.radians(60.0)
p_tab = np.asarray(lindal.label_pressure_Pa, dtype="float64")
gauge_level = lindal.gauge_level_index
p_b = tr.boundary_pressure(anchors)
closure_field = wf.WindField(inputs.wind)
Omega = float(inputs.rotation["angular_rate_rad_s"])
CONSTANTS = (Omega, float(inputs.gravity["GM_m3s2"]), np.asarray(inputs.gravity["J"].values),
             np.asarray(inputs.gravity["degree"].values),
             float(inputs.gravity["normalization_radius_m"]))

TOP, BOTTOM = 0, p_tab.size - 1
TEN_MBAR = int(np.argmin(np.abs(p_tab - 1.0e3)))
NAMED = (("the top level", TOP), ("10 mbar", TEN_MBAR), ("the gauge", gauge_level),
         ("the bottom", BOTTOM))


def spacings(latitude_scale=1.0, geopotential_scale=None):
    geopotential_scale = latitude_scale if geopotential_scale is None else geopotential_scale
    return (f"{math.degrees(LATITUDE_SPACING * latitude_scale):.4g} degrees by "
            f"{GEOPOTENTIAL_SPACING * geopotential_scale:.4g} m2/s2")


def loop_to(target_rad, *, run_anchors=None, field=None, scale=1.0, geopotential_scale=None,
            gauge_latitude_rad=None, mesh=None):
    """One converged outer loop, at the namelist's spacings times the scales given."""
    run_anchors = list(anchors if run_anchors is None else run_anchors)
    geopotential_scale = scale if geopotential_scale is None else geopotential_scale
    phi_r = (fe.gauge_latitude(run_anchors) if gauge_latitude_rad is None
             else float(gauge_latitude_rad))
    return tr.outer_loop(
        inputs, run_anchors,
        gauge_latitude_rad=phi_r, target_latitude_rad=target_rad,
        gauge_isobar_Pa=namelist.gauge_isobar_Pa, p_b=tr.boundary_pressure(run_anchors),
        latitude_spacing_rad=LATITUDE_SPACING * scale,
        geopotential_spacing_m2s2=GEOPOTENTIAL_SPACING * geopotential_scale,
        relative_tolerance_ln_p=ACCEPTED_RELATIVE_TOLERANCE_LN_P,
        max_iterations=LOOP["max_iterations"],
        field=closure_field if field is None else field, mesh=mesh)


def along(state, index=0, anchor=None):
    """One anchor's curves and its `ln N` along them, with the run's composition term."""
    curves = state.curves[index]
    labels = state.placements[index].label_pressure_Pa
    latitude_c, slope = lk.composition_term(inputs.composition, labels)
    ln_N = tr.transfer(state.isobar_kernel, curves,
                       (anchors[index] if anchor is None else anchor).ln_N, labels, latitude_c, slope)
    return curves, ln_N


def mean_properties_on(labels, latitude_rad):
    """`(R_bar, m_bar)` at one latitude, on the isobar labels rather than on the file's levels.

    SPEC_04 section 10 finding 5: transfer mode reads the composition log-linearly in pressure
    onto each isobar's label. `forward.production.produce` takes the column level by level, which
    is right when the levels are the anchor's own and wrong when they are the union's, so the
    production at the target takes the composition this way.
    """
    R_bar, m_bar = fp.mean_properties(inputs.composition, math.degrees(float(latitude_rad)))
    levels = np.asarray(inputs.composition.to_dataset(inherit=False)["pressure_Pa"].values,
                        dtype="float64")
    rising = np.argsort(np.log(levels))
    x = np.log(levels)[rising]
    return (np.interp(np.log(labels), x, np.asarray(R_bar, dtype="float64")[rising]),
            np.interp(np.log(labels), x, np.asarray(m_bar, dtype="float64")[rising]))


def produce_at(state, phi, geopotential_m2s2, ln_N, label_pressure_Pa):
    """`p` and `T` on arrival levels at one latitude, top down. Beyond this step.

    SPEC_04 Step 5 puts the production at the target in the library. The M = 2 identity test
    compares `p` and `T` there, and decision H reports the pressure identity at every level of
    every run, so both are formed here: B4, B5 and B6 through `lib.hydrostatic`, on the
    geopotential the tracing gives and with the composition on the isobar labels. The radius and
    the altitude along the local vertical come from the state's own column and are what Step 5
    delivers; they are read off here so that the profile is the one the production describes.
    """
    column = state.columns[state.mesh.latitude_index(phi)]
    order = np.argsort(np.asarray(geopotential_m2s2, dtype="float64"))[::-1]
    Phi = np.asarray(geopotential_m2s2, dtype="float64")[order]
    labels = np.asarray(label_pressure_Pa, dtype="float64")[order]
    N = np.exp(np.asarray(ln_N, dtype="float64")[order])
    radius = np.interp(Phi[::-1], column.geopotential_m2s2, column.radius_m)[::-1]
    z_lv = np.interp(Phi[::-1], column.geopotential_m2s2, column.z_local_vertical_m)[::-1]
    R_bar, m_bar = mean_properties_on(labels, phi)
    pressure = hs.pressure_from_top(p_b, hs.layer_mass(hs.density(N, R_bar, m_bar), Phi))
    return order, labels, TargetProfile(
        geopotential_m2s2=Phi, radius_m=radius, z_local_vertical_m=z_lv, refractivity=N,
        pressure_Pa=pressure, temperature_K=hs.temperature(pressure, N, R_bar))


def crossing(state, phi, values):
    """One column of a field on the curves' nodes, at a latitude that is a mesh node."""
    return values[:, state.mesh.latitude_index(phi)]


def relative(got, want):
    return abs(float(got) - float(want)) / abs(float(want))


print(f"SPEC_04 Step 4 acceptance. Mesh spacings {spacings()}"
      f"{'' if QUICK == 1.0 else f' (development quick factor {QUICK})'}", flush=True)

# ---------------------------------------------------------------------------------------------
# 1. The closure wind to 10 N at the namelist's spacings
# ---------------------------------------------------------------------------------------------
STATED_10 = {"shift": {TOP: -31258.1, gauge_level: 0.0, BOTTOM: 11339.7},
             "d_ln_N": {TOP: 1.1090e-2, TEN_MBAR: 1.0988e-2, gauge_level: 1.0950e-2,
                        BOTTOM: 1.0899e-2}}
STATED_60 = {"shift": {TOP: -7117.0, BOTTOM: 2585.1},
             "d_ln_N": {TOP: 2.5119e-3, TEN_MBAR: 2.4922e-3, gauge_level: 2.4840e-3,
                        BOTTOM: 2.4749e-3}}


def measure(state, target_rad, field=None):
    """The isobar shift and `d ln N` from the anchor to a target, and what they produce there."""
    curves, ln_N = along(state)
    shift = crossing(state, target_rad, curves.geopotential_m2s2) - crossing(
        state, phi_a, curves.geopotential_m2s2)
    d_ln_N = crossing(state, target_rad, ln_N) - crossing(state, phi_a, ln_N)
    return curves, ln_N, shift, d_ln_N


def report_lines(shift, d_ln_N, stated):
    lines = []
    for name, index in NAMED:
        want = stated["d_ln_N"].get(index)
        lines.append(f"    {name:14s} d ln N {d_ln_N[index]:+.6e}"
                     + (f", stated {want:+.4e}, difference {d_ln_N[index] - want:+.2e}"
                        if want is not None else ""))
    for name, index in NAMED:
        want = stated["shift"].get(index)
        if want is None:
            continue
        lines.append(f"    {name:14s} shift  {shift[index]:+12.1f} m2/s2, stated {want:+.1f}"
                     + (f", {100.0 * relative(shift[index], want):.3f} percent"
                        if want != 0.0 else f", difference {shift[index] - want:+.3e}"))
    return lines


state10 = loop_to(target10)
curves10, lnN10, shift10, d10 = measure(state10, target10)
worst_d10 = max(abs(d10[i] - w) for i, w in STATED_10["d_ln_N"].items())
worst_s10 = max(relative(shift10[i], w) for i, w in STATED_10["shift"].items() if w != 0.0)
gauge_shift10 = float(shift10[gauge_level])

# Beyond the acceptance list: the estimate at M = 1 is the identity (deliverable 2), and the
# pressure identity of decision H is reported at every level of every run.
phi_r10 = fe.gauge_latitude(anchors)
arrival10 = fe.Arrival(
    slug=lindal.slug, latitude_rad=phi_a, weight=lindal.weight,
    geopotential_m2s2=crossing(state10, phi_r10, curves10.geopotential_m2s2),
    C_i=crossing(state10, phi_r10, lnN10), sigma_ln_N=fe.sigma_ln_N(lindal),
    label_pressure_Pa=state10.placements[0].label_pressure_Pa)
estimate10 = fe.estimate([arrival10], SIGMA_K, phi_r10)
rising = np.argsort(arrival10.geopotential_m2s2)
identity10 = float(np.max(np.abs(estimate10.C - arrival10.C_i[rising])))
order10, labels10, produced10 = produce_at(
    state10, target10, crossing(state10, target10, curves10.geopotential_m2s2),
    crossing(state10, target10, lnN10), state10.placements[0].label_pressure_Pa)
identity_p10 = float(np.max(np.abs(produced10.pressure_Pa / labels10 - 1.0)))

record(1, "the closure wind to 10 N at the namelist's spacings reproduces d ln N to 1e-4 and the "
          "isobar shift to 1 percent",
       worst_d10 <= 1e-4 and worst_s10 <= 0.01 and gauge_shift10 == 0.0,
       f"spacings {spacings()}, mesh {state10.mesh.shape[0]} latitudes by "
       f"{state10.mesh.shape[1]} geopotential nodes, outer loop {state10.record['passes']} "
       f"passes, residual {state10.record['residual_ln_p']:.3e} in ln p\n"
       + "\n".join(report_lines(shift10, d10, STATED_10))
       + f"\nlargest |d ln N - stated| {worst_d10:.3e}, bound 1e-4; largest shift departure "
       f"{100.0 * worst_s10:.3f} percent, bound 1 percent; the gauge isobar does not move, "
       f"{gauge_shift10!r} exactly, because I = 0 on it at every latitude\n"
       f"beyond the acceptance list, deliverable 2 at M = 1: phi_r is the anchor's own latitude "
       f"({math.degrees(phi_r10):.6f} degrees), the union is its own {estimate10.C.size} levels, "
       f"and C - ln N is {identity10:.1e}, one unit in the last place of ln N; no D exists "
       f"({len(estimate10.differences)} pairs) and the reduced chi-square is absent at every "
       f"level ({bool(np.all(np.isnan(estimate10.chi_square_reduced)))})\n"
       f"beyond the acceptance list, decision H: the pressure identity at the target, the "
       f"produced p against the label each isobar carries, is {identity_p10:.3e} at worst over "
       f"{labels10.size} levels; the production at the target is Step 5's deliverable and this "
       f"is formed here from the state's own column")

# ---------------------------------------------------------------------------------------------
# 2. The closure wind to 60 N
# ---------------------------------------------------------------------------------------------
state60 = loop_to(target60)
curves60, lnN60, shift60, d60 = measure(state60, target60)
worst_d60 = max(abs(d60[i] - w) for i, w in STATED_60["d_ln_N"].items())
worst_s60 = max(relative(shift60[i], w) for i, w in STATED_60["shift"].items())
record(2, "the closure wind to 60 N at the namelist's spacings reproduces d ln N to 1e-4 and the "
          "isobar shift to 1 percent",
       worst_d60 <= 1e-4 and worst_s60 <= 0.01,
       f"spacings {spacings()}, mesh {state60.mesh.shape[0]} latitudes by "
       f"{state60.mesh.shape[1]} geopotential nodes, outer loop {state60.record['passes']} "
       f"passes, residual {state60.record['residual_ln_p']:.3e} in ln p\n"
       + "\n".join(report_lines(shift60, d60, STATED_60))
       + f"\nlargest |d ln N - stated| {worst_d60:.3e}, bound 1e-4; largest shift departure "
       f"{100.0 * worst_s60:.3f} percent, bound 1 percent")

# ---------------------------------------------------------------------------------------------
# 3. Both spacings halved
# ---------------------------------------------------------------------------------------------
if 3 in CARRIED:
    carry(3)
else:
    half10 = loop_to(target10, scale=0.5)
    _, _, shift10h, d10h = measure(half10, target10)
    half60 = loop_to(target60, scale=0.5)
    _, _, shift60h, d60h = measure(half60, target60)
    worst_d10h = max(abs(d10h[i] - w) for i, w in STATED_10["d_ln_N"].items())
    worst_s10h = max(relative(shift10h[i], w) for i, w in STATED_10["shift"].items() if w != 0.0)
    worst_d60h = max(abs(d60h[i] - w) for i, w in STATED_60["d_ln_N"].items())
    worst_s60h = max(relative(shift60h[i], w) for i, w in STATED_60["shift"].items())
    record(3, "with both spacings halved the same values are reproduced to 2e-5 in d ln N and 0.2 "
              "percent in the shift",
           max(worst_d10h, worst_d60h) <= 2e-5 and max(worst_s10h, worst_s60h) <= 0.002,
           f"spacings {spacings(0.5)}, meshes {half10.mesh.shape[0]} by {half10.mesh.shape[1]} to "
           f"10 N and {half60.mesh.shape[0]} by {half60.mesh.shape[1]} to 60 N\n"
           + "\n".join(report_lines(shift10h, d10h, STATED_10))
           + "\n" + "\n".join(report_lines(shift60h, d60h, STATED_60))
           + f"\nto 10 N: largest |d ln N - stated| {worst_d10h:.3e} against {worst_d10:.3e} at the "
           f"namelist's spacings, ratio {worst_d10 / worst_d10h:.2f}; shift "
           f"{100.0 * worst_s10h:.4f} percent against {100.0 * worst_s10:.4f}, ratio "
           f"{worst_s10 / worst_s10h:.2f}\n"
           f"to 60 N: largest |d ln N - stated| {worst_d60h:.3e} against {worst_d60:.3e}, ratio "
           f"{worst_d60 / worst_d60h:.2f}; shift {100.0 * worst_s60h:.4f} percent against "
           f"{100.0 * worst_s60:.4f}, ratio {worst_s60 / worst_s60h:.2f}\n"
           f"the ratios are reported and no order is claimed: the stated values are themselves a "
           f"measurement, so what falls here is the departure from them")

# ---------------------------------------------------------------------------------------------
# 4. The cylinder-extended wind (decision P), built here by the Step 2 construction
# ---------------------------------------------------------------------------------------------
if 4 in CARRIED:
    carry(4)
else:
    lat_file = np.radians(np.asarray(inputs.wind["latitude_planetocentric_deg"].values,
                                     dtype="float64"))
    p_file = np.asarray(inputs.wind["pressure_Pa"].values, dtype="float64")
    u_reference_file = np.asarray(inputs.wind["u_reference_ms"].values, dtype="float64")
    p_ref = float(inputs.wind["reference_level_pressure_Pa"])
    ref_col = int(np.flatnonzero(p_file == p_ref)[0])
    profile_a = fp.profile_from_anchor(lindal.tree)
    Phi_closure = tr.place(lindal, inputs, closure_field, p_b).geopotential_m2s2

    _x = np.log(p_tab)[np.argsort(np.log(p_tab))]
    _y = Phi_closure[np.argsort(np.log(p_tab))]
    _flat = tr.IsobarMap(geopotential_m2s2=np.sort(Phi_closure)[None, :],
                         ln_pressure=np.log(p_tab)[np.argsort(Phi_closure)][None, :])


    def Phi_of_p(p):
        """The anchor's flat-isobar map, linear in `ln p`, continued at the slope of its last
        interval, which is the rule `forward.transfer.IsobarMap` applies in the other direction."""
        lp = np.log(np.asarray(p, dtype="float64"))
        low = (_y[1] - _y[0]) / (_x[1] - _x[0])
        high = (_y[-1] - _y[-2]) / (_x[-1] - _x[-2])
        return np.where(lp < _x[0], _y[0] + low * (lp - _x[0]),
                        np.where(lp > _x[-1], _y[-1] + high * (lp - _x[-1]), np.interp(lp, _x, _y)))


    def on_file_grid(Phi):
        """`p(Phi)` of the flat map, held inside the wind file's own pressure range.

        The round trip through the map is not exact to the bit at its ends, so the pressure can land
        one ulp outside the grid and `lib.windfield` refuses it, correctly. The excursion is asserted
        to be at that scale, so this is a clamp of one ulp and nothing structural.
        """
        p = np.exp(_flat.ln_p_at(0, Phi))
        low, high = float(p_file[0]), float(p_file[-1])
        excursion = max(float(np.max(low - np.minimum(p, low))) / low,
                        float(np.max(np.maximum(p, high) - high)) / high)
        assert excursion < 1e-12, f"the flat map leaves the file's grid by {excursion:.3e} relative"
        return np.clip(p, low, high)


    fine = np.arange(float(lat_file[0]), float(lat_file[-1]) + 0.5 * LATITUDE_SPACING,
                     LATITUDE_SPACING)
    lat_inv = np.unique(np.clip(np.concatenate([fine, lat_file, [phi_a, target10, target60]]),
                                float(lat_file[0]), float(lat_file[-1])))
    on_file = np.searchsorted(lat_inv, lat_file)

    Phi_k, field_now, cylinder_passes = Phi_closure.copy(), closure_field, 0
    while True:
        Phi_p = Phi_of_p(p_file)
        nodes = np.unique(np.concatenate([Phi_p, Phi_k, [0.0]]))
        at_pressure = np.searchsorted(nodes, Phi_p)
        r0_inv = geoid.through_anchor(
            lat_inv, phi_a, float(np.asarray(lindal.radius_m)[gauge_level]),
            lambda x: field_now.wind_at(np.asarray(x, dtype="float64"), namelist.gauge_isobar_Pa),
            *CONSTANTS).radius
        at_level = np.searchsorted(nodes, Phi_k)
        radius_inv = np.empty((lat_inv.size, p_file.size))
        anchor_column = None
        for i, latitude in enumerate(lat_inv):
            built_column = lm.column(
                latitude, nodes, float(r0_inv[i]),
                (lambda Phi, la=latitude: field_now.wind_at(la, on_file_grid(Phi))),
                *CONSTANTS)
            radius_inv[i] = built_column.radius_m[at_pressure]
            if latitude == phi_a:
                anchor_column = built_column
        s_reference = radius_inv[:, ref_col] * np.cos(lat_inv)

        def branch(mask):
            s, q = s_reference[mask], lat_inv[mask]
            rising_s = np.argsort(s)
            return s[rising_s], q[rising_s]

        s_north, phi_north = branch(lat_inv >= 0.0)
        s_south, phi_south = branch(lat_inv <= 0.0)
        s_grid = radius_inv[on_file] * np.cos(lat_file)[:, None]
        in_north = lat_file >= 0.0
        phi_star = np.empty_like(s_grid)
        phi_star[in_north] = np.interp(s_grid[in_north], s_north, phi_north)
        phi_star[~in_north] = np.interp(s_grid[~in_north], s_south, phi_south)
        u_cylinder = np.interp(phi_star, lat_file, u_reference_file)
        u_cylinder[0], u_cylinder[-1] = 0.0, 0.0
        # The wind on the construction itself, at the anchor's own levels: the anchor column's radius
        # there turned into a cylinder radius and inverted, with no file grid in between. This is the
        # quantity decision P's stated values are, and it is what the Step 2 acceptance compares with
        # them; reading the written file back at the same points is a different quantity, and the
        # difference between the two is the file grid's own truncation.
        on_construction = np.interp(
            np.interp(anchor_column.radius_m[at_level] * math.cos(phi_a), s_north, phi_north),
            lat_file, u_reference_file)
        change = float(np.max(np.abs(u_cylinder - field_now.u_total_ms)))
        cylinder_passes += 1
        built = inputs.wind.copy(deep=True)
        built["u_total_ms"] = (("latitude_planetocentric", "pressure"), u_cylinder,
                               dict(inputs.wind["u_total_ms"].attrs))
        field_now = wf.WindField(built)
        Phi_k = tr.place(lindal, inputs, field_now, p_b).geopotential_m2s2
        if change < 1.0e-6 or cylinder_passes >= 10:
            break

    cylinder_field = field_now
    u_on_column = np.asarray(cylinder_field.wind_at(np.full(p_tab.shape, phi_a), p_tab),
                             dtype="float64")
    one_bar = int(np.argmin(np.abs(p_tab - 1.0e5)))
    stated_column = {TOP: 9.490, gauge_level: 3.717, one_bar: 2.168, BOTTOM: 1.926}
    worst_column = max(abs(float(on_construction[i]) - w) for i, w in stated_column.items())
    worst_read_back = max(abs(u_on_column[i] - float(on_construction[i]))
                          for i in stated_column)

    state_cyl = loop_to(target10, field=cylinder_field)
    curves_cyl, lnN_cyl = along(state_cyl)
    span = ((state_cyl.mesh.latitude_rad >= target10) & (state_cyl.mesh.latitude_rad <= phi_a))
    Phi_anchor = crossing(state_cyl, phi_a, curves_cyl.geopotential_m2s2)
    lnN_anchor = crossing(state_cyl, phi_a, lnN_cyl)
    shift_cyl = np.abs(curves_cyl.geopotential_m2s2[:, span] - Phi_anchor[:, None])
    d_cyl = np.abs(lnN_cyl[:, span] - lnN_anchor[:, None])
    worst_shift_cyl = float(np.max(shift_cyl))
    worst_d_cyl = float(np.max(d_cyl))
    record(4, "under the cylinder-extended wind the traced isobars and the transferred ln N are the "
              "anchor's own to 300 m2/s2 and 2e-3 at every node (SPEC_04 v0.12 decision Q)",
           worst_shift_cyl <= 300.0 and worst_d_cyl <= 2e-3 and worst_column <= 1e-2,
           f"spacings {spacings()}, mesh {state_cyl.mesh.shape[0]} by {state_cyl.mesh.shape[1]}, "
           f"outer loop {state_cyl.record['passes']} passes, residual "
           f"{state_cyl.record['residual_ln_p']:.3e} in ln p\n"
           f"the cylinder file is rebuilt here by the Step 2 construction (decision P, per "
           f"hemisphere, the inversion curve sampled at the mesh spacing over {lat_inv.size} "
           f"latitudes): the fixed point converged in {cylinder_passes} passes, and on the "
           f"construction, at the anchor's own levels with no file grid in between, u is "
           + ", ".join(f"{float(on_construction[i]):.3f} (stated {w:.3f})"
                       for i, w in stated_column.items())
           + f" m/s at the top, gauge, 1 bar and bottom levels, largest departure "
           f"{worst_column:.3e} m/s, bound 1e-2, which is the quantity the Step 2 acceptance "
           f"compares with decision P\n"
           f"    read back through the written file at the same points it is "
           + ", ".join(f"{u_on_column[i]:.3f}" for i in stated_column)
           + f" m/s, up to {worst_read_back:.3e} m/s from the construction. That difference is the "
           f"file grid's own: a field constant on cylinders, stored on {lat_file.size} latitudes by "
           f"{p_file.size} pressures and read back by decision L's interpolant, is not constant on "
           f"cylinders any more. It is decision Q's truncation measured at the wind rather than at "
           f"the kernel, and it is reported and not bounded\n"
           f"over the {int(np.sum(span))} latitude nodes from the anchor to 10 N and all "
           f"{p_tab.size} levels: largest |Phi_k(phi) - Phi_k| {worst_shift_cyl:.1f} m2/s2, bound "
           f"300; largest |ln N_k(phi) - ln N_k| {worst_d_cyl:.3e}, bound 2e-3\n"
           f"the Step 3 line integrals of the same state estimated 130.3 m2/s2 and 4.746e-04; the "
           f"traced values are {worst_shift_cyl / 130.3:.2f} and {worst_d_cyl / 4.746e-04:.2f} times "
           f"them. The two kernels' floors together are 130 to 224 m2/s2 and 4.7e-4 to 1.1e-3 "
           f"(SPEC_04 v0.13), and the traced values sit inside that class")

# ---------------------------------------------------------------------------------------------
# 5. The crossing refusal, on constructed inputs
# ---------------------------------------------------------------------------------------------
small = lm.build_mesh([target10, math.radians(30.0)], np.array([-2.0e5, 0.0, 2.0e5]),
                      math.radians(5.0), 5.0e4)
wave = 4.0e6 * np.sin(2.0 * np.pi * small.geopotential_m2s2 / 2.0e5)
I_wave = np.tile(wave, (small.latitude_rad.size, 1))
crossed = refusal(tr.trace, small, I_wave, np.array([-2.5e4, 2.5e4]), math.radians(30.0),
                  target10)
I_flat = np.zeros(small.shape)
coincident = refusal(tr.trace, small, I_flat, np.array([2.5e4, 2.5e4]), math.radians(30.0),
                     target10)
if 5 in CARRIED:
    carry(5)
else:
    def stays_ordered(curves):
        """Whether one anchor's curves keep their order in geopotential at every node."""
        Phi = curves.geopotential_m2s2
        return bool(np.all(np.diff(Phi[np.argsort(Phi[:, 0])], axis=0) > 0.0))


    ordered = [stays_ordered(c) for c in (curves10, curves60, curves_cyl)]
    record(5, "a pair of curves that cross is refused, and the runs' own curves stay ordered",
           crossed is not None and "cross at" in crossed
           and coincident is not None and "cross at" in coincident and all(ordered),
           f"constructed kernel: I = 4e6 sin(2 pi Phi / 2e5) per radian on a "
           f"{small.shape[0]} by {small.shape[1]} mesh at 5 degrees by 5e4 m2/s2, two levels at "
           f"-2.5e4 and +2.5e4 m2/s2. The kernel reverses sign twice between the levels and the "
           f"scheme's step carries the lower curve past the upper one:\n"
           f"    {crossed}\n"
           f"constructed anchor: two levels at the same geopotential, which is a profile with a "
           f"layer of no thickness; the refusal fires at the starting node:\n"
           f"    {coincident}\n"
           f"the closure runs to 10 N and 60 N and the cylinder run keep their {p_tab.size} curves "
           f"strictly ordered at every node: {ordered}")

# ---------------------------------------------------------------------------------------------
# 6. A curve reaching the mesh edge: the side and the excess, and the extension
# ---------------------------------------------------------------------------------------------
ran_out = tr.trace(small, np.full(small.shape, 3.0e6), np.array([0.0]),
                   math.radians(30.0), target10)
other_side = tr.trace(small, np.full(small.shape, -3.0e6), np.array([0.0]),
                      math.radians(30.0), target10)
# The caller extends and repeats until the curve stays inside, which is what the outer loop
# does; a constant kernel climbs at every node, so one extension covers one node's overshoot.
extended, after, by_hand = small, ran_out, 0
while after.reached is not None and by_hand < 20:
    extended = extended.extend(after.reached[0], after.reached[1])
    after = tr.trace(extended, np.full(extended.shape, 3.0e6), np.array([0.0]),
                     math.radians(30.0), target10)
    by_hand += 1

reversed_total = inputs.wind.copy(deep=True)
reversed_total["u_total_ms"] = (("latitude_planetocentric", "pressure"),
                                -np.asarray(inputs.wind["u_total_ms"].values, dtype="float64"),
                                dict(inputs.wind["u_total_ms"].attrs))
reversed_field = wf.WindField(reversed_total)
REVERSED_SCALE = 10.0
state_rev = loop_to(target10, field=reversed_field, scale=REVERSED_SCALE,
                    geopotential_scale=1.0)
extensions = [dict(e) for e in state_rev.record["extensions"]]
record(6, "a curve that reaches the mesh's geopotential edge stops the trace and returns the "
          "side and the excess, and the outer loop extends the mesh by it and repeats, recorded",
       (ran_out.reached is not None and after.reached is None and len(extensions) > 0
        and other_side.reached is not None and other_side.reached[0] == "below"),
       f"constructed kernel: I = 3e6 per radian everywhere on the same "
       f"{small.shape[0]} by {small.shape[1]} mesh, one level at Phi = 0. The trace stops at "
       f"{ran_out.latitude_rad.size} of {small.latitude_rad.size} nodes and reports side "
       f"{ran_out.reached[0]!r}, excess {ran_out.reached[1]:.1f} m2/s2; extending by what each "
       f"stop reports and repeating, as the loop does, takes {by_hand} extensions to carry the "
       f"curve across, the mesh growing from {small.geopotential_m2s2.size} to "
       f"{extended.geopotential_m2s2.size} geopotential nodes, and it arrives at "
       f"{after.geopotential_m2s2[0, 0]:+.1f} m2/s2\n"
       f"with the sign of the kernel reversed the curve leaves the other edge: side "
       f"{other_side.reached[0]!r}, excess {other_side.reached[1]:.1f} m2/s2. Both sides are "
       f"exercised because the two are separate tests in the code and the first filing of this "
       f"script had the lower one written with the least overshoot instead of the largest, which "
       f"cannot fire\n"
       f"in the outer loop: the closure wind with u_total reversed in sign, which turns the "
       f"isobar shift outward instead of inward, at {spacings(REVERSED_SCALE, 1.0)} (the "
       f"latitude spacing coarse so that the case is a mechanism and not a value, the "
       f"geopotential spacing the namelist's so that the margin is the run's); the loop recorded "
       f"{len(extensions)} extension(s) "
       f"{extensions} and converged in {state_rev.record['passes']} passes on a mesh of "
       f"{state_rev.mesh.shape[0]} by {state_rev.mesh.shape[1]}\n"
       f"the closure and cylinder runs need none: their isobars move inward, the top level down "
       f"and the bottom level up, so the anchor's own levels stay the outermost")

# ---------------------------------------------------------------------------------------------
# 7. The outer loop's pass count for every wind file
# ---------------------------------------------------------------------------------------------
def pass_line(name, state):
    history = ", ".join(f"{h['residual_ln_p']:.3e}" for h in state.record["history"])
    return (f"    {name:34s} {state.record['passes']} passes, residuals in ln p {history}, "
            f"tolerance {state.record['relative_tolerance_ln_p']:.1e}")


# The runs a carried check made are not in this pass, so their rows are in that check's own
# carried block and not here.
loop_report = [pass_line("closure wind to 10 N", state10),
               pass_line("closure wind to 60 N", state60)]
loop_states = [state10, state60]
if 3 not in CARRIED:
    loop_report += [pass_line("closure wind to 10 N, halved", half10),
                    pass_line("closure wind to 60 N, halved", half60)]
    loop_states += [half10, half60]
if 4 not in CARRIED:
    loop_report.append(pass_line("cylinder wind to 10 N", state_cyl))
    loop_states.append(state_cyl)
loop_report.append(pass_line("reversed closure wind to 10 N", state_rev))
loop_states.append(state_rev)

# ---------------------------------------------------------------------------------------------
# 8. M = 2: the synthetic anchor at 60 N, the identity test
# ---------------------------------------------------------------------------------------------
# SPEC_04 v0.15 Step 4, section 15 ruling 1: the synthetic anchor is written from the Lindal
# anchor's trace to 60 N **on the M = 2 run's own mesh**, so that the check is an identity of the
# transfer and not a comparison of two lattices. The mesh is therefore built first, for the M = 2
# run's required latitudes, and every M = 2 run is formed on it.
#
# It carries both gauge latitudes as nodes: the M = 2 run's, which decision 11 makes the midpoint
# of the two anchors' latitudes because their `sigma_ln_N` is identical, and the one the season
# variant of check 10 moves to, which is computable before anything is traced because the season
# column is the script's own and both anchors share the measurement column. Three runs that are
# compared level for level have to be on one mesh; giving check 10 its own would put back a
# smaller copy of the lattice difference ruling 1 removes.
sigma_shared = fe.sigma_ln_N(lindal)
SEASON_TERM = 2.0e-3
weight_plain = 1.0 / float(np.mean(sigma_shared ** 2))
weight_seasoned = 1.0 / float(np.mean(sigma_shared ** 2 + SEASON_TERM ** 2))
phi_r2 = 0.5 * (phi_a + target60)
phi_r_season = ((weight_plain * phi_a + weight_seasoned * target60)
                / (weight_plain + weight_seasoned))
Phi_lindal = tr.place(lindal, inputs, closure_field, p_b).geopotential_m2s2
mesh_m2 = lm.build_mesh([target10, phi_a, phi_r2, phi_r_season, target60], Phi_lindal,
                        LATITUDE_SPACING, GEOPOTENTIAL_SPACING)

# The Lindal anchor alone on that mesh, which is what the anchor is written from.
state_pre = loop_to(target60, run_anchors=[lindal], gauge_latitude_rad=phi_r2, mesh=mesh_m2)
curves_pre, lnN_pre = along(state_pre)


def write_synthetic(path, state, curves, ln_N_along, note):
    """A synthetic 60 N anchor written from one state's trace of the Lindal anchor to 60 N.

    Decision M: `radius_m` and `height_above_anchor_isobar_m` from the 60 N column at the arrival
    levels, `refractivity` the transferred `N`, the thermo group's pressures the Lindal anchor's
    tabulated ones and its temperatures `p R_bar / (k_B N)`, the anchor radius and `psi` the
    column's at the gauge, the companions the Lindal anchor's and its uncertainty copied as a
    relative one so that `sigma_ln_N` is the Lindal anchor's (v0.15, decision 11 in the text).

    The column is integrated on the traced `Phi` itself, not read off the state's column by
    interpolating from the mesh's geopotential nodes (v0.16 ruling 1). That interpolation was the
    smaller of the two parts of the identity's floor, 1.0 m2/s2 in `Phi` at 5e4 spacing and 0.004
    at 5,000; integrating on the traced nodes removes it, and what is left is decision G's
    placement. `lib.mesh.column` takes the nodes it is given, which is what makes this a change of
    nodes and not of method.
    """
    shutil.copy(lindal.path, path)
    index60 = state.mesh.latitude_index(target60)
    Phi60 = curves.geopotential_m2s2[:, index60]
    N60 = np.exp(ln_N_along[:, index60])
    nodes = np.unique(np.concatenate([Phi60, [0.0]]))
    column = lm.column(
        target60, nodes, float(state.reference_radius_m[index60]),
        (lambda Phi: closure_field.wind_at(
            target60, np.exp(state.isobar_map.ln_p_at(index60, Phi)))), *CONSTANTS)
    at_traced = np.searchsorted(nodes, Phi60)
    radius = column.radius_m[at_traced]
    z_lv = column.z_local_vertical_m[at_traced]
    R_bar60, _ = fp.mean_properties(inputs.composition, 60.0)
    gauge_node = int(np.flatnonzero(nodes == 0.0)[0])
    psi60 = float(np.degrees(column.psi_rad[gauge_node]))
    relative_sigma = (np.asarray(lindal.tree.to_dataset(inherit=False)["refractivity_uncertainty"]
                                 .values, dtype="float64") / np.exp(lindal.ln_N))
    with netCDF4.Dataset(path, "a") as handle:
        handle["radius_m"][:] = radius
        handle["height_above_anchor_isobar_m"][:] = z_lv
        handle["refractivity"][:] = N60
        handle["refractivity_uncertainty"][:] = relative_sigma * N60
        handle["latitude_planetocentric_deg"][...] = 60.0
        # phi_g = phi_c + psi at the anchor isobar (lib.latitude), and psi there is the column's.
        handle["latitude_planetographic_deg"][...] = 60.0 + psi60
        handle["psi_deg"][...] = psi60
        handle["anchor_isobar_radius_m"][...] = float(column.radius_m[gauge_node])
        handle["inputs/thermo"]["temperature_K"][:] = hs.temperature(p_tab, N60, R_bar60)
        handle.setncattr("profile_or_run", path.name.split("_refractivity")[0])
        handle.setncattr("title", f"Synthetic 60 N anchor, {note}, written by the SPEC_04 "
                                  "Step 4 acceptance")
        return (str(handle.getncattr("casspian_git_commit")), Phi60, N60, psi60,
                float(column.radius_m[gauge_node]))


synthetic_path = HERE / "synthetic60_refractivity.nc"
synthetic_commit, Phi60_own, N60_own, psi60, anchor_radius60 = write_synthetic(
    synthetic_path, state_pre, curves_pre, lnN_pre, "traced on the M = 2 run's own mesh")

# The variant of ruling 1, reported and not bounded: the same anchor written from the separate
# M = 1 run to 60 N of check 2, whose mesh is uniform from one spacing below the Lindal latitude.
lattice_path = HERE / "synthetic60lattice_refractivity.nc"
_, Phi60_other, N60_other, _, _ = write_synthetic(
    lattice_path, state60, curves60, lnN60, "traced on a separate run's mesh")
lattice_lnN = float(np.max(np.abs(np.log(N60_own) - np.log(N60_other))))
lattice_Phi = float(np.max(np.abs(Phi60_own - Phi60_other)))

_original_refuse = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
    attrs.get("casspian_git_commit", ""))

M2_NAMELIST = """# Written by the SPEC_04 Step 4 acceptance. The transfer namelist with a second
# anchor, the synthetic 60 N profile the closure-wind run produced (decision M). The run keeps
# the name `lindal_transfer`, because a namelist is named after its run and its inputs carry the
# run's prefix, and the two weights are two run directories beside each other.
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
slug   = "{slug}"
path   = "../{anchor_file}"
weight = {weight}
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


def m2_namelist(weight, anchor_path=None, slug="synthetic60"):
    """One M = 2 run directory: its own copy of the four inputs and a namelist.

    A namelist is named after its run (SPEC_00 section 7.2) and a run's inputs carry the run's
    prefix, so each case is a run directory of its own rather than a second file in one.
    """
    anchor_path = synthetic_path if anchor_path is None else Path(anchor_path)
    directory = HERE / f"m2_{slug}_weight{int(weight)}"
    (directory / "inputs").mkdir(parents=True, exist_ok=True)
    (directory / "output").mkdir(parents=True, exist_ok=True)
    for source in namelist.inputs.values():
        shutil.copy(source, directory / "inputs" / source.name)
    path = directory / "lindal_transfer.toml"
    path.write_text(M2_NAMELIST.format(weight=f"{float(weight):.1f}", slug=slug,
                                       anchor_file=anchor_path.name,
                                       geopotential=GEOPOTENTIAL_SPACING,
                                       latitude=math.degrees(LATITUDE_SPACING)),
                    encoding="utf-8")
    return ctl.read_run_namelist(path)


def m2_anchors(weight, path=None, slug="synthetic60"):
    loaded = ctl.load_run_inputs(m2_namelist(weight, path, slug))
    return tuple(fprop.propagate(a, 18.2) for a in loaded.anchors)


def estimate_from(state, run_anchors, phi_r):
    """Every anchor traced to `phi_r`, the estimate there, and `C` carried to the target."""
    arrivals = []
    for index, anchor in enumerate(run_anchors):
        curves, ln_N = along(state, index, anchor)
        arrivals.append(fe.Arrival(
            slug=anchor.slug,
            latitude_rad=float(np.radians(anchor.latitude_planetocentric_deg)),
            weight=anchor.weight,
            geopotential_m2s2=crossing(state, phi_r, curves.geopotential_m2s2),
            C_i=crossing(state, phi_r, ln_N), sigma_ln_N=fe.sigma_ln_N(anchor),
            label_pressure_Pa=state.placements[index].label_pressure_Pa))
    E = fe.estimate(arrivals, SIGMA_K, phi_r)
    to_target = tr.trace(state.mesh, state.shear_integral, E.geopotential_m2s2, phi_r,
                         target10)
    latitude_c, slope = lk.composition_term(inputs.composition, E.label_pressure_Pa)
    carried = tr.transfer(state.isobar_kernel, to_target, E.C, E.label_pressure_Pa, latitude_c, slope)
    column = int(np.flatnonzero(to_target.latitude_rad == target10)[0])
    return E, arrivals, to_target.geopotential_m2s2[:, column], carried[:, column]


pair = m2_anchors(1.0)
# The two gauge latitudes the mesh was built for are what the estimate computes from the anchors
# as they arrive; if they were not, the mesh would be the wrong one and the runs would refuse.
gauge_m2 = fe.gauge_latitude(pair)
if abs(gauge_m2 - phi_r2) > 1.0e-12:
    raise AssertionError(f"phi_r from the anchors is {math.degrees(gauge_m2)} degrees and the "
                         f"mesh was built for {math.degrees(phi_r2)}")
state_m2 = loop_to(target10, run_anchors=pair, gauge_latitude_rad=phi_r2, mesh=mesh_m2)
E2, arrivals2, Phi_target2, C_target2 = estimate_from(state_m2, pair, phi_r2)
order2, labels2, produced2 = produce_at(state_m2, target10, Phi_target2, C_target2,
                                        E2.label_pressure_Pa)

# The M = 1 product on the same levels, for the identity. The union at M = 2 is the two anchors'
# arrival levels, so the M = 1 profile is interpolated in Phi onto the M = 2 levels, in ln N and
# in ln p, which is the one interpolation decision C names.
m1_Phi = crossing(state10, target10, curves10.geopotential_m2s2)
m1_lnN = crossing(state10, target10, lnN10)
m1_order = np.argsort(m1_Phi)
same_lnN = np.interp(Phi_target2[order2], m1_Phi[m1_order], m1_lnN[m1_order])
m1_p = produced10.pressure_Pa
m1_T = produced10.temperature_K
m1_Phi_top_down = m1_Phi[order10]
same_p = np.exp(np.interp(Phi_target2[order2], m1_Phi_top_down[::-1],
                          np.log(m1_p)[::-1]))
same_T = np.interp(Phi_target2[order2], m1_Phi_top_down[::-1], m1_T[::-1])
inside = ((Phi_target2[order2] >= m1_Phi.min()) & (Phi_target2[order2] <= m1_Phi.max()))
d_lnN_target = float(np.max(np.abs(C_target2[order2][inside] - same_lnN[inside])))
d_p_target = float(np.max(np.abs(produced2.pressure_Pa[inside] / same_p[inside] - 1.0)))
d_T_target = float(np.max(np.abs(produced2.temperature_K[inside] / same_T[inside] - 1.0)))

# The variant of ruling 1, run beside the identity and reported, not bounded: the same estimate
# with the second anchor written from the separate M = 1 run to 60 N, whose lattice starts one
# spacing below the Lindal latitude. The state is this run's, so no outer loop is repeated: the
# anchor's own levels and `ln N` are the only thing that changes, and what it measures is the
# tracing's discretization between two lattices.
def estimate_with(arrival_second):
    """The estimate and the target profile with a substituted second anchor, on this state."""
    E = fe.estimate([arrivals2[0], arrival_second], SIGMA_K, phi_r2)
    to_target = tr.trace(state_m2.mesh, state_m2.shear_integral, E.geopotential_m2s2, phi_r2,
                         target10)
    latitude_c, slope = lk.composition_term(inputs.composition, E.label_pressure_Pa)
    at_target = int(np.flatnonzero(to_target.latitude_rad == target10)[0])
    carried = tr.transfer(state_m2.isobar_kernel, to_target, E.C, E.label_pressure_Pa, latitude_c,
                          slope)[:, at_target]
    Phi_here = to_target.geopotential_m2s2[:, at_target]
    order_here = np.argsort(Phi_here)[::-1]
    span = ((Phi_here[order_here] >= m1_Phi.min()) & (Phi_here[order_here] <= m1_Phi.max()))
    against = np.interp(Phi_here[order_here], m1_Phi[m1_order], m1_lnN[m1_order])
    difference = np.asarray(list(E.differences.values())[0], dtype="float64")
    return (float(np.max(np.abs(difference[np.isfinite(difference)]))),
            float(np.max(np.abs(carried[order_here][span] - against[span]))))


lattice_anchors = m2_anchors(1.0, path=lattice_path, slug='synthetic60lattice')
lattice_placed = tr.place(lattice_anchors[1], inputs, closure_field, p_b)
lattice_curves = tr.trace(state_m2.mesh, state_m2.shear_integral,
                          lattice_placed.geopotential_m2s2, target60, phi_r2)
lattice_at_gauge = int(np.flatnonzero(lattice_curves.latitude_rad == phi_r2)[0])
lattice_latitude_c, lattice_slope = lk.composition_term(inputs.composition,
                                                        lattice_placed.label_pressure_Pa)
lattice_lnN_at_gauge = tr.transfer(state_m2.isobar_kernel, lattice_curves,
                                   np.asarray(lattice_anchors[1].ln_N, dtype='float64'),
                                   lattice_placed.label_pressure_Pa,
                                   lattice_latitude_c, lattice_slope)
worst_D_lattice, d_lnN_lattice = estimate_with(fe.Arrival(
    slug='synthetic60lattice', latitude_rad=target60, weight=1.0,
    geopotential_m2s2=lattice_curves.geopotential_m2s2[:, lattice_at_gauge],
    C_i=lattice_lnN_at_gauge[:, lattice_at_gauge],
    sigma_ln_N=fe.sigma_ln_N(lattice_anchors[1]),
    label_pressure_Pa=lattice_placed.label_pressure_Pa))

# The tracing's own identity, bounded separately at 1e-12 (v0.16 ruling 1). No anchor file and no
# placement: the Lindal anchor's own curves are carried to 60 N on this mesh and the traced `Phi`
# and `ln N` there are taken straight back to the gauge, which is the transfer out and back with
# nothing in between. What separates this from check 8's bound is the anchor's placement, and
# bounding the two apart is what keeps either from hiding the other.
curves_lindal_m2, lnN_lindal_m2 = along(state_m2, 0, pair[0])
at_60 = state_m2.mesh.latitude_index(target60)
traced_Phi_60 = curves_lindal_m2.geopotential_m2s2[:, at_60]
traced_lnN_60 = lnN_lindal_m2[:, at_60]
traced_back = tr.trace(state_m2.mesh, state_m2.shear_integral, traced_Phi_60, target60, phi_r2)
traced_back_at_gauge = int(np.flatnonzero(traced_back.latitude_rad == phi_r2)[0])
traced_latitude_c, traced_slope = lk.composition_term(
    inputs.composition, state_m2.placements[0].label_pressure_Pa)
traced_lnN_back = tr.transfer(state_m2.isobar_kernel, traced_back, traced_lnN_60,
                              state_m2.placements[0].label_pressure_Pa,
                              traced_latitude_c, traced_slope)
tracing_identity = float(np.max(np.abs(
    traced_lnN_back[:, traced_back_at_gauge]
    - crossing(state_m2, phi_r2, lnN_lindal_m2))))
tracing_identity_Phi = float(np.max(np.abs(
    traced_back.geopotential_m2s2[:, traced_back_at_gauge]
    - crossing(state_m2, phi_r2, curves_lindal_m2.geopotential_m2s2))))

# Beyond the acceptance list: the map takes the nearest anchor's knots, so it has a seam at the
# midpoint of the two anchors' latitudes, where the knots change hands. The step there is the two
# anchors' disagreement about where their isobars sit, which is what `D_12` measures in `ln N`.
seam_index = int(np.argmin(np.abs(state_m2.mesh.latitude_rad
                                  - 0.5 * (phi_a + target60))))
seam_maps = []
for which in (0, 1):
    seam_Phi = state_m2.curves[which].geopotential_m2s2[:, seam_index]
    seam_labels = np.log(np.asarray(state_m2.placements[which].label_pressure_Pa,
                                    dtype="float64"))
    seam_rising = np.argsort(seam_Phi)
    seam_maps.append(tr.IsobarMap((seam_Phi[seam_rising],), (seam_labels[seam_rising],))
                     .ln_p_at(0, state_m2.mesh.geopotential_m2s2))
seam_step = float(np.max(np.abs(seam_maps[0] - seam_maps[1])))

D12 = np.asarray(list(E2.differences.values())[0], dtype="float64")
both = np.isfinite(D12)
worst_D12 = float(np.max(np.abs(D12[both])))
chi2 = E2.chi_square_reduced[np.isfinite(E2.chi_square_reduced)]
worst_chi = float(np.max(chi2)) if chi2.size else float("nan")
surface = geoid.through_anchor(
    state_m2.mesh.latitude_rad, phi_a, float(np.asarray(lindal.radius_m)[gauge_level]),
    lambda x: closure_field.wind_at(np.asarray(x, dtype="float64"), namelist.gauge_isobar_Pa),
    *CONSTANTS)
synthetic_radius = float(np.asarray(pair[1].radius_m)[pair[1].gauge_level_index])
surface_residual = abs(float(surface.residual(target60, synthetic_radius)))
midpoint = 0.5 * (lindal.latitude_planetocentric_deg + 60.0)

record(8, "M = 2 with the synthetic 60 N anchor traced on this run's own mesh: phi_r at the "
          "weighted centroid, D_12 and the target profile below 1e-7 through a written anchor "
          "and the tracing's own identity below 1e-12, the reduced chi-square below one, and the "
          "synthetic anchor on the reference surface to 1 m (SPEC_04 v0.16 section 16 ruling 1)",
       (abs(math.degrees(phi_r2) - midpoint) <= 1e-9 and worst_D12 <= 1e-7
        and (not np.isfinite(worst_chi) or worst_chi <= 1.0)
        and max(d_lnN_target, d_p_target, d_T_target) <= 1e-7 and surface_residual <= 1.0
        and tracing_identity <= 1e-12),
       f"spacings {spacings()}, mesh {state_m2.mesh.shape[0]} by {state_m2.mesh.shape[1]}, "
       f"outer loop {state_m2.record['passes']} passes\n"
       f"RELAXATION: casspian.lib.control._refuse_dirty_commit is replaced for this script only, "
       f"so that the synthetic anchors this script writes on a working tree that is not clean can "
       f"be loaded; their commit is {synthetic_commit!r}\n"
       f"the mesh was built first, for this run's required latitudes, and every M = 2 run is "
       f"formed on it (ruling 1): 10 N the target, {math.degrees(phi_a):.6f} the Lindal anchor, "
       f"{math.degrees(phi_r2):.6f} the gauge, {math.degrees(phi_r_season):.6f} the gauge the "
       f"season variant of check 10 moves to, and 60 N the second anchor. Both gauge latitudes "
       f"are nodes because three runs that are compared level for level have to be on one mesh; "
       f"the second is computable before anything is traced, the season column being the script's "
       f"own and the measurement column shared\n"
       f"the M = 2 namelist's four inputs are byte copies of the run's, so the anchors it loads "
       f"are placed and traced with the run's own loaded inputs\n"
       f"the synthetic anchor (decision M): written from the Lindal anchor's trace to 60 N on "
       f"this mesh, radius_m and height_above_anchor_isobar_m from the 60 N column at the arrival "
       f"levels, refractivity the transferred N, the thermo group's pressures the Lindal anchor's "
       f"tabulated ones and its temperatures p R_bar / (k_B N), anchor radius "
       f"{anchor_radius60:.1f} m and psi "
       f"{psi60:.4f} degrees at 60 N, slug 'synthetic60'\n"
       f"    its refractivity_uncertainty is the Lindal anchor's relative uncertainty, so that "
       f"sigma_ln_N is identical on the two anchors and phi_r is the midpoint the specification "
       f"states; copying the absolute column instead would move phi_r off it\n"
       f"    the derived companions it carries (number density, mean refractivity, mean molar "
       f"mass) are the Lindal anchor's and are read by nothing in the transfer\n"
       f"phi_r = {math.degrees(phi_r2):.6f} degrees, the midpoint {midpoint:.6f} of "
       f"{lindal.latitude_planetocentric_deg:.6f} and 60, to "
       f"{abs(math.degrees(phi_r2) - midpoint):.1e} degrees, and what the anchors as they arrive "
       f"give through gauge_latitude to {abs(gauge_m2 - phi_r2):.1e} radians\n"
       f"the union at phi_r has {E2.geopotential_m2s2.size} levels from the two anchors' "
       f"{p_tab.size} each; largest |D_12| over the {int(np.sum(both))} common levels "
       f"{worst_D12:.3e}, bound 1e-7; reduced chi-square largest {worst_chi:.3e}, bound 1\n"
       f"the target profile against the M = 1 product on the same levels: ln N "
       f"{d_lnN_target:.3e}, p {d_p_target:.3e} relative, T {d_T_target:.3e} relative, bound "
       f"1e-7 on each over the {int(np.sum(inside))} levels the M = 1 profile spans\n"
       f"the tracing's own identity, bounded apart from the anchor's placement (ruling 1): the "
       f"Lindal anchor's curves carried to 60 N on this mesh and taken straight back to the "
       f"gauge, with no anchor file and no placement in between, return its own C_i to "
       f"{tracing_identity:.3e}, bound 1e-12, and its own Phi_k to {tracing_identity_Phi:.3e} "
       f"m2/s2. The difference between that and the {worst_D12:.3e} above is decision G's "
       f"placement of an anchor read from a file, which is what the floor is\n"
       f"the synthetic anchor's radius at the gauge isobar against the reference surface marched "
       f"from the Lindal anchor: {surface_residual:.3e} m, bound 1 m\n"
       f"the isobar labels at phi_r agree between the anchors to "
       f"{float(np.max(E2.label_spread_ln_p)):.3e} in ln p\n"
       f"the variant of ruling 1, reported and not bounded: the same anchor written from the "
       f"separate M = 1 run to 60 N of check 2, whose mesh is uniform from one spacing below "
       f"{lindal.latitude_planetocentric_deg:.6f} degrees instead of below 10, differs from this "
       f"one by {lattice_lnN:.3e} in ln N and {lattice_Phi:.1f} m2/s2 in Phi_k at 60 N, and "
       f"substituted into this run's estimate it gives |D_12| {worst_D_lattice:.3e} and a target "
       f"difference of {d_lnN_lattice:.3e}. That is the tracing's discretization between two "
       f"lattices at the namelist's spacings, which is what the construction as posed measured "
       f"and called an identity\n"
       f"beyond the acceptance list, the map's seam: the knots change hands at the midpoint of "
       f"the two anchors' latitudes, and the two anchors' maps differ there by "
       f"{seam_step:.3e} in ln p. Under the closure wind this reaches nothing, since "
       f"du/dln p is exactly zero and the map's slopes multiply it")
# ---------------------------------------------------------------------------------------------
# 9. The validation anchor
# ---------------------------------------------------------------------------------------------
validation = m2_anchors(0.0)
phi_r_v = fe.gauge_latitude(validation)
state_v = loop_to(target10, run_anchors=validation, gauge_latitude_rad=phi_r_v,
                  mesh=mesh_m2)
Ev, arrivals_v, Phi_target_v, C_target_v = estimate_from(state_v, validation, phi_r_v)
own = np.argsort(arrivals_v[0].geopotential_m2s2)
identity_v = float(np.max(np.abs(Ev.C - arrivals_v[0].C_i[own])))
Dv = np.asarray(list(Ev.differences.values())[0], dtype="float64")
finite_v = np.isfinite(Dv)
record(9, "a validation anchor (weight 0) leaves phi_r and C at the M = 1 values and reports its "
          "own C_i and D_12",
       (abs(phi_r_v - phi_a) == 0.0 and identity_v <= 1e-14
        and Ev.geopotential_m2s2.size == Ev.C.size and np.all(np.isnan(Ev.chi_square_reduced))),
       f"spacings {spacings()}, mesh {state_v.mesh.shape[0]} by {state_v.mesh.shape[1]}, "
       f"outer loop {state_v.record['passes']} passes; the namelist is the M = 2 one with the "
       f"synthetic anchor's weight set to 0\n"
       f"phi_r = {math.degrees(phi_r_v):.6f} degrees, the construction anchor's own, against "
       f"{math.degrees(phi_r2):.6f} at M = 2; the difference from the M = 1 value is "
       f"{abs(phi_r_v - phi_a):.1e} radians\n"
       f"C on the union of {Ev.geopotential_m2s2.size} levels equals the construction anchor's "
       f"own C_i to {identity_v:.1e}: the validation anchor enters neither the weights nor the "
       f"mean, and the reduced chi-square is absent at every level\n"
       f"the validation anchor is still placed, traced and reported: its C_i runs "
       f"{float(np.nanmin(Ev.C_i[1])):+.6f} to {float(np.nanmax(Ev.C_i[1])):+.6f} and D_12 over "
       f"the {int(np.sum(finite_v))} common levels runs {float(np.nanmin(Dv[finite_v])):+.3e} to "
       f"{float(np.nanmax(Dv[finite_v])):+.3e}")

# ---------------------------------------------------------------------------------------------
# 10. A nonzero season column, set after the hook
# ---------------------------------------------------------------------------------------------
seasoned = (pair[0], dataclasses.replace(
    pair[1], sigma_ln_N_season=np.full(pair[1].ln_N.shape, SEASON_TERM),
    season_term="set by the SPEC_04 Step 4 acceptance, after the hook"))
phi_r_s = fe.gauge_latitude(seasoned)
if abs(phi_r_s - phi_r_season) > 1.0e-12:
    raise AssertionError(f'the season gauge is {phi_r_s} rad and the mesh was built for '
                         f'{phi_r_season}')
state_s = loop_to(target10, run_anchors=seasoned, gauge_latitude_rad=phi_r_s, mesh=mesh_m2)
Es, arrivals_s, Phi_target_s, C_target_s = estimate_from(state_s, seasoned, phi_r_s)
order_s = np.argsort(Phi_target_s)[::-1]
_, labels_s, produced_s = produce_at(state_s, target10, Phi_target_s, C_target_s,
                                     Es.label_pressure_Pa)
inside_s = ((Phi_target_s[order_s] >= m1_Phi.min()) & (Phi_target_s[order_s] <= m1_Phi.max()))
same_lnN_s = np.interp(Phi_target_s[order_s], m1_Phi[m1_order], m1_lnN[m1_order])
same_p_s = np.exp(np.interp(Phi_target_s[order_s], m1_Phi_top_down[::-1], np.log(m1_p)[::-1]))
same_T_s = np.interp(Phi_target_s[order_s], m1_Phi_top_down[::-1], m1_T[::-1])
d_lnN_s = float(np.max(np.abs(C_target_s[order_s][inside_s] - same_lnN_s[inside_s])))
d_p_s = float(np.max(np.abs(produced_s.pressure_Pa[inside_s] / same_p_s[inside_s] - 1.0)))
d_T_s = float(np.max(np.abs(produced_s.temperature_K[inside_s] / same_T_s[inside_s] - 1.0)))


def mean_weight(E, index):
    """The synthetic anchor's mean weight where it is present. The two runs place the gauge at
    different latitudes, so their unions differ and the comparison is of means, not of levels."""
    w = E.estimate_weight[index][E.present[index]]
    return float(np.mean(w))


weight_ratio = mean_weight(Es, 1) / mean_weight(E2, 1)
record(10, "a nonzero season column set on the synthetic anchor after the hook moves phi_r and "
           "the weights and changes nothing else in the identity test",
       (abs(phi_r_s - phi_r2) > 1e-6 and weight_ratio < 1.0
        and max(d_lnN_s, d_p_s, d_T_s) <= 1e-7 and pair[1].sigma_ln_N_season.max() == 0.0),
       f"spacings {spacings()}, mesh {state_s.mesh.shape[0]} by {state_s.mesh.shape[1]}, the M = 2 run's own (ruling 1)\n"
       f"the hook leaves the column zero and says so: the anchor as it arrives carries "
       f"sigma_ln_N_season max {float(pair[1].sigma_ln_N_season.max())!r} and season_term "
       f"{pair[1].season_term!r}, and its propagation record is "
       f"{pair[1].propagation['propagation']!r}\n"
       f"the script then sets {SEASON_TERM:g} on every level of the synthetic anchor, which is "
       f"{SEASON_TERM / float(np.mean(fe.sigma_ln_N(pair[1]))):.2f} times its mean measurement "
       f"term\n"
       f"phi_r moves from {math.degrees(phi_r2):.6f} to {math.degrees(phi_r_s):.6f} degrees, "
       f"toward the anchor with the smaller uncertainty; the synthetic anchor's mean weight in "
       f"the estimate falls to {weight_ratio:.4f} of its M = 2 value\n"
       f"the target profile is still the M = 1 product's: ln N to {d_lnN_s:.3e}, p to "
       f"{d_p_s:.3e} relative and T to {d_T_s:.3e} relative, bound 1e-7 on each over the "
       f"{int(np.sum(inside_s))} levels the M = 1 profile spans. The weights move where the "
       f"estimate is placed, not what the two anchors say there")

loop_report.append(pass_line("closure wind, M = 2 to 10 N", state_m2))
loop_report.append(pass_line("closure wind, M = 2 with a validation anchor", state_v))
loop_report.append(pass_line("closure wind, M = 2 with a season column", state_s))
loop_states.extend([state_m2, state_v, state_s])

# ---------------------------------------------------------------------------------------------
# 11. A sheared synthetic state: out and back, ordering, convergence
# ---------------------------------------------------------------------------------------------
if 11 in CARRIED:
    carry(11)
else:
    BETA = 0.1
    sheared = inputs.wind.copy(deep=True)
    u_reference_grid = np.broadcast_to(u_reference_file[:, None], (lat_file.size, p_file.size))
    factor = np.where(p_file < p_ref, 1.0 + BETA * np.log(p_ref / p_file), 1.0)
    u_sheared = u_reference_grid * factor[None, :]
    sheared["u_total_ms"] = (("latitude_planetocentric", "pressure"), u_sheared,
                             dict(inputs.wind["u_total_ms"].attrs))
    sheared["u_shear_ms"] = (("latitude_planetocentric", "pressure"),
                             u_sheared - u_reference_grid,
                             dict(inputs.wind["u_shear_ms"].attrs))
    sheared.attrs["title"] = "Sheared synthetic wind for the SPEC_04 Step 4 acceptance"
    sheared_path = cio.write(HERE / "step04_4_sheared_wind.nc", sheared, "wind",
                             created_by="accept_step04_4")
    sheared_back = read_closed(sheared_path, "wind")
    sheared_field = wf.WindField(sheared_back)


    def out_and_back(scale):
        state = loop_to(target10, field=sheared_field, scale=scale)
        curves, ln_N = along(state)
        Phi_out = crossing(state, target10, curves.geopotential_m2s2)
        lnN_out = crossing(state, target10, ln_N)
        back = tr.trace(state.mesh, state.shear_integral, Phi_out, target10, phi_a)
        latitude_c, slope = lk.composition_term(inputs.composition,
                                                state.placements[0].label_pressure_Pa)
        lnN_back = tr.transfer(state.isobar_kernel, back, lnN_out,
                               state.placements[0].label_pressure_Pa, latitude_c, slope)
        home = int(np.flatnonzero(back.latitude_rad == phi_a)[0])
        start_Phi = state.placements[0].geopotential_m2s2
        d_Phi = float(np.max(np.abs(back.geopotential_m2s2[:, home] - start_Phi)))
        d_lnN = float(np.max(np.abs(lnN_back[:, home] - lindal.ln_N)))
        return state, d_Phi, d_lnN, stays_ordered(curves) and stays_ordered(back)


    # The pair of spacings is the namelist's own and twice it, rather than the namelist's and half:
    # this check has no absolute tolerance to meet, only the fall under halving, and the finer member
    # of the pair is the run's own spacing either way. Check 3, which does have absolute bounds at
    # both, runs at the namelist's spacings and at half of them.
    SHEARED_COARSE = 2.0

    def ratio(coarse, fine):
        """`coarse / fine` for the report line; a fine value of exactly zero is named, not divided by."""
        if fine != 0.0:
            return f"{coarse / fine:.2f}"
        return "none (both zero)" if coarse == 0.0 else "inf (the finer is zero)"

    state_shc, dPhi_shc, dlnN_shc, ordered_shc = out_and_back(SHEARED_COARSE)
    state_sh, dPhi_sh, dlnN_sh, ordered_sh = out_and_back(1.0)
    loop_report.append(pass_line("sheared wind to 10 N, doubled", state_shc))
    loop_report.append(pass_line("sheared wind to 10 N", state_sh))
    loop_states.extend([state_shc, state_sh])
    sum_identity = float(np.max(np.abs(
        np.asarray(sheared_back["u_total_ms"].values)
        - np.asarray(sheared_back["u_reference_ms"].values)[:, None]
        - np.asarray(sheared_back["u_shear_ms"].values))))
    record(11, "a sheared synthetic wind file returns Phi_k and ln N_k on the way out and back, to a "
               "tolerance that falls under halving both spacings, with the curves ordered and the "
               "outer loop converging in under ten passes",
           (dPhi_sh < dPhi_shc and dlnN_sh < dlnN_shc and ordered_sh and ordered_shc
            and state_sh.record["passes"] < 10 and state_shc.record["passes"] < 10),
           f"the file: u_total(phi, p) = u_reference(phi) [1 + {BETA} ln(p_ref / p)] above "
           f"p_ref = {p_ref:g} Pa and u_reference below, written to "
           f"{sheared_path.name} in the three parts and read back through the kind W reader, which "
           f"checks the sum identity ({sum_identity:.1e} m/s) and the poles\n"
           f"    on the anchor's column it runs {float(sheared_field.wind_at(phi_a, p_tab[TOP])):.3f} "
           f"m/s at the top level against "
           f"{float(closure_field.wind_at(phi_a, p_tab[TOP])):.3f} under the closure wind\n"
           f"out and back, phi_c to 10 N and back, at {spacings(SHEARED_COARSE)}: largest "
           f"|Phi_k - Phi_k(out and back)| {dPhi_shc:.3e} m2/s2, largest |ln N_k - ln N_k(out and "
           f"back)| {dlnN_shc:.3e}\n"
           f"with both spacings halved, at {spacings()}, the namelist's own: {dPhi_sh:.3e} m2/s2 and "
           f"{dlnN_sh:.3e}, ratios {ratio(dPhi_shc, dPhi_sh)} and {ratio(dlnN_shc, dlnN_sh)}; no order "
           f"is claimed, since u' is piecewise constant and both the trapezoid and the scheme are "
           f"first order across the file's nodes\n"
           f"the pair is the namelist's spacings and twice them rather than the namelist's and half: "
           f"this check has no absolute tolerance to meet, only the fall under halving, and the finer "
           f"member is the run's own spacing either way. Check 3, which does have bounds at both, "
           f"runs at the namelist's spacings and at half of them\n"
           f"the curves stay ordered at every node on both meshes: {ordered_shc} and {ordered_sh}\n"
           f"the outer loop took {state_shc.record['passes']} passes at the coarser spacings and "
           f"{state_sh.record['passes']} at the namelist's, bound 10")

# ---------------------------------------------------------------------------------------------
# 7 (reported last, since it gathers every run)
# ---------------------------------------------------------------------------------------------
record(7, "the outer loop's pass count and residuals are reported for every wind file",
       all(s.record["residual_ln_p"] < s.record["relative_tolerance_ln_p"]
           for s in loop_states),
       "\n".join(loop_report)
       + ("" if not CARRIED else
          f"\n    the runs of the carried checks {sorted(CARRIED)} are not in this pass, and "
          f"their pass counts are in those checks' own carried rows")
       + "\na wind whose columns do not vary with pressure has a wind on the mesh that does not "
       "depend on the map, so the second pass reproduces the first exactly and the residual is "
       "zero; the code does not special-case it, and the first pass's residual is the departure "
       "of the traced map from the flat guess")

ctl._refuse_dirty_commit = _original_refuse

# ---------------------------------------------------------------------------------------------
results.sort(key=lambda row: row[0])
write_output()
print(flush=True)
passed = sum(1 for _, _, ok, _ in results if ok)
print(f"{passed} of {len(results)} checks pass", flush=True)
with (HERE / "output.txt").open("a", encoding="utf-8") as handle:
    handle.write(f"\n{passed} of {len(results)} checks pass\n")
sys.exit(0 if passed == len(results) else 1)
