"""Acceptance checks for SPEC_07 v0.3 Step 1: the wind interpolant with a continuous slope.

Decision L2 (PCHIP in latitude on every pressure row, then PCHIP in `ln p` at the point) in
`lib.windfield.WindField`, read by the model and, through `refrac.anchor.wind_of_latitude`, by the
reduction. Every run is made in a copy under `reports/step07_1/`, from the working tree with the
`-dirty` refusal relaxed in this process only, as in SPEC_05 and SPEC_06; the SPEC_05 experiments
are rerun in place by `tests/step05_4/run_experiments.py`, as SPEC_06 did. The SPEC_06 values are
the products set aside before this step under `reports/step07_1/spec06/`.

"Every wind file" is the registered source and run wind files and the eleven experiment run
files. The nine checks of v0.3 §3; `--no-run` reuses the candidates and the experiments already
made. Writes `reports/step07_1/output.txt` and two figures. Run from the repository root.
"""

import json
import shutil
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import PchipInterpolator

from casspian.lib import control as ctl
from casspian.lib import kernel as lk
from casspian.lib.windfield import WindField

HERE = Path("reports/step07_1")
SPEC06 = HERE / "spec06"
TREE = HERE / "chain" / "tree"
RUN = "--no-run" not in sys.argv
NEWLINE = chr(10)
WIND_FILES = [
    "occul_data/lindal/lindal_wind.nc",
    "forward/lindal_closure/inputs/lindal_closure_wind_source.nc",
    "forward/lindal_closure/inputs/lindal_closure_wind.nc",
    "forward/lindal_transfer/inputs/lindal_transfer_wind_source.nc",
    "forward/lindal_transfer/inputs/lindal_transfer_wind.nc",
] + sorted(p.as_posix() for p in Path("forward").glob("shear_*/inputs/shear_*_wind.nc"))
results = []
_started = time.time()


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}", flush=True)
    for line in str(detail).splitlines():
        print(f"        {line}", flush=True)


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


def relaxed():
    """The `-dirty` refusal relaxed for this process, returning the original."""
    original = ctl._refuse_dirty_commit
    ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
        attrs.get("casspian_git_commit", ""))
    return original


def turnover_wind():
    """A synthetic kind W field whose vertical profile turns over at a pressure node at some
    latitudes and not at others: `u = U cos(phi) (x - c(phi))^2` with the vertex `c(phi)` swept
    across a node in `x = ln p` as latitude changes. Read by `WindField` like a file."""
    latitude = np.linspace(-90.0, 90.0, 37)
    pressure = 10.0 ** np.linspace(2.0, 5.0, 16)
    x = np.log(pressure)
    vertex = x[7] + 1.2 * np.sin(np.radians(latitude))
    u = 100.0 * np.cos(np.radians(latitude))[:, None] * (x[None, :] - vertex[:, None]) ** 2 / 4.0
    u[0] = u[-1] = 0.0
    return xr.Dataset(
        {"u_total_ms": (("latitude_planetocentric", "pressure"), u),
         "latitude_planetocentric_deg": (("latitude_planetocentric",), latitude),
         "pressure_Pa": (("pressure",), pressure),
         "reference_level_pressure_Pa": ((), float(pressure[7]))})


def load_field(name):
    if name == "synthetic turnover":
        return WindField(turnover_wind())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(name, engine="netcdf4").load()
    tree.close()
    return WindField(tree.to_dataset(inherit=False))


def inside(nodes, per_interval):
    t = (np.arange(per_interval) + 0.5) / per_interval
    return (nodes[:-1, None] + np.diff(nodes)[:, None] * t[None, :]).ravel()


FIELDS = {name: load_field(name) for name in WIND_FILES}

# ---------------------------------------------------------------------------------------------
# 1. The rule: every row and every column against scipy's PchipInterpolator
# ---------------------------------------------------------------------------------------------
lines, ok = [], True
for name, W in FIELDS.items():
    lat, x, u = W.latitude_rad, W.ln_pressure, W.u_total_ms
    phis, xs = inside(lat, 5), inside(x, 5)
    worst = {"row value": 0.0, "row du/dphi": 0.0, "column value": 0.0, "column du/dln p": 0.0}
    rows = PchipInterpolator(lat, u, axis=0)
    row_values, row_slopes = rows(phis), rows.derivative()(phis)
    for k in range(x.size):
        got = W.wind_at(phis, np.full(phis.shape, np.exp(x[k])))
        d_phi, _ = W.wind_derivatives(phis, np.full(phis.shape, np.exp(x[k])))
        worst["row value"] = max(worst["row value"], float(np.max(np.abs(got - row_values[:, k]))))
        worst["row du/dphi"] = max(worst["row du/dphi"],
                                   float(np.max(np.abs(d_phi - row_slopes[:, k]))))
    for j in range(lat.size):
        column = PchipInterpolator(x, u[j])
        got = W.wind_at(np.full(xs.shape, lat[j]), np.exp(xs))
        _, d_x = W.wind_derivatives(np.full(xs.shape, lat[j]), np.exp(xs))
        worst["column value"] = max(worst["column value"], float(np.max(np.abs(got - column(xs)))))
        worst["column du/dln p"] = max(worst["column du/dln p"],
                                       float(np.max(np.abs(d_x - column.derivative()(xs)))))
    scale_u = float(np.max(np.abs(u))) or 1.0
    scale_phi = float(np.max(np.abs(row_slopes))) or 1.0
    scale_x = float(np.max(np.abs(np.diff(u, axis=1) / np.diff(x)))) or 1.0
    rel = {"row value": worst["row value"] / scale_u, "column value": worst["column value"] / scale_u,
           "row du/dphi": worst["row du/dphi"] / scale_phi, "column du/dln p": worst["column du/dln p"] / scale_x}
    ok = ok and max(rel.values()) <= 1e-14
    lines.append(f"{name}: " + ", ".join(f"{k} {v:.1e}" for k, v in rel.items()))
record(1, "on every row and column of every wind file, values and first derivatives agree with "
          "scipy's PchipInterpolator to 1e-14 relative (of the file's largest value of each)",
       ok, NEWLINE.join(lines))

# ---------------------------------------------------------------------------------------------
# 2. Nodes exact
# ---------------------------------------------------------------------------------------------
lines, ok = [], True
for name, W in FIELDS.items():
    L, X = np.meshgrid(W.latitude_rad, W.ln_pressure, indexing="ij")
    nodes = bool(np.array_equal(W.wind_at(L, np.exp(X)), W.u_total_ms))
    column = int(np.flatnonzero(W.pressure_Pa == W.reference_pressure_Pa)[0])
    reference = bool(np.array_equal(W.reference_wind(W.latitude_rad), W.u_total_ms[:, column]))
    ok = ok and nodes and reference
    lines.append(f"{name}: every node {nodes}, reference_wind at every latitude node {reference}")
record(2, "wind_at at every node of every wind file returns u_total bit for bit, and reference_wind "
          "at every latitude node the reference row", ok, NEWLINE.join(lines))

# ---------------------------------------------------------------------------------------------
# 3. No overshoot, and the sample for checks 4 and 5
# ---------------------------------------------------------------------------------------------
lines, ok = [], True
zero_vertical = {}
for name, W in FIELDS.items():
    lat, x, u = W.latitude_rad, W.ln_pressure, W.u_total_ms
    t = (np.arange(20) + 0.5) / 20
    xs = (x[:-1, None] + np.diff(x)[:, None] * t[None, :]).ravel()
    cell_k = np.repeat(np.arange(x.size - 1), 20)
    worst, tolerance = 0.0, 1e-12 * (float(np.max(np.abs(u))) or 1.0)
    largest_dx = 0.0
    for j in range(lat.size - 1):
        phis = lat[j] + (lat[j + 1] - lat[j]) * t
        P, XS = np.meshgrid(phis, xs, indexing="ij")
        got = W.wind_at(P, np.exp(XS))
        corners = np.stack([u[j, cell_k], u[j + 1, cell_k], u[j, cell_k + 1], u[j + 1, cell_k + 1]])
        low, high = corners.min(axis=0)[None, :], corners.max(axis=0)[None, :]
        worst = max(worst, float(np.max(np.maximum(low - got, got - high))))
        if name == "forward/lindal_closure/inputs/lindal_closure_wind.nc":
            largest_dx = max(largest_dx, float(np.max(np.abs(W.wind_derivatives(P, np.exp(XS))[1]))))
    zero_vertical[name] = largest_dx
    ok = ok and worst <= tolerance
    lines.append(f"{name}: largest excursion beyond a cell's corner values {max(worst, 0.0):.2e} m/s "
                 f"(tolerance {tolerance:.1e})")
record(3, "on 20 points per cell along each axis of every wind file, every value within its cell's "
          "four corner values, to 1e-12 of the file's largest |u|", ok, NEWLINE.join(lines))

# ---------------------------------------------------------------------------------------------
# 4. Continuous slopes, the loci, and the interpolant's own derivatives
# ---------------------------------------------------------------------------------------------
FIELDS_4 = dict(FIELDS, **{"synthetic turnover": load_field("synthetic turnover")})


def loci(W):
    """Latitudes off the poles where two adjacent vertical secants of the latitude interpolated
    rows change from agreeing in sign to not, at an interior pressure node; and the jump of
    `(du/dphi)` across each, measured on the interval below and above that node."""
    lat, x = W.latitude_rad, W.ln_pressure
    fine = np.linspace(lat[0], lat[-1], 40 * (lat.size - 1) + 1)
    fine = fine[np.abs(fine) < np.radians(89.0)]
    rows = W._rows(fine)
    secants = np.diff(rows, axis=1)
    agree = secants[:, :-1] * secants[:, 1:] > 0
    found = []
    for k in range(agree.shape[1]):
        for j in np.flatnonzero(agree[1:, k] != agree[:-1, k]):
            lo, hi = fine[j], fine[j + 1]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                s = np.diff(W._rows(np.array([lo, mid]))[:, k:k + 3], axis=1)
                (lo, hi) = (mid, hi) if (s[0, 0] * s[0, 1] > 0) == (s[1, 0] * s[1, 1] > 0) else (lo, mid)
            at = 0.5 * (lo + hi)
            xs = np.concatenate([x[k] + (x[k + 1] - x[k]) * np.linspace(0.05, 0.95, 7),
                                 x[k + 1] + (x[k + 2] - x[k + 1]) * np.linspace(0.05, 0.95, 7)])
            eps = 1e-9
            left = W.wind_derivatives(np.full(xs.shape, at - eps), np.exp(xs))[0]
            right = W.wind_derivatives(np.full(xs.shape, at + eps), np.exp(xs))[0]
            rows_d = W._rows_dphi(np.array([at]))[0]
            found.append((float(np.degrees(at)), float(np.exp(x[k + 1])),
                          float(np.max(np.abs(right - left))),
                          float(np.max(np.abs(np.diff(rows_d[k:k + 3]))))))
    return found


lines, ok = [], True
for name, W in FIELDS_4.items():
    lat, x, u = W.latitude_rad, W.ln_pressure, W.u_total_ms
    phis, xs = inside(lat, 3), inside(x, 3)
    # The one sided limits are read 1e-12 from the node: PCHIP's second derivative jumps at a node,
    # so a difference taken farther out measures that jump times the distance (at 1e-10 it is
    # 3.8e-8 of the largest du/dphi, and falls tenfold with each tenfold smaller step).
    eps = 1e-12
    # One-sided limits across every pressure node, at the latitude samples, and across every
    # latitude node, at the pressure samples.
    P, XN = np.meshgrid(phis, x[1:-1], indexing="ij")
    below = W.wind_derivatives(P, np.exp(XN - eps))
    above = W.wind_derivatives(P, np.exp(XN + eps))
    LN, XS = np.meshgrid(lat[1:-1], xs, indexing="ij")
    south = W.wind_derivatives(LN - eps, np.exp(XS))
    north = W.wind_derivatives(LN + eps, np.exp(XS))
    scale_phi = max(float(np.max(np.abs(below[0]))), float(np.max(np.abs(south[0])))) or 1.0
    scale_x = max(float(np.max(np.abs(below[1]))), float(np.max(np.abs(south[1])))) or 1.0
    jump = {"du/dln p across pressure nodes": float(np.max(np.abs(above[1] - below[1]))) / scale_x,
            "du/dphi across pressure nodes": float(np.max(np.abs(above[0] - below[0]))) / scale_phi,
            "du/dln p across latitude nodes": float(np.max(np.abs(north[1] - south[1]))) / scale_x,
            "du/dphi across latitude nodes": float(np.max(np.abs(north[0] - south[0]))) / scale_phi}
    # Decision L's largest jump of du/dphi across a latitude node, on the node rows.
    slopes = np.diff(u, axis=0) / np.diff(lat)[:, None]
    linear_jump = float(np.max(np.abs(np.diff(slopes, axis=0))))
    found = loci(W)
    largest_locus = max((f[2] for f in found), default=0.0)
    # Centered differences away from the loci.
    rng = np.random.default_rng(7)
    sp = rng.uniform(lat[1], lat[-2], 400)
    sx = rng.uniform(x[1], x[-2], 400)
    if found:
        keep = np.min(np.abs(sp[:, None] - np.radians([f[0] for f in found])[None, :]), axis=1) > 1e-3
        sp, sx = sp[keep], sx[keep]
    d_phi, d_x = W.wind_derivatives(sp, np.exp(sx))
    h = 1e-6
    fd_phi = (W.wind_at(sp + h, np.exp(sx)) - W.wind_at(sp - h, np.exp(sx))) / (2 * h)
    fd_x = (W.wind_at(sp, np.exp(sx + h)) - W.wind_at(sp, np.exp(sx - h))) / (2 * h)
    centred = max(float(np.max(np.abs(d_phi - fd_phi))) / scale_phi,
                  float(np.max(np.abs(d_x - fd_x))) / scale_x)
    passed = (max(jump.values()) <= 1e-9 and centred <= 1e-6
              and (largest_locus < linear_jump or not found))
    ok = ok and passed
    text = f"{name}: " + ", ".join(f"{k} {v:.1e}" for k, v in jump.items())
    text += f"; centred differences {centred:.1e}; loci off the poles {len(found)}"
    if found:
        worst = max(found, key=lambda f: f[2])
        text += (f", largest jump of du/dphi {worst[2]:.3e} m/s per rad at {worst[0]:.3f} deg and the "
                 f"node {worst[1]:.4g} Pa, against the difference of the rows' du/dphi across the "
                 f"interval pair {worst[3]:.3e} (ratio {worst[2] / worst[3]:.3f}); decision L's largest "
                 f"jump across a latitude node {linear_jump:.3e}")
    else:
        text += f"; decision L's largest jump of du/dphi across a latitude node {linear_jump:.3e}"
    lines.append(text)
record(4, "across every pressure and latitude node of every wind file the one sided limits of both "
          "derivatives agree to 1e-9 of the file's largest; off the poles the loci's jumps, if any, are "
          "below decision L's latitude node jumps (a synthetic file with a turnover included); both "
          "derivatives agree with centred differences to 1e-6 away from the loci",
       ok, NEWLINE.join(lines))

# ---------------------------------------------------------------------------------------------
# 5. No vertical shear, no vertical derivative
# ---------------------------------------------------------------------------------------------
closure_name = "forward/lindal_closure/inputs/lindal_closure_wind.nc"
record(5, "under the closure wind du/dln p is exactly zero at every sampled point (the SPEC_04 Step 3 "
          "closed form check is step04_3's, in the regression)",
       zero_vertical[closure_name] == 0.0,
       f"largest |du/dln p| over check 3's sample of the closure wind: {zero_vertical[closure_name]!r}")

# ---------------------------------------------------------------------------------------------
# 6. Nearly the same answer, the chain
# ---------------------------------------------------------------------------------------------


def copy_tree():
    """The registered files' sources and control files, in the repository's layout (step05_2)."""
    if TREE.exists():
        shutil.rmtree(TREE)
    shutil.copytree("data_static", TREE / "data_static")
    lindal = TREE / "occul_data" / "lindal"
    (lindal / "raw").mkdir(parents=True)
    for name in ("lindal_build.toml", "lindal_reduction.toml"):
        shutil.copy2(Path("occul_data/lindal") / name, lindal / name)
    for path in Path("occul_data/lindal/raw").iterdir():
        if path.suffix != ".nc":
            shutil.copy2(path, lindal / "raw" / path.name)
    for run in ("lindal_closure", "lindal_transfer"):
        directory = TREE / "forward" / run
        directory.mkdir(parents=True)
        for name in (f"{run}.toml", f"{run}_build.toml"):
            shutil.copy2(Path("forward") / run / name, directory / name)


def build_chain():
    """The order of `tests/step04_0/sweep.py`, on the copy (step05_2's build_candidates)."""
    from casspian.forward import production as forward_production
    from casspian.forward import transfer as forward_transfer
    from casspian.refrac import product as refrac_product
    from casspian.tools.composition import build_composition
    from casspian.tools.gravity import build_gravity, build_rotation
    from casspian.tools.lindal import build_inputs, build_raw
    from casspian.tools.run import run_inputs
    from casspian.tools.wind import build_wind

    lindal = TREE / "occul_data" / "lindal"
    closure, transfer = TREE / "forward" / "lindal_closure", TREE / "forward" / "lindal_transfer"
    build = lindal / "lindal_build.toml"
    build_raw.build(lindal / "raw", lindal / "raw" / "lindal_raw.nc")
    build_gravity.build(build), build_rotation.build(build)
    build_wind.build(build), build_composition.build(build), build_inputs.build(build)
    refrac_product.build_product(lindal / "lindal_reduction.toml")
    run_inputs.build(closure / "lindal_closure_build.toml")
    forward_production.run(closure / "lindal_closure.toml")
    run_inputs.build(transfer / "lindal_transfer_build.toml")
    forward_transfer.run(transfer / "lindal_transfer.toml")


if RUN:
    copy_tree()
    original = relaxed()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            started = time.time()
            build_chain()
            print(f"the chain rebuilt in a copy in {time.time() - started:.0f} s", flush=True)
    finally:
        ctl._refuse_dirty_commit = original

lines, ok = [], True
REDUCTION = "occul_data/lindal/lindal_refractivity.nc"
old_tree, new_tree = tree_of(REDUCTION), tree_of(TREE / REDUCTION)
old_root, new_root = old_tree.to_dataset(inherit=False), new_tree.to_dataset(inherit=False)
d_phi_c = float(new_root["latitude_planetocentric_deg"].values - old_root["latitude_planetocentric_deg"].values)
d_r0 = float(new_root["anchor_isobar_radius_m"].values - old_root["anchor_isobar_radius_m"].values)
record_old, record_new = old_tree["reduction_record"].attrs, new_tree["reduction_record"].attrs
polar = {name: float(record_new[name]) - float(record_old[name])
         for name in ("polar_radius_north_m", "polar_radius_south_m", "polar_radius_mean_m",
                      "polar_asymmetry_m")}
ok = abs(d_phi_c) <= 1e-4 and abs(d_r0) <= 100.0 and all(abs(v) <= 100.0 for v in polar.values())
lines.append(f"the reduction: phi_c {d_phi_c:+.3e} deg, r0 {d_r0:+.3f} m; "
             + ", ".join(f"{k} {v:+.3f} m" for k, v in polar.items()))
for product in ("forward/lindal_closure/output/lindal_closure_profile.nc",
                "forward/lindal_transfer/output/lindal_transfer_profile.nc"):
    a, b = root_of(product), root_of(TREE / product)
    dT = float(np.max(np.abs(values(b, "temperature_K") - values(a, "temperature_K"))))
    dN = float(np.max(np.abs(np.log(values(b, "refractivity") / values(a, "refractivity")))))
    p_name = "pressure_label_Pa" if "pressure_label_Pa" in a else "pressure_Pa"
    dp = float(np.max(np.abs(values(b, p_name) / values(a, p_name) - 1.0)))
    lengths = [n for n in ("radius_m", "altitude_m", "height_above_anchor_isobar_m") if n in a]
    dz = {n: float(np.max(np.abs(values(b, n) - values(a, n)))) for n in lengths}
    passed = dT <= 0.1 and dN <= 1e-3 and dp <= 1e-3 and all(v <= 100.0 for v in dz.values())
    ok = ok and passed
    lines.append(f"{Path(product).name}: T {dT:.3e} K, ln N {dN:.3e}, {p_name} {dp:.3e} relative, "
                 + ", ".join(f"{n} {v:.3f} m" for n, v in dz.items()))
predictor = []
for name in ("occul_data/lindal/lindal_wind.nc",):
    W = FIELDS[name]
    phis = inside(W.latitude_rad, 8)
    column = int(np.flatnonzero(W.pressure_Pa == W.reference_pressure_Pa)[0])
    linear = np.interp(phis, W.latitude_rad, W.u_total_ms[:, column])
    predictor.append(f"{name}: largest |u_L2 - u_L| on the reference level "
                     f"{float(np.max(np.abs(W.reference_wind(phis) - linear))):.3f} m/s; at phi_c "
                     f"{float(W.reference_wind(np.radians(old_root['latitude_planetocentric_deg'].values)) - np.interp(np.radians(float(old_root['latitude_planetocentric_deg'].values)), W.latitude_rad, W.u_total_ms[:, column])):+.4f} m/s")
record(6, "the chain against the registered files before the sweep: phi_c within 1e-4 deg, r0 and the "
          "record's polar radii and asymmetry within 100 m; the closure and transfer products' T within "
          "0.1 K, ln N and pressure within 1e-3 relative, radii and altitudes within 100 m (the at-anchor "
          "identity is step04_5's, in the regression)",
       ok, NEWLINE.join(lines + ["the predictor: " + p for p in predictor]))

# ---------------------------------------------------------------------------------------------
# 7. No steps: the SPEC_05 experiments rerun
# ---------------------------------------------------------------------------------------------
if RUN:
    sys.path.insert(0, str(Path("tests/step05_4").resolve()))
    import run_experiments  # noqa: E402
    run_experiments.main()
now = json.loads(Path("reports/step05_4/results.json").read_text(encoding="utf-8"))


def spec06_product(run, spacing="5e4"):
    return SPEC06 / (f"{run}_profile.nc" if spacing == "5e4" else f"{run}_2p5e4_profile.nc")


def along_isobar(path, level_pressure=1.0e5):
    """`(latitude_deg, S/g, wind latitude nodes)` along the lindal anchor's isobar nearest 1 bar."""
    tree = tree_of(path)
    iso = tree["isobars"].to_dataset(inherit=False)
    est = tree["estimate"].to_dataset(inherit=False)
    root = tree.to_dataset(inherit=False)
    latitude = values(iso, "latitude_planetocentric_deg")
    gauge = int(np.argmin(np.abs(latitude - float(root["gauge_latitude_planetocentric_deg"].values))))
    union = int(np.argmin(np.abs(values(est, "label_pressure_Pa") - level_pressure)))
    k = int(np.argmin(np.abs(values(iso, "geopotential_lindal_m2s2")[:, gauge]
                             - values(est, "geopotential_m2s2")[union])))
    wind = tree["inputs/wind"].to_dataset(inherit=False)
    return latitude, values(iso, "shear_kernel_lindal_per_rad")[k], values(wind, "latitude_planetocentric_deg")


lines, ok = [], True
R7F = "shear_r7f_decay_linp_fine"
lat_old, s_old, wind_nodes = along_isobar(spec06_product(R7F))
lat_new, s_new, _ = along_isobar(now[R7F]["5e4"]["product"])
steps_old, steps_new = np.abs(np.diff(s_old)), np.abs(np.diff(s_new))
at_node = np.array([np.any((wind_nodes > a - 1e-9) & (wind_nodes <= b + 1e-9))
                    for a, b in zip(lat_new[:-1], lat_new[1:])])
pattern = float(np.mean(steps_new[at_node]) / np.mean(steps_new[~at_node]))
first = steps_new.max() * 10.0 <= steps_old.max() and pattern < 2.0
ok = ok and first
lines.append(f"run 7f at 5e4, the isobar nearest 1 bar: largest step of S/g between consecutive curve "
             f"nodes {steps_new.max():.3e} per rad, SPEC_06 {steps_old.max():.3e} (ratio "
             f"{steps_old.max() / steps_new.max():.1f}); mean step at the wind's latitude nodes against "
             f"elsewhere {pattern:.2f}")

from matplotlib.figure import Figure  # noqa: E402

fig = Figure(figsize=(8.0, 5.0))
ax = fig.subplots()
ax.plot(lat_old, s_old, linewidth=0.9, label="SPEC_06 (decision L)")
ax.plot(lat_new, s_new, linewidth=0.9, label="SPEC_07 (decision L2)")
ax.set_xlabel("planetocentric latitude (deg)")
ax.set_ylabel("S/g along the isobar (per rad)")
ax.set_title("Run 7f at 5e4, the isobar nearest 1 bar", fontsize=9)
ax.legend(fontsize=8)
fig.savefig(HERE / "F8_isobar_L_against_L2.png", dpi=130)


def layer_slopes(path):
    d = root_of(path)
    T, p = values(d, "temperature_K"), values(d, "pressure_label_Pa")
    z = values(d, "altitude_m")
    return T, p, z, np.diff(T) / np.diff(np.log(p))


T_old, p_lab, z_old, g_old = layer_slopes(spec06_product(R7F))
T_new, _, z_new, g_new = layer_slopes(now[R7F]["5e4"]["product"])
bar = int(np.argmin(np.abs(p_lab - 1.0e5)))
low = np.flatnonzero((z_new - z_new[bar]) <= 50.0e3)
wind_p = FIELDS[f"forward/{R7F}/inputs/{R7F}_wind.nc"].pressure_Pa
checked, misses = 0, []
for i in range(max(low.min(), 1), p_lab.size - 2):
    a, b = sorted((p_lab[i], p_lab[i + 1]))
    if not np.any((wind_p > a) & (wind_p < b)):
        continue
    neighbours = np.array([g_new[i - 1], g_new[i + 1]])
    widen = 0.1 * float(np.max(np.abs(neighbours)))
    checked += 1
    if not (neighbours.min() - widen <= g_new[i] <= neighbours.max() + widen):
        misses.append(f"layer {i} to {i + 1}: {g_new[i]:.2f} against {neighbours.min():.2f} to "
                      f"{neighbours.max():.2f} K per unit ln p")
second = not misses and checked > 0
ok = ok and second
old_misses = 0
for i in range(max(low.min(), 1), p_lab.size - 2):
    a, b = sorted((p_lab[i], p_lab[i + 1]))
    if np.any((wind_p > a) & (wind_p < b)):
        nb = np.array([g_old[i - 1], g_old[i + 1]])
        w = 0.1 * float(np.max(np.abs(nb)))
        old_misses += not (nb.min() - w <= g_old[i] <= nb.max() + w)
lines.append(f"run 7f's lowest 50 km: {checked} layers straddle a wind pressure node; dT/dln p within its "
             f"neighbours' range widened by 10 percent of their larger magnitude in "
             f"{checked - len(misses)} (SPEC_06: {checked - old_misses})"
             + ("" if not misses else "; outside: " + "; ".join(misses)))

fig = Figure(figsize=(6.0, 6.0))
ax = fig.subplots()
ax.plot(T_old, (z_old - z_old[bar]) / 1e3, linewidth=1.0, label="SPEC_06 (decision L)")
ax.plot(T_new, (z_new - z_new[bar]) / 1e3, linewidth=1.0, label="SPEC_07 (decision L2)")
ax.set_ylim(-20, 120)
ax.set_xlabel("delivered temperature at 10 N (K)")
ax.set_ylabel("height above the 1000 mbar level (km)")
ax.set_title("Run 7f at 5e4, F7's profile, old and new", fontsize=9)
ax.legend(fontsize=8)
fig.savefig(HERE / "F7_profile_L_against_L2.png", dpi=130)

before = json.loads((SPEC06 / "results.json").read_text(encoding="utf-8"))
identity_lines, third = [], True
for run, entry in now.items():
    cells, largest = [], {}
    for spacing in ("5e4", "2.5e4"):
        e, b = entry.get(spacing), before.get(run, {}).get(spacing)
        if e is None:
            continue
        if e["status"] != "completed":
            third = False
            cells.append(f"{spacing} {e['status']}")
            continue
        largest[spacing] = e["pressure_identity_largest_abs"]
        old = b["pressure_identity_largest_abs"]
        unsheared = run in ("shear_r2_uniform", "shear_r3a_nowind", "shear_r3b_nowind_anchor", "shear_r3c_half")
        not_above = largest[spacing] <= old * (1.0 + 1e-3 if unsheared else 1.0)
        third = third and not_above
        cells.append(f"{spacing} {largest[spacing]:.3e} at level {e['pressure_identity_largest_level']} "
                     f"(SPEC_06 {old:.3e}){'' if not_above else ''}{'' if not_above else ' ABOVE'}")
    if len(largest) == 2:
        apart = abs(largest["5e4"] / largest["2.5e4"] - 1.0)
        third = third and apart <= 0.10
        cells.append(f"apart {100 * apart:.1f} percent")
    identity_lines.append(f"    {run}: " + "; ".join(cells))
ok = ok and third
lines.append("the pressure identity, every run, not above SPEC_06 (the unsheared runs to 1e-3 of their "
             "value), the two spacings within 10 percent:")
lines.extend(identity_lines)

# The largest remaining identity, read layer by layer: where its offset p - p_label arises, and the
# layer law of SPEC_06 check 7 (the layer's mass times its jump in ln T relative to run 2) beside it.
worst_run = max(now, key=lambda r: now[r]["5e4"].get("pressure_identity_largest_abs", 0.0))
d = root_of(now[worst_run]["5e4"]["product"])
run2 = root_of(now["shear_r2_uniform"]["5e4"]["product"])
offset = values(d, "pressure_Pa") - values(d, "pressure_label_Pa")
gain = np.diff(offset)
jump = np.diff(np.log(values(d, "temperature_K"))) - np.diff(np.log(values(run2, "temperature_K")))
mass = np.diff(values(run2, "pressure_label_Pa"))
top = np.argsort(np.abs(gain))[::-1][:4]
old_d = root_of(spec06_product(worst_run))
old_offset = values(old_d, "pressure_Pa") - values(old_d, "pressure_label_Pa")
lines.append(f"the largest remaining identity, {worst_run} at 5e4 "
             f"({now[worst_run]['5e4']['pressure_identity_largest_abs']:.3e}): the layers that gain most of "
             f"its offset, with the layer law's scale M |jump in ln T| beside each (SPEC_06's offset at the "
             f"same layer in brackets):")
for i in sorted(top):
    lines.append(f"    layer {i} to {i + 1} ({values(run2, 'pressure_label_Pa')[i]:.4g} to "
                 f"{values(run2, 'pressure_label_Pa')[i + 1]:.4g} Pa): gains {gain[i]:+.3f} Pa "
                 f"[{np.diff(old_offset)[i]:+.3f}], M |jump| {mass[i] * abs(jump[i]):.3f} Pa")

# Away from slope changes, the delivered temperature against SPEC_06.
SLOPE_CHANGES = {"shear_r4_decay12": (1.0e5, 20.0), "shear_r5_decay20": (1.0e5, 700.0),
                 "shear_r6_decay40": (1.0e5, 8000.0), "shear_r7_decay_linp": (1.0e5,),
                 "shear_r7f_decay_linp_fine": (1.0e5,), "shear_r8_increase25": (1.0e5,),
                 "shear_r9_increase50": (1.0e5,)}
fourth, far_lines = True, []
for run, entry in now.items():
    for spacing in ("5e4", "2.5e4"):
        if entry.get(spacing, {}).get("status") != "completed":
            continue
        new_T = values(root_of(entry[spacing]["product"]), "temperature_K")
        old_T = values(root_of(spec06_product(run, spacing)), "temperature_K")
        diff = np.abs(new_T - old_T)
        if run in SLOPE_CHANGES:
            grid = FIELDS[f"forward/{run}/inputs/{run}_wind.nc"].ln_pressure
            width = float(np.max(np.diff(grid)))
            near = np.zeros(diff.shape, dtype=bool)
            for change in SLOPE_CHANGES[run]:
                near |= np.abs(np.log(p_lab) - np.log(change)) <= width
            far = float(np.max(diff[~near])) if np.any(~near) else 0.0
            inner = float(np.max(diff[near])) if np.any(near) else 0.0
            fourth = fourth and far <= 0.5
            far_lines.append(f"    {run} {spacing}: farther than one wind interval from a slope change "
                             f"{far:.3f} K (bound 0.5), within one interval {inner:.3f} K (reported)")
        else:
            fourth = fourth and float(np.max(diff)) <= 0.1
            far_lines.append(f"    {run} {spacing}: every level {float(np.max(diff)):.3e} K (bound 0.1)")
ok = ok and fourth
lines.append("the delivered temperature against SPEC_06:")
lines.extend(far_lines)
record(7, "no steps: run 7f's isobar ten times smoother with no pattern at the wind's latitude nodes; "
          "run 7f's lowest 50 km continuous at every layer straddling a wind node; every identity not "
          "above SPEC_06's and the spacings within 10 percent, the largest read by the layer law; the "
          "delivered temperature within 0.1 K of SPEC_06 without slope changes and 0.5 K away from them",
       ok, NEWLINE.join(lines))

# ---------------------------------------------------------------------------------------------
# 8. The quadrature, three ways
# ---------------------------------------------------------------------------------------------
_two_point = lk._piecewise


def _gauss(kernel, isobar_map, index, start, stop, offsets, weights, extra=None):
    """`_piecewise` with a Gauss rule given by `offsets` and `weights` on [-1, 1], and, with
    `extra`, further breakpoints (the column's map knots)."""
    start, stop = np.broadcast_arrays(np.asarray(start, dtype="float64"),
                                      np.asarray(stop, dtype="float64"))
    phi = float(kernel.mesh.latitude_rad[index])
    nodes = kernel.field.ln_pressure
    breaks = [isobar_map.geopotential_at(index, n) for n in nodes]
    if extra is not None:
        breaks += list(extra(index))
    breaks = np.asarray(breaks, dtype="float64")
    a, b = np.minimum(start, stop), np.maximum(start, stop)
    inner = np.where((breaks[:, None] > a.ravel()[None, :]) & (breaks[:, None] < b.ravel()[None, :]),
                     breaks[:, None], stop.ravel()[None, :]).reshape((breaks.size,) + start.shape)
    inner = np.sort(inner, axis=0)
    inner = np.where(start > stop, inner[::-1], inner)
    edges = np.concatenate([start[None], inner, stop[None]])
    mid, half = 0.5 * (edges[1:] + edges[:-1]), 0.5 * (edges[1:] - edges[:-1])
    total = 0.0
    for o, w in zip(offsets, weights):
        points = mid + half * o
        total = total + w * half * lk.shear_at(kernel, phi, points, np.exp(isobar_map.ln_p_at(index, points)))
    return total.sum(axis=0)


GL4 = np.polynomial.legendre.leggauss(4)
GL2 = np.polynomial.legendre.leggauss(2)
variants = {
    "two points (as delivered)": None,
    "four points": lambda k, m, i, s, e: _gauss(k, m, i, s, e, GL4[0], GL4[1]),
    "two points and the map's knots": lambda k, m, i, s, e: _gauss(
        k, m, i, s, e, GL2[0], GL2[1], extra=lambda j: m.geopotential_m2s2[j]),
}
from casspian.forward import transfer as forward_transfer  # noqa: E402

namelist = Path(f"forward/{R7F}/{R7F}.toml")
quad = {}
original = relaxed()
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for label, variant in variants.items():
            lk._piecewise = _two_point if variant is None else variant
            copy = HERE / "quadrature" / label.split()[0] / R7F
            if copy.exists():
                shutil.rmtree(copy)
            (copy / "inputs").mkdir(parents=True)
            for path in (Path("forward") / R7F / "inputs").glob("*.nc"):
                shutil.copy2(path, copy / "inputs" / path.name)
            text = namelist.read_text(encoding="utf-8").replace('"../../', '"../../../../../')
            text = text.replace("figures = true", "figures = false")
            (copy / namelist.name).write_bytes(text.encode("utf-8"))
            started = time.time()
            product = forward_transfer.run(copy / namelist.name).product
            quad[label] = (values(root_of(product), "temperature_K"), time.time() - started)
finally:
    lk._piecewise = _two_point
    ctl._refuse_dirty_commit = original
base = quad["two points (as delivered)"][0]
lines = [f"{label}: largest |dT| from the model as delivered {float(np.max(np.abs(T - base))):.3e} K, "
         f"{seconds:.0f} s" for label, (T, seconds) in quad.items()]
record(8, "run 7f at 5e4 three ways: the model as delivered within 0.1 K of four Gauss points and of two "
          "with the map's knots added, at every level (the knots measured, not adopted)",
       all(float(np.max(np.abs(T - base))) <= 0.1 for T, _ in quad.values()), NEWLINE.join(lines))

# ---------------------------------------------------------------------------------------------
# 9. The cost
# ---------------------------------------------------------------------------------------------
times = {}
original = relaxed()
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for label, run, baseline in (("the Lindal transfer run", "lindal_transfer", 14.0),
                                     ("run 8", "shear_r8_increase25", 16.5)):
            copy = HERE / "cost" / run
            if copy.exists():
                shutil.rmtree(copy)
            (copy / "inputs").mkdir(parents=True)
            for path in (Path("forward") / run / "inputs").glob("*.nc"):
                shutil.copy2(path, copy / "inputs" / path.name)
            text = (Path("forward") / run / f"{run}.toml").read_text(encoding="utf-8")
            text = text.replace('"../../', '"../../../../').replace("figures = true", "figures = false")
            (copy / f"{run}.toml").write_bytes(text.encode("utf-8"))
            started = time.time()
            forward_transfer.run(copy / f"{run}.toml")
            times[label] = (time.time() - started, baseline)
finally:
    ctl._refuse_dirty_commit = original
record(9, "the Lindal transfer run and run 8 each no more than twice SPEC_06's time (14 s and 16.5 s), "
          "figures off",
       all(t <= 2 * b for t, b in times.values()),
       "; ".join(f"{label} {t:.1f} s against {b} s, ratio {t / b:.2f}" for label, (t, b) in times.items()))

failed = [r for r in results if not r[2]]
summary = f"{len(results) - len(failed)} of {len(results)} checks pass, {time.time() - _started:.0f} s"
print(summary)
(HERE / "output.txt").write_text(
    NEWLINE.join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}" + NEWLINE + "        "
                 + str(detail).replace(NEWLINE, NEWLINE + "        ")
                 for n, d, ok, detail in results) + NEWLINE + NEWLINE + summary + NEWLINE,
    encoding="utf-8")
sys.exit(1 if failed else 0)
