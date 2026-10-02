"""Acceptance checks for SPEC_05 v0.6 Step 2: the column march, faster (decision R, first part).

This step changes speed and nothing else, so its acceptance is that nothing else changed. The
script builds candidate copies of all 18 registered files from the working tree, which carries
the step's uncommitted code, and compares every variable of every group with the committed files.

**The candidates, and the relaxation named here** (SPEC_05 §2a, as SPEC_03 and SPEC_04 built
theirs). The registered files and everything they are built from (`data_static/`, the Lindal
reduction's control files and raw bundle, the two runs' control files) are copied into
`reports/step05_2/<stage>/tree/` in the repository's own layout, so every relative path resolves
inside the copy, and the chain is rebuilt there in the order `tests/step04_0/sweep.py` uses. The
working tree is not clean, so every candidate carries `casspian_git_commit` ending in `-dirty`, and
`refrac` and `forward` refuse a `-dirty` input. The single refusal point,
`casspian.lib.control._refuse_dirty_commit`, is relaxed for this script and nowhere else; nothing
under `occul_data/` or `forward/` is written. The comparison is by value, so the stamps do not
matter.

`--stage` names the output directory and the record: `fix1` after fix 1, `fix2` (the default) after
fix 2. Run from the repository root.
"""

import shutil
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import control as ctl

STAGE = sys.argv[sys.argv.index("--stage") + 1] if "--stage" in sys.argv else "fix2"
HERE = Path("reports/step05_2") / STAGE
TREE = HERE / "tree"
#: REVIEW_05_step0, carried into SPEC_05 §2a check 2.
RELATIVE_BOUND = 1e-14
#: The registered files: the reduction chain, the two runs' inputs and both products.
REGISTERED = (
    "occul_data/lindal/raw/lindal_raw.nc",
    "occul_data/lindal/lindal_gravity.nc",
    "occul_data/lindal/lindal_rotation.nc",
    "occul_data/lindal/lindal_wind.nc",
    "occul_data/lindal/lindal_composition.nc",
    "occul_data/lindal/lindal_thermo.nc",
    "occul_data/lindal/lindal_geodesy.nc",
    "occul_data/lindal/lindal_refractivity.nc",
    "forward/lindal_closure/inputs/lindal_closure_gravity.nc",
    "forward/lindal_closure/inputs/lindal_closure_rotation.nc",
    "forward/lindal_closure/inputs/lindal_closure_wind.nc",
    "forward/lindal_closure/inputs/lindal_closure_composition.nc",
    "forward/lindal_closure/output/lindal_closure_profile.nc",
    "forward/lindal_transfer/inputs/lindal_transfer_gravity.nc",
    "forward/lindal_transfer/inputs/lindal_transfer_rotation.nc",
    "forward/lindal_transfer/inputs/lindal_transfer_wind.nc",
    "forward/lindal_transfer/inputs/lindal_transfer_composition.nc",
    "forward/lindal_transfer/output/lindal_transfer_profile.nc",
)
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def copy_tree():
    """The sources and control files of the registered files, in the repository's layout."""
    if HERE.exists():
        shutil.rmtree(HERE)
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


def build_candidates():
    """The order of `tests/step04_0/sweep.py`, on the copy. Returns seconds per stage."""
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
    stages = [
        ("reduction inputs", lambda: (build_raw.build(lindal / "raw", lindal / "raw" / "lindal_raw.nc"),
                                      build_gravity.build(build), build_rotation.build(build),
                                      build_wind.build(build), build_composition.build(build),
                                      build_inputs.build(build))),
        ("kind N (refrac)", lambda: refrac_product.build_product(lindal / "lindal_reduction.toml")),
        ("closure inputs", lambda: run_inputs.build(closure / "lindal_closure_build.toml")),
        ("closure production", lambda: forward_production.run(closure / "lindal_closure.toml")),
        ("transfer inputs", lambda: run_inputs.build(transfer / "lindal_transfer_build.toml")),
        ("transfer run", lambda: forward_transfer.run(transfer / "lindal_transfer.toml")),
    ]
    seconds = {}
    for name, action in stages:
        started = time.perf_counter()
        action()
        seconds[name] = time.perf_counter() - started
    return seconds


def compare(registered: Path, candidate: Path):
    """Every variable of every group: `(variables, bit-identical, outside the bound, worst)`."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = xr.open_datatree(registered, engine="netcdf4").load()
        b = xr.open_datatree(candidate, engine="netcdf4").load()
    variables, identical, outside, worst = 0, 0, [], []
    for node in a.subtree:
        x = node.to_dataset(inherit=False)
        y = b[node.path].to_dataset(inherit=False) if node.path in [n.path for n in b.subtree] else None
        for name in x.variables:
            variables += 1
            if y is None or name not in y.variables:
                outside.append(f"{node.path}:{name} missing from the candidate")
                continue
            xa, ya = np.asarray(x[name].values), np.asarray(y[name].values)
            if xa.shape != ya.shape:
                outside.append(f"{node.path}:{name} shape {ya.shape} against {xa.shape}")
                continue
            if xa.dtype.kind not in "fc":
                if np.array_equal(xa, ya):
                    identical += 1
                else:
                    outside.append(f"{node.path}:{name} differs")
                continue
            if np.array_equal(xa, ya, equal_nan=True):
                identical += 1
                continue
            xa, ya = xa.astype("float64"), ya.astype("float64")
            d = np.abs(ya - xa)
            d[np.isnan(xa) & np.isnan(ya)] = 0.0
            nonzero = xa != 0
            rel = float(np.nanmax(d[nonzero] / np.abs(xa[nonzero]))) if nonzero.any() else 0.0
            ulp = float(np.nanmax(d / np.spacing(np.abs(xa))))
            at = np.unravel_index(int(np.nanargmax(d)), d.shape)
            worst.append(f"{node.path}:{name} at index {tuple(int(i) for i in at)}, largest relative "
                         f"{rel:.3e}, {ulp:g} ulp")
            if rel > RELATIVE_BOUND or np.any(d[~nonzero] != 0) or np.any(np.isnan(xa) != np.isnan(ya)):
                outside.append(worst[-1])
    a.close()
    b.close()
    return variables, identical, outside, worst


print(f"stage {STAGE}: candidates under {TREE}")
copy_tree()
_original_refuse = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(attrs.get("casspian_git_commit", ""))
print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this script only\n")
started = time.perf_counter()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    seconds = build_candidates()
ctl._refuse_dirty_commit = _original_refuse
built = [name for name in REGISTERED if (TREE / name).exists()]
record(1, "the 18 registered files built as candidates from the working tree, in a copy, the -dirty "
          "refusal relaxed for this script only",
       len(built) == len(REGISTERED),
       f"built {len(built)} of {len(REGISTERED)} in {time.perf_counter() - started:.1f} s; "
       + "; ".join(f"{name} {value:.1f} s" for name, value in seconds.items()))

lines, total, total_identical, total_outside = [], 0, 0, []
for name in REGISTERED:
    variables, identical, outside, worst = compare(Path(name), TREE / name)
    total += variables
    total_identical += identical
    total_outside += outside
    lines.append(f"{name}: {variables} variables, {identical} bit-identical"
                 + (f"; not bit-identical: {'; '.join(worst)}" if worst else ""))
record(2, "every variable of every group of the 18 is bit-identical to the committed file, or within "
          f"{RELATIVE_BOUND:g} relative where it is not (where, the largest relative difference and "
          "the largest in units in the last place printed)",
       not total_outside,
       f"{total} variables compared, {total_identical} bit-identical, {total - total_identical} not; "
       f"outside the bound: {total_outside or 'none'}\n" + "\n".join(lines))

# ---------------------------------------------------------------------------
# 3. A refusal inside the all-columns march names the latitude and the level
# ---------------------------------------------------------------------------
from casspian.forward import transfer as tr  # noqa: E402
from casspian.lib import io as cio  # noqa: E402
from casspian.lib import mesh as lm  # noqa: E402
from casspian.lib.windfield import WindField  # noqa: E402

handle = cio.read(Path("forward/lindal_transfer/inputs/lindal_transfer_wind.nc"), "wind")
try:
    field = WindField(handle.load())
finally:
    handle.close()
mesh = lm.build_mesh(np.radians([10.0, 30.0]), np.array([-1.0e5, 1.0e5]), np.radians(1.0), 5.0e4)
n = mesh.latitude_rad.size
# An isobar map whose pressure passes the wind grid's bottom, 1e6 Pa, below Phi = -5e4 m2/s2.
knots = np.array([-2.0e5, 2.0e5])
crossing = tr.IsobarMap(geopotential_m2s2=tuple(knots for _ in range(n)),
                        ln_pressure=tuple(np.log([4.0e6, 1.0e3]) for _ in range(n)))
message = None
try:
    tr._wind_table(mesh, crossing, field)
except ValueError as exc:
    message = str(exc)
names = message is not None and "latitude" in message and "geopotential" in message and " Pa" in message
record(3, "a point outside the wind grid, met inside the all-columns march, is refused with the "
          "latitude and the level named",
       names, f"mesh {n} latitudes by {mesh.geopotential_m2s2.size} geopotential nodes; an isobar map "
              f"running from 4e6 to 1e3 Pa\nrefusal: {message}")

failed = [r for r in results if not r[2]]
summary = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(summary)
(Path("reports/step05_2") / f"output_{STAGE}.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{summary}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
