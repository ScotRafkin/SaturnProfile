"""Acceptance checks for SPEC_08 v0.3 Step 1: kind W holds what the data hold.

The registered winds and products are read, never written. Every file this script makes is under
`reports/step08_1/`: the check 2 files and the check 5 sources through `lib.io.write`, and the
check 3 and 4 copies of the registered winds through `xarray.to_netcdf` with their attributes kept
(SPEC_08 section 6 ruling 5), so that they carry the registered commit and the dirty refusal does
not fire before the refusal under test. The loaders are pointed at the copies through the loaded
manifest and namelists (`dataclasses.replace`), since the copies cannot sit in the run directories.

Every check prints its measured values. Run from the repository root.
"""

import dataclasses
import subprocess
import sys
import warnings
from pathlib import Path
from types import MappingProxyType

import matplotlib
matplotlib.use("Agg")
import numpy as np
import xarray as xr

from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib.control import ControlFileError
from casspian.lib.schema import CasspianSchemaError
from casspian.tools.plots import figures_profile
from casspian.tools.plots.render import render as render_file
from casspian.tools.wind import shear

HERE = Path("reports/step08_1")
COPIES = HERE / "copies"
FILES = HERE / "files"
CASES = HERE / "cases"
FIGURES = HERE / "figures"
for directory in (COPIES, FILES, CASES, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)
REDUCTION_WIND = Path("occul_data/lindal/lindal_wind.nc")
WINDS = [REDUCTION_WIND,
         Path("forward/lindal_closure/inputs/lindal_closure_wind.nc"),
         Path("forward/lindal_closure/inputs/lindal_closure_wind_source.nc"),
         Path("forward/lindal_transfer/inputs/lindal_transfer_wind.nc"),
         Path("forward/lindal_transfer/inputs/lindal_transfer_wind_source.nc")]
TRANSFER_PRODUCT = Path("forward/lindal_transfer/output/lindal_transfer_profile.nc")
RELATIVE = 1e-12
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def read_closed(path, kind="wind"):
    handle = cio.read(path, kind)
    try:
        return handle.load()
    finally:
        handle.close()


def raw(path):
    """The file's values and attributes as on disk, with no schema check, and the file released."""
    handle = xr.open_dataset(path, engine="netcdf4")
    try:
        return handle.load()
    finally:
        handle.close()


def write_raw(dataset, path):
    """Write as `lib.io.write` encodes, without validating or stamping: the attributes are kept."""
    encoding = {str(n): {"_FillValue": np.nan} for n, v in dataset.data_vars.items()
                if np.issubdtype(v.dtype, np.floating)}
    dataset.to_netcdf(path, mode="w", format="NETCDF4", engine="netcdf4", encoding=encoding)
    return path


def array_equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and bool(np.array_equal(a, b, equal_nan=a.dtype.kind in "fc"))


def refusal(action):
    try:
        action()
    except (ControlFileError, CasspianSchemaError) as exc:
        return str(exc)
    return None


manifest = ctl.read_reduction_manifest("occul_data/lindal/lindal_reduction.toml")
closure = ctl.read_run_namelist("forward/lindal_closure/lindal_closure.toml")
transfer = ctl.read_run_namelist("forward/lindal_transfer/lindal_transfer.toml")
LOADERS = (("load_reduction_inputs", manifest, ctl.load_reduction_inputs),
           ("load_run_inputs, closure", closure, ctl.load_run_inputs),
           ("load_run_inputs, transfer", transfer, ctl.load_run_inputs))


def pointed_at(obj, wind_path):
    """The loaded manifest or namelist with its wind replaced by `wind_path`."""
    return dataclasses.replace(
        obj, inputs=MappingProxyType({**obj.inputs, "wind": Path(wind_path).resolve()}))


# ---------------------------------------------------------------------------
# 1. Nothing moves
# ---------------------------------------------------------------------------
lines, ok = [], True
tracked = [str(p) for p in WINDS] + [str(TRANSFER_PRODUCT),
                                     "forward/lindal_closure/output/lindal_closure_profile.nc",
                                     "occul_data/lindal/lindal_refractivity.nc"]
status = subprocess.run(["git", "status", "--porcelain", "--", *tracked], capture_output=True,
                        text=True).stdout.strip()
ok &= status == ""
lines.append(f"git status of the {len(tracked)} registered winds and products: "
             f"{status or 'unmodified'}")
for path in WINDS:
    message = refusal(lambda: read_closed(path))
    ok &= message is None
    lines.append(f"{path}: reads as kind W {message is None}{'' if message is None else ': ' + message}")
for label, obj, loader in LOADERS:
    message = refusal(lambda: loader(obj))
    ok &= message is None
    lines.append(f"{label}: the registered inputs admitted {message is None}"
                 f"{'' if message is None else ': ' + message}")
record(1, "the registered winds and products are unmodified, read as kind W, and every loader admits "
          "the registered inputs", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 2. A file holds what the data hold (section 6 ruling 11)
# ---------------------------------------------------------------------------
base = read_closed(REDUCTION_WIND)
latitude = np.asarray(base["latitude_planetocentric_deg"].values, dtype="float64")
pressure = np.asarray(base["pressure_Pa"].values, dtype="float64")
p_ref = float(base["reference_level_pressure_Pa"])
k_ref = int(np.flatnonzero(pressure == p_ref)[0])
inner = np.flatnonzero(np.abs(latitude) < 85.0)
one = base.isel(latitude_planetocentric=inner, pressure=[k_ref])
three = base.isel(latitude_planetocentric=inner, pressure=[k_ref - 1, k_ref, k_ref + 1])
# A shear that is not zero, so that array equality on read means something: u_reference times
# 0.1 ln(p / p_ref), zero at the reference level.
three = three.copy(deep=True)
three["u_shear_ms"].values[...] = (np.asarray(three["u_reference_ms"].values)[:, None]
                                   * 0.1 * np.log(np.asarray(three["pressure_Pa"].values) / p_ref)[None, :])
FORMS = {
    "reference_only": one.drop_vars(["u_total_ms", "u_total_uncertainty_ms", "u_shear_ms"]),
    "reference_and_shear": three.drop_vars(["u_total_ms", "u_total_uncertainty_ms"]),
    "reference_and_total": one.drop_vars(["u_shear_ms"]),
}
lines, ok, written2 = [], True, {}
for name, dataset in FORMS.items():
    path = FILES / f"step08_1_{name}_wind.nc"
    message = refusal(lambda: cio.write(path, dataset, "wind", created_by="tests/step08_1"))
    if message is None:
        back = read_closed(path)
        unequal = [v for v in dataset.variables if not array_equal(back[v].values, dataset[v].values)]
        good = not unequal and set(back.variables) == set(dataset.variables)
        written2[name] = path
    else:
        good, unequal = False, message
    ok &= good
    lines.append(f"{name}: {latitude[inner].size} latitudes {latitude[inner].min():g} to "
                 f"{latitude[inner].max():g} deg, {dataset.sizes['pressure']} pressure level(s); parts "
                 f"{sorted(v for v in ('u_total_ms', 'u_reference_ms', 'u_shear_ms') if v in dataset)}"
                 f"; written and read back, differing {unequal or 'none'}")
no_reference = write_raw(raw(written2.get("reference_and_shear", REDUCTION_WIND)).drop_vars("u_reference_ms"),
                         FILES / "step08_1_shear_without_reference_wind.nc")
message = refusal(lambda: read_closed(no_reference))
refused = message is not None and "u_reference_ms" in message
ok &= refused
lines.append(f"the shear and the level without the reference wind, on read: "
             f"{'refused' if refused else 'NOT REFUSED BY NAME'}\n    {message}")
record(2, "a file holds what the data hold: the reference wind alone, with the shear on three levels, "
          "and with the total on one level, each short of both poles, read back array-equal; a file "
          "without the reference wind refused on read, by name", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 3. The model takes only what it can run on (section 6 ruling 4)
# ---------------------------------------------------------------------------
DROPS = (("u_total_ms", ["u_total_ms", "u_total_uncertainty_ms"], ControlFileError),
         ("reference_level_pressure_Pa", ["reference_level_pressure_Pa"], CasspianSchemaError),
         ("u_reference_ms", ["u_reference_ms"], CasspianSchemaError))
lines, ok = [], True
for label, obj, loader in LOADERS:
    source = raw(obj.inputs["wind"])
    for name, drop, expected in DROPS:
        copy = write_raw(source.drop_vars(drop), COPIES / f"{obj.inputs['wind'].stem}_no_{name}.nc")
        try:
            loader(pointed_at(obj, copy))
            kind, message = None, "ADMITTED"
        except (ControlFileError, CasspianSchemaError) as exc:
            kind, message = type(exc), str(exc)
        good = kind is expected and name in message
        ok &= good
        where = "by the loader" if expected is ControlFileError else "on read, by the schema"
        lines.append(f"{label}, without {name}: {'refused ' + where if good else 'FAIL'}\n    {message}")
record(3, "each loader refuses a copy without u_total_ms, naming it; a copy without the level or the "
          "reference wind is refused on read by the schema, naming it", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 4. The polar rule at the model
# ---------------------------------------------------------------------------
lines, ok = [], True
for label, obj, loader in LOADERS:
    source = raw(obj.inputs["wind"])
    north = np.abs(np.asarray(source["latitude_planetocentric_deg"].values) - 90.0) <= 1e-9
    source["u_total_ms"].values[north, :] = 1e-6
    copy = write_raw(source, COPIES / f"{obj.inputs['wind'].stem}_north_pole.nc")
    read_message = refusal(lambda: read_closed(copy))
    message = refusal(lambda: loader(pointed_at(obj, copy)))
    good = read_message is None and message is not None and "+90" in message
    ok &= good
    lines.append(f"{label}: the copy reads through lib.io {read_message is None}; the loader "
                 f"{'refuses it' if message else 'ADMITS IT'}\n    {message}")
record(4, "a copy with the north pole at 1e-6 m/s reads through lib.io and is refused by each loader",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 5. The shear tool (section 2 deliverable 5, section 6 ruling 2)
# ---------------------------------------------------------------------------
SHEAR_SOURCE = Path("forward/lindal_transfer/inputs/lindal_transfer_wind.nc")
source = read_closed(SHEAR_SOURCE)
s_pressure = np.asarray(source["pressure_Pa"].values, dtype="float64")
u_s = np.asarray(source["u_total_ms"].values, dtype="float64")[:, int(np.flatnonzero(s_pressure == 1e5)[0])]
no_shear = cio.write(CASES / "no_shear_source_wind.nc", source.drop_vars("u_shear_ms"), "wind",
                     created_by="tests/step08_1")
no_total = cio.write(CASES / "no_total_source_wind.nc",
                     source.drop_vars(["u_total_ms", "u_total_uncertainty_ms"]), "wind",
                     created_by="tests/step08_1")


def control_file(name, source_path):
    path = CASES / f"{name}.toml"
    path.write_text("\n".join([
        "[shear]", f'source = "{Path(source_path).name}"', f'output = "{name}_wind.nc"',
        'case = "decay_above"', "shear_reference_pressure_Pa = 100000.0", 'shape = "linear_ln_p"',
        "stop_pressure_Pa = 700.0", "stop_fraction = 0.0"]) + "\n", encoding="utf-8")
    return path


lines, ok = [], True
out = read_closed(shear.build(control_file("decay_no_shear", no_shear)))
x = np.log(s_pressure / 1e5) / np.log(700.0 / 1e5)
expected = u_s[:, None] * np.where(x <= 0, 1.0, np.where(x >= 1, 0.0, 1.0 - x))[None, :]
total = np.asarray(out["u_total_ms"].values, dtype="float64")
rel = float(np.max(np.abs(total - expected)) / max(float(np.max(np.abs(expected))), 1.0))
no_part = "u_shear_ms" not in out.variables
kept = array_equal(out["u_reference_ms"].values, source["u_reference_ms"].values)
ok &= rel <= RELATIVE and no_part and kept
lines.append(f"input without u_shear_ms, the case of SPEC_05 Step 1 check 5 (decay_above, linear_ln_p, "
             f"p_s = 1e5 Pa, p_stop = 700 Pa, f = 0): largest departure of u_total from u_s F(p) "
             f"{rel:.3e} relative (bound 1e-12); no u_shear_ms in the output {no_part}; u_reference_ms "
             f"array-equal to the input's {kept}")
target = CASES / "decay_no_total_wind.nc"
message = refusal(lambda: shear.build(control_file("decay_no_total", no_total)))
refused = message is not None and "u_total_ms" in message and not target.exists()
ok &= refused
lines.append(f"input without u_total_ms: {'refused by name, nothing written' if refused else 'FAIL'}"
             f"\n    {message}")
record(5, "the shear tool: on an input without u_shear_ms the closed form to 1e-12 relative and no "
          "shear part written; an input without u_total_ms refused by name", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 6. Figures draw what a file has (section 6 rulings 6, 10 and 12)
# ---------------------------------------------------------------------------
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    tree = xr.open_datatree(TRANSFER_PRODUCT, engine="netcdf4").load()
tree.close()
lines, ok = [], True
reference_Pa = float(tree["inputs/wind"].to_dataset(inherit=False)["reference_level_pressure_Pa"])
source_label = f"source wind, assigned to {reference_Pa / 100:g} mbar"
total_label = f"u_total at {reference_Pa / 100:g} mbar"
for name, drop in (("F9_transfer", False), ("F9_without_reference", True)):
    if drop:
        node = tree["inputs/wind"]
        node.dataset = node.to_dataset(inherit=False).drop_vars("u_reference_ms")
    fig, data = figures_profile.figure_9(figures_profile._Profile(tree))
    out_path = FIGURES / f"{name}.png"
    fig.savefig(out_path, dpi=110)
    labels = [line.get_label() for line in fig.axes[0].get_lines()]
    good = (total_label in labels and (source_label in labels) != drop
            and ("u_reference_ms" in data) != drop)
    ok &= good
    lines.append(f"{name}: lines {[l for l in labels if l in (source_label, total_label)]}; "
                 f"u_reference_ms in the figure's data {'u_reference_ms' in data}; written {out_path}")
EXPECT = {"reference_only": None, "reference_and_shear": "u_shear_profile",
          "reference_and_total": "u_profile"}
for name, path in written2.items():
    result = render_file(path, out_dir=FIGURES / name)
    keys = sorted(result.data["wind"])
    profiles = {k.rsplit("_", 2)[0] for k in keys if "_profile_" in k}
    want = EXPECT[name]
    good = len(result.written) == 1 and profiles == ({want} if want else set())
    ok &= good
    lines.append(f"casspian-render {path.name}: right panel {sorted(profiles) or 'blank'} "
                 f"(expected {want or 'blank'}); written {[str(p) for p in result.written]}")
record(6, "F9 for the transfer run (both lines) and without u_reference_ms (the total line only); "
          "casspian-render of check 2's three files draws the parts each carries; for the author to view",
       ok, "\n".join(lines))

failed = [r for r in results if not r[2]]
summary = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(summary)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{summary}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
