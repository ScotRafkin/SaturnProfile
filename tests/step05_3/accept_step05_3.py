"""Acceptance checks for SPEC_05 v0.6 Step 3: the shear step in every forward run, the new-run tool.

Everything here runs in a copy, `reports/step05_3/tree/`, in the repository's own layout:
`data_static/`, the registered Lindal reduction (`occul_data/lindal/`, which this step does not
reach, copied as it is) and the two runs' migrated control files. The runs' inputs and products
are rebuilt there from the working tree, and compared by value with the registered ones; nothing
under `occul_data/` or `forward/` is written. The working tree is not clean, so every file built
here carries `-dirty`, and the one refusal point, `casspian.lib.control._refuse_dirty_commit`, is
relaxed for this script only, as SPEC_05 §2a's candidates were. The registered files themselves
are rebuilt by `tests/step04_0/sweep.py --runs` after the acceptance commit.

Run from the repository root.
"""

import re
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import production as forward_production
from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib.control import ControlFileError
from casspian.tools.plots import figures_profile
from casspian.tools.run import new_run, run_inputs

HERE = Path("reports/step05_3")
TREE = HERE / "tree"
FORWARD = TREE / "forward"
CLOSURE, TRANSFER = FORWARD / "lindal_closure", FORWARD / "lindal_transfer"
REGISTERED_CLOSURE = Path("forward/lindal_closure/output/lindal_closure_profile.nc")
REGISTERED_TRANSFER = Path("forward/lindal_transfer/output/lindal_transfer_profile.nc")
#: `step04_5` check 11's comparison of the closure production with the registered product.
CLOSURE_COMPARED = ("pressure_Pa", "temperature_K", "number_density_m3", "mean_refractivity_m3",
                    "mean_molar_mass_kg_mol", "refractivity", "radius_m",
                    "height_above_anchor_isobar_m", "geopotential_m2s2")
#: Step 3 check 5: the transfer product's delivered N, p, T, altitude, r0 and pressure identity.
TRANSFER_COMPARED = ("refractivity", "pressure_Pa", "temperature_K", "altitude_m",
                     "reference_surface_radius_m", "pressure_identity_residual")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def root(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(path, engine="netcdf4").load()
    tree.close()
    return tree


def equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and bool(np.array_equal(a, b, equal_nan=a.dtype.kind in "fc"))


def compare(product, registered, names):
    a, b = root(product).to_dataset(inherit=False), root(registered).to_dataset(inherit=False)
    return {name: equal(a[name].values, b[name].values) for name in names}


def copy_tree():
    if HERE.exists():
        shutil.rmtree(HERE)
    shutil.copytree("data_static", TREE / "data_static")
    shutil.copytree("occul_data/lindal", TREE / "occul_data" / "lindal",
                    ignore=shutil.ignore_patterns("figures"))
    for run in ("lindal_closure", "lindal_transfer"):
        (FORWARD / run).mkdir(parents=True)
        for name in (f"{run}.toml", f"{run}_build.toml"):
            shutil.copy2(Path("forward") / run / name, FORWARD / run / name)


def build_and_run(directory, mode):
    run = directory.name
    written = run_inputs.build(directory / f"{run}_build.toml")
    driver = forward_transfer.run if mode == "transfer" else forward_production.run
    return written, Path(driver(directory / f"{run}.toml").product)


copy_tree()
_original_refuse = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(attrs.get("casspian_git_commit", ""))
print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this script only\n")
warnings.simplefilter("ignore")

# ---------------------------------------------------------------------------
# 1. A build file without [shear] is refused, naming the five sections
# ---------------------------------------------------------------------------
case = HERE / "case_no_shear" / "lindal_transfer"
case.mkdir(parents=True)
text = (TRANSFER / "lindal_transfer_build.toml").read_text(encoding="utf-8")
without = re.sub(r"\n\[shear\]\n(?:.*\n)*?\n(?=\[composition\])", "\n", text)
(case / "lindal_transfer_build.toml").write_bytes(without.encode("utf-8"))
try:
    run_inputs.build(case / "lindal_transfer_build.toml")
    message = None
except ControlFileError as exc:
    message = str(exc)
five = "['gravity', 'rotation', 'wind', 'composition', 'shear']"
record(1, "a build file without [shear] is refused, naming the five sections",
       message is not None and five in message and "missing ['shear']" in message
       and "[shear]" not in without,
       f"refusal: {message}")

# ---------------------------------------------------------------------------
# 2. Both migrated build files write both wind files; _wind.nc array-equal to _wind_source.nc
# ---------------------------------------------------------------------------
built = {}
lines, ok = [], True
for directory, mode in ((CLOSURE, "closure"), (TRANSFER, "transfer")):
    run = directory.name
    written = run_inputs.build(directory / f"{run}_build.toml")
    built[run] = written
    source, wind = directory / "inputs" / f"{run}_wind_source.nc", directory / "inputs" / f"{run}_wind.nc"
    a, b = root(source).to_dataset(inherit=False), root(wind).to_dataset(inherit=False)
    unequal = [v for v in a.variables if v not in b.variables or not equal(a[v].values, b[v].values)]
    good = source.exists() and wind.exists() and not unequal and set(a.variables) == set(b.variables)
    ok = ok and good
    lines.append(f"{run}: wrote {[Path(p).name for p in written]}; {len(a.variables)} variables of "
                 f"{wind.name} against {source.name}, differing {unequal or 'none'}")
record(2, "casspian-run-inputs on both migrated build files writes both wind files, the _wind.nc "
          "file array-equal to _wind_source.nc in every variable", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 3. The closure run's content comparison against the anchor's embedded copies
# ---------------------------------------------------------------------------
loaded = {}
for kind in ("composition", "gravity", "rotation", "wind"):
    handle = cio.read(CLOSURE / "inputs" / f"lindal_closure_{kind}.nc", kind)
    loaded[kind] = handle.load()
    handle.close()
anchor = cio.read(TREE / "occul_data/lindal/lindal_refractivity.nc", "refractivity")
anchor_tree = anchor.load()
anchor.close()
comparison = ctl.check_closure_inputs(loaded, anchor_tree)
record(3, "the closure run's content comparison against the anchor's embedded copies passes",
       all(c.identical for c in comparison),
       "\n".join(f"{c.kind}: identical {c.identical}"
                 + (f"; {list(c.differences)}" if c.differences else "")
                 + (f"; dropped {list(c.dropped.get('/', ()))}" if c.dropped else "")
                 for c in comparison))

# ---------------------------------------------------------------------------
# 4. The closure product equals the registered product (step04_5 check 11's comparison)
# ---------------------------------------------------------------------------
closure_product = Path(forward_production.run(CLOSURE / "lindal_closure.toml").product)
same4 = compare(closure_product, REGISTERED_CLOSURE, CLOSURE_COMPARED)
record(4, "the closure product equals the registered product under step04_5 check 11 (the computed "
          "variables by array_equal)",
       all(same4.values()), "; ".join(f"{name} {value}" for name, value in same4.items()))

# ---------------------------------------------------------------------------
# 5. The transfer product equals the registered transfer product
# ---------------------------------------------------------------------------
transfer_product = Path(forward_transfer.run(TRANSFER / "lindal_transfer.toml").product)
same5 = compare(transfer_product, REGISTERED_TRANSFER, TRANSFER_COMPARED)
record(5, "the transfer product's delivered N, p, T, altitude, r0 and pressure identity equal the "
          "registered transfer product's (array_equal)",
       all(same5.values()), "; ".join(f"{name} {value}" for name, value in same5.items()))

# ---------------------------------------------------------------------------
# 6. casspian-new-run
# ---------------------------------------------------------------------------
made = new_run.new_run("step05_3_new", TRANSFER)
only_names, lines = True, []
for old_name, new_name in (("lindal_transfer_build.toml", "step05_3_new_build.toml"),
                           ("lindal_transfer.toml", "step05_3_new.toml")):
    old_text = (TRANSFER / old_name).read_text(encoding="utf-8")
    new_text = (made / new_name).read_text(encoding="utf-8")
    back = new_run.renamed(new_text, "step05_3_new", "lindal_transfer")
    changed = sum(1 for a, b in zip(old_text.splitlines(), new_text.splitlines()) if a != b)
    only_names = only_names and back == old_text and len(old_text.splitlines()) == len(new_text.splitlines())
    lines.append(f"{new_name}: {changed} lines differ from {old_name}, every difference the name "
                 f"(renaming back gives the source exactly: {back == old_text})")
try:
    new_run.new_run("step05_3_new", TRANSFER)
    second = None
except ControlFileError as exc:
    second = str(exc)
_, new_product = build_and_run(made, "transfer")
same6 = compare(new_product, REGISTERED_TRANSFER, TRANSFER_COMPARED)
record(6, "casspian-new-run makes a directory whose two files differ from the source only in the name "
          "fields; a second call is refused; the new run builds and runs to check 5's values",
       only_names and second is not None and all(same6.values()),
       "\n".join(lines) + f"\nsecond call: {second}\n"
       + "the new run against the registered transfer product: "
       + "; ".join(f"{name} {value}" for name, value in same6.items()))

# ---------------------------------------------------------------------------
# 7. F9's left panel: both lines
# ---------------------------------------------------------------------------
zero = new_run.new_run("step05_3_zero", TRANSFER)
build_file = zero / "step05_3_zero_build.toml"
text = build_file.read_text(encoding="utf-8")
text = text.replace('case   = "identity"\n',
                    'case   = "uniform"\nshear_reference_pressure_Pa = 1.0e5\nscale = 0.0\n')
build_file.write_bytes(text.encode("utf-8"))
_, zero_product = build_and_run(zero, "transfer")
lines, ok = [], True
for label, product, expect_zero in (("migrated transfer", transfer_product, False),
                                    ("uniform, c = 0", zero_product, True)):
    fig, data = figures_profile.figure_9(figures_profile._Profile(root(product)))
    out = HERE / f"F9_{label.split(',')[0].replace(' ', '_')}.png"
    fig.savefig(out, dpi=110)
    texts = [line.get_label() for line in fig.axes[0].get_lines()]
    p_mbar = data["reference_level_pressure_Pa"] / 100.0
    labelled = (f"source wind, assigned to {p_mbar:g} mbar" in texts
                and f"u_total at {p_mbar:g} mbar" in texts)
    total, source = data["u_total_at_reference_ms"], data["u_reference_ms"]
    shape = bool(np.all(total == 0.0)) if expect_zero else equal(total, source)
    good = labelled and shape
    ok = ok and good
    lines.append(f"{label}: both lines labelled {labelled}; "
                 + (f"u_total at {p_mbar:g} mbar exactly zero {shape}, the source wind up to "
                    f"{np.max(np.abs(source)):.1f} m/s" if expect_zero else
                    f"u_total at {p_mbar:g} mbar array-equal to the source wind (coincident) {shape}")
                 + f"; written {out}")
record(7, "F9's left panel carries both lines, for the migrated transfer run (coincident) and for a "
          "uniform, c = 0 run (the total on zero)", ok, "\n".join(lines))

ctl._refuse_dirty_commit = _original_refuse
failed = [r for r in results if not r[2]]
summary = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(summary)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{summary}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
