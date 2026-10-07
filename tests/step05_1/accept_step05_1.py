"""Acceptance checks for SPEC_05 v0.4 Step 1: the shear tool, `casspian-wind-shear`.

The input is the transfer run's own wind file, `forward/lindal_transfer/inputs/
lindal_transfer_wind.nc`, which is committed: this script reads it and writes only under
`reports/step05_1/`, never over it, and prints its SHA-256 before and after. Each case is built
from a control file this script writes under `reports/step05_1/cases/`, through `build`, the path
the command line and the run build take; check 5's file is written by the console entry itself.
Check 4 needs an input that varies in pressure, which the Lindal wind does not; it is made here,
in memory and then on disk under `reports/step05_1/`.

Every check prints its measured values. Run from the repository root.
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib.control import ControlFileError
from casspian.lib.schema import CasspianSchemaError
from casspian.tools.wind import shear

HERE = Path("reports/step05_1")
HERE.mkdir(parents=True, exist_ok=True)
CASES = HERE / "cases"
SOURCE = Path("forward/lindal_transfer/inputs/lindal_transfer_wind.nc")
#: What `lib.io.write` stamps on every file and so differs from the input in every output.
STAMPED = ("created_by", "created_at", "casspian_git_commit", "history")
REPLACED = ("vertical_structure", "input_hashes")
RELATIVE = 1e-12
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def read_closed(path):
    handle = cio.read(path, "wind")
    try:
        return handle.load()
    finally:
        handle.close()


def toml_value(value):
    if isinstance(value, str):
        return f'"{value}"'
    return repr(float(value)) if isinstance(value, float) else repr(value)


def control_file(name, case, source=SOURCE, **parameters):
    """A control file with one [shear] section, under cases/<name>/, and its output path."""
    directory = CASES / name
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}.toml"
    lines = ["[shear]",
             f'source = "{Path(os.path.relpath(source.resolve(), directory.resolve())).as_posix()}"',
             f'output = "{name}_wind.nc"', f'case = "{case}"']
    lines += [f"{key} = {toml_value(value)}" for key, value in parameters.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_case(name, case, source=SOURCE, **parameters):
    written = shear.build(control_file(name, case, source, **parameters))
    return written, read_closed(written)


def array_equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    return bool(np.array_equal(a, b, equal_nan=a.dtype.kind in "fc"))


def worst_relative(a, b):
    """Largest |a - b| over the largest |b|: relative to the field's size, as the sum identity is."""
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1.0))


def refusal(action):
    try:
        action()
    except (ControlFileError, CasspianSchemaError) as exc:
        return str(exc)
    return None


hash_before = cio.sha256(SOURCE)
source = read_closed(SOURCE)
pressure = np.asarray(source["pressure_Pa"].values, dtype="float64")
latitude = np.asarray(source["latitude_planetocentric_deg"].values, dtype="float64")
column = int(np.flatnonzero(pressure == 1e5)[0])
u_source = np.asarray(source["u_total_ms"].values, dtype="float64")
u_s = u_source[:, column]
print(f"input {SOURCE}, sha256 {hash_before}")
print(f"wind grid {latitude.size} latitudes x {pressure.size} pressures, {pressure.min():g} to "
      f"{pressure.max():g} Pa; vertical_structure {source.attrs['vertical_structure']!r}\n")

# ---------------------------------------------------------------------------
# 1. identity
# ---------------------------------------------------------------------------
path1, out1 = run_case("identity", "identity")
unequal = [v for v in source.variables if not array_equal(out1[v].values, source[v].values)]
stripped_out, dropped = ctl._strip(out1)
stripped_in, _ = ctl._strip(source)
content = stripped_out.identical(stripped_in)
record(1, "identity: every variable array-equal to the input, and the closure comparison's content "
          "check passes against the input",
       not unequal and set(out1.variables) == set(source.variables) and content,
       f"{len(source.variables)} variables compared, differing {unequal or 'none'}\n"
       f"closure content check (CLOSURE_DROPPED_ATTRIBUTES and every sha256 value dropped, then "
       f"identical): {content}; dropped from the output {list(dropped)}")

# ---------------------------------------------------------------------------
# 2. uniform, c = 1, p_s on a node
# ---------------------------------------------------------------------------
path2, out2 = run_case("uniform_c1", "uniform", shear_reference_pressure_Pa=1e5, scale=1.0)
same_total = array_equal(out2["u_total_ms"].values, u_source)
same_reference = (array_equal(out2["u_reference_ms"].values, source["u_reference_ms"].values)
                  and array_equal(out2["reference_level_pressure_Pa"].values,
                                  source["reference_level_pressure_Pa"].values))
record(2, "uniform, c = 1, p_s = 1e5 Pa: u_total array-equal to the input's; u_reference and the "
          "reference pressure array-equal to the input's",
       same_total and same_reference,
       f"u_total array-equal {same_total} (1e5 Pa is node {column} of the wind grid, so this is also "
       f"the statement that WindField returns a node value exactly); u_reference and "
       f"reference_level_pressure_Pa array-equal {same_reference}")

# ---------------------------------------------------------------------------
# 3. uniform, c = 0
# ---------------------------------------------------------------------------
path3, out3 = run_case("uniform_c0", "uniform", shear_reference_pressure_Pa=1e5, scale=0.0)
total3 = np.asarray(out3["u_total_ms"].values)
zero = bool(np.all(total3 == 0.0))
minus = array_equal(out3["u_shear_ms"].values,
                    np.broadcast_to(-np.asarray(source["u_reference_ms"].values)[:, None], total3.shape))
record(3, "uniform, c = 0: u_total exactly zero; u_shear = -u_reference exactly; the file reads back "
          "through the schema",
       zero and minus,
       f"u_total exactly zero at all {total3.size} nodes: {zero}; u_shear equal to -u_reference "
       f"exactly: {minus}; read back by lib.io.read without refusal: True")

# ---------------------------------------------------------------------------
# 4. uniform, c = 1, p_s between two nodes, on an input sheared in pressure
# ---------------------------------------------------------------------------
reference = np.asarray(source["u_reference_ms"].values, dtype="float64")
synthetic = source.copy(deep=True)
g = 1.0 + 0.1 * np.log(pressure / 1e5) ** 2
synthetic["u_total_ms"].values[...] = reference[:, None] * g[None, :]
synthetic["u_shear_ms"].values[...] = synthetic["u_total_ms"].values - reference[:, None]
synthetic.attrs["vertical_structure"] = "sheared in pressure for the Step 1 acceptance"
synthetic_path = cio.write(HERE / "synthetic_sheared_wind.nc", synthetic, "wind",
                           created_by="tests/step05_1/accept_step05_1.py")
synthetic = read_closed(synthetic_path)
p_between = 1.1e5
k = int(np.searchsorted(pressure, p_between) - 1)
w = (np.log(p_between) - np.log(pressure[k])) / (np.log(pressure[k + 1]) - np.log(pressure[k]))
u_syn = np.asarray(synthetic["u_total_ms"].values, dtype="float64")
expected4 = u_syn[:, k] * (1.0 - w) + u_syn[:, k + 1] * w
path4, out4 = run_case("uniform_between", "uniform", synthetic_path,
                       shear_reference_pressure_Pa=p_between, scale=1.0)
total4 = np.asarray(out4["u_total_ms"].values, dtype="float64")
rel4 = max(worst_relative(total4[:, j], expected4) for j in range(pressure.size))
record(4, "uniform, c = 1, p_s = 1.1e5 Pa between two nodes, on an input sheared in pressure made "
          "here: u_s equals the ln p interpolant at every latitude to 1e-12 relative",
       rel4 <= RELATIVE,
       f"input u_total = u_reference (1 + 0.1 ln(p / 1e5)^2); p_s between nodes {pressure[k]:g} and "
       f"{pressure[k + 1]:g} Pa, ln p weight {w:.6f}\n"
       f"largest departure of u_total from the interpolant, over every column: {rel4:.3e} "
       f"relative (bound 1e-12)")

# ---------------------------------------------------------------------------
# 5. decay_above, linear_ln_p, p_stop = 700 Pa, f = 0, written by the console entry
# ---------------------------------------------------------------------------
parameters5 = dict(shear_reference_pressure_Pa=1e5, shape="linear_ln_p", stop_pressure_Pa=700.0,
                   stop_fraction=0.0)
control5 = control_file("decay_ln_p", "decay_above", **parameters5)
cli = subprocess.run(["casspian-wind-shear", str(control5)], capture_output=True, text=True)
path5 = CASES / "decay_ln_p" / "decay_ln_p_wind.nc"
out5 = read_closed(path5)
x5 = np.log(pressure / 1e5) / np.log(700.0 / 1e5)
F5 = np.where(x5 <= 0, 1.0, np.where(x5 >= 1, 0.0, 1.0 - x5))
expected5 = u_s[:, None] * F5[None, :]
total5 = np.asarray(out5["u_total_ms"].values, dtype="float64")
rel5 = worst_relative(total5, expected5)
below = pressure >= 1e5
above = pressure <= 700.0
one_below = bool(np.all(total5[:, below] == u_s[:, None]))
zero_above = bool(np.all(total5[:, above] == 0.0))
record(5, "decay_above, linear_ln_p, p_s = 1e5, p_stop = 700 Pa, f = 0: u_total equals u_s F(p) to "
          "1e-12 relative; F is 1 at and below p_s and 0 at and above p_stop",
       cli.returncode == 0 and rel5 <= RELATIVE and one_below and zero_above,
       f"casspian-wind-shear exit {cli.returncode}: {cli.stdout.strip()} {cli.stderr.strip()}\n"
       f"largest departure from the closed form {rel5:.3e} relative (bound 1e-12)\n"
       f"u_total = u_s exactly at the {int(below.sum())} levels at and below 1e5 Pa: {one_below}; "
       f"exactly zero at the {int(above.sum())} levels at and above 700 Pa: {zero_above}")

# ---------------------------------------------------------------------------
# 6. decay_above, linear_p, p_stop = 0, f = 0
# ---------------------------------------------------------------------------
path6, out6 = run_case("decay_p", "decay_above", shear_reference_pressure_Pa=1e5, shape="linear_p",
                       stop_pressure_Pa=0.0, stop_fraction=0.0)
total6 = np.asarray(out6["u_total_ms"].values, dtype="float64")
up = pressure < 1e5
rel6 = worst_relative(total6[:, up], u_s[:, None] * (pressure[up] / 1e5)[None, :])
unchanged6 = bool(np.all(total6[:, ~up] == u_s[:, None]))
record(6, "decay_above, linear_p, p_stop = 0, f = 0: u_total = u_s p / p_s above p_s to 1e-12 "
          "relative",
       rel6 <= RELATIVE,
       f"largest departure from u_s p / p_s at the {int(up.sum())} levels above 1e5 Pa: {rel6:.3e} "
       f"relative (bound 1e-12); at and below p_s u_total = u_s exactly: {unchanged6}")

# ---------------------------------------------------------------------------
# 7. increase_below, linear_ln_p, p_stop = 1e6 Pa, f = 1.5
# ---------------------------------------------------------------------------
path7, out7 = run_case("increase_ln_p", "increase_below", shear_reference_pressure_Pa=1e5,
                       shape="linear_ln_p", stop_pressure_Pa=1e6, stop_fraction=1.5)
total7 = np.asarray(out7["u_total_ms"].values, dtype="float64")
x7 = np.log(pressure / 1e5) / np.log(1e6 / 1e5)
F7 = np.where(x7 <= 0, 1.0, np.where(x7 >= 1, 1.5, 1.0 + 0.5 * x7))
rel7 = worst_relative(total7, u_s[:, None] * F7[None, :])
bottom = int(np.flatnonzero(pressure == 1e6)[0])
at_bottom = worst_relative(total7[:, bottom], 1.5 * u_s)
record(7, "increase_below, linear_ln_p, p_stop = 1e6 Pa, f = 1.5: the closed form to 1e-12 relative; "
          "u_total = 1.5 u_s at 1e6 Pa",
       rel7 <= RELATIVE and at_bottom <= RELATIVE,
       f"largest departure from the closed form {rel7:.3e} relative (bound 1e-12); at 1e6 Pa "
       f"{at_bottom:.3e} relative, exactly equal: {bool(np.all(total7[:, bottom] == 1.5 * u_s))}")

# ---------------------------------------------------------------------------
# 8. Every output of checks 2 to 7
# ---------------------------------------------------------------------------
OUTPUTS = {2: (out2, source), 3: (out3, source), 4: (out4, synthetic), 5: (out5, source),
           6: (out6, source), 7: (out7, source)}
poles = np.flatnonzero(np.abs(np.abs(latitude) - 90.0) <= 1e-9)
lines, ok = [], True
for number, (out, given) in OUTPUTS.items():
    pole_zero = bool(np.all(np.asarray(out["u_total_ms"].values)[poles] == 0.0))
    provenance = bool(np.all(np.asarray(out["value_provenance"].values) == shear.PARAMETERIZED))
    flags_kept = array_equal(out["value_provenance"].attrs["flag_values"],
                             given["value_provenance"].attrs["flag_values"])
    source_data = (array_equal(out["u_reference_ms"].values, given["u_reference_ms"].values)
                   and array_equal(out["reference_level_pressure_Pa"].values,
                                   given["reference_level_pressure_Pa"].values))
    differ = sorted(key for key in set(out.attrs) | set(given.attrs)
                    if key not in out.attrs or key not in given.attrs
                    or not array_equal(out.attrs[key], given.attrs[key]))
    unexpected = [key for key in differ if key not in STAMPED + REPLACED]
    variable_attrs = sorted(f"{v}.{key}" for v in given.variables
                            for key in set(out[v].attrs) | set(given[v].attrs)
                            if not array_equal(out[v].attrs.get(key), given[v].attrs.get(key)))
    good = pole_zero and provenance and flags_kept and source_data and not unexpected and not variable_attrs
    ok = ok and good
    lines.append(f"check {number}: read back through the schema; u_total zero at both "
                 f"poles {pole_zero}; value_provenance 2 everywhere {provenance}, flag_values kept "
                 f"{flags_kept}; source data array-equal {source_data}; globals differing {differ}, of "
                 f"which unexpected {unexpected or 'none'}; variable attributes differing "
                 f"{variable_attrs or 'none'}; vertical_structure {out.attrs['vertical_structure']!r}")
record(8, "every output of checks 2 to 7: read back through the schema, u_total exactly zero at both poles, "
          "value_provenance 2 everywhere, source data array-equal, every attribute equal to the "
          "input's except vertical_structure, input_hashes and what lib.io.write stamps",
       ok, f"excepted, as stamped by lib.io.write: {list(STAMPED)}; as replaced by the tool: "
           f"{list(REPLACED)}\n" + "\n".join(lines))

# ---------------------------------------------------------------------------
# 9. construct in memory against the file build writes
# ---------------------------------------------------------------------------
in_memory = shear.construct(source, "decay_above", parameters5)
unequal9 = [v for v in out5.variables if not array_equal(in_memory[v].values, out5[v].values)]
record(9, "construct on the in-memory input returns a dataset whose every variable is array-equal "
          "to the file build writes and reads back, for the case of check 5",
       not unequal9 and set(in_memory.variables) == set(out5.variables),
       f"{len(out5.variables)} variables compared, differing {unequal9 or 'none'}")

# ---------------------------------------------------------------------------
# 9a. The uncertainty
# ---------------------------------------------------------------------------
kept = {number: array_equal(out["u_total_uncertainty_ms"].values,
                            given["u_total_uncertainty_ms"].values)
        for number, (out, given) in OUTPUTS.items()}
path9a, out9a = run_case("decay_ln_p_uncertainty", "decay_above", uncertainty_ms=25.0, **parameters5)
set25 = bool(np.all(np.asarray(out9a["u_total_uncertainty_ms"].values) == 25.0))
record("9a", "u_total_uncertainty_ms array-equal to the input's in every output of checks 2 to 7; "
             "with uncertainty_ms = 25 on the case of check 5 it is 25 at every node",
       all(kept.values()) and set25,
       f"array-equal to the input's, by check: {kept}\n"
       f"with uncertainty_ms = 25: 25 at every node {set25}; long_name "
       f"{out9a['u_total_uncertainty_ms'].attrs.get('long_name')!r}, uncertainty_method "
       f"{out9a['u_total_uncertainty_ms'].attrs.get('uncertainty_method')!r}")

# ---------------------------------------------------------------------------
# 10. The named refusals
# ---------------------------------------------------------------------------
ramp = dict(shear_reference_pressure_Pa=1e5, shape="linear_ln_p", stop_pressure_Pa=700.0,
            stop_fraction=0.0)
REFUSALS = [
    ("unknown case", "r_case", "cylinders", {}),
    ("unknown shape", "r_shape", "decay_above", {**ramp, "shape": "linear_lnp"}),
    ("missing parameter", "r_missing", "decay_above",
     {k: v for k, v in ramp.items() if k != "stop_fraction"}),
    ("a parameter the case does not use", "r_unused", "uniform",
     {"shear_reference_pressure_Pa": 1e5, "stop_fraction": 0.0}),
    ("uncertainty_ms with identity", "r_identity_unc", "identity", {"uncertainty_ms": 25.0}),
    ("p_s outside the wind grid, below its top", "r_ps_top", "uniform",
     {"shear_reference_pressure_Pa": 0.5}),
    ("p_s outside the wind grid, beyond its bottom", "r_ps_bottom", "uniform",
     {"shear_reference_pressure_Pa": 2e6}),
    ("decay_above with p_stop below p_s", "r_side_decay", "decay_above",
     {**ramp, "stop_pressure_Pa": 2e5}),
    ("increase_below with p_stop above p_s", "r_side_increase", "increase_below",
     {**ramp, "stop_pressure_Pa": 5e4}),
    ("linear_ln_p with p_stop = 0", "r_ln_p_zero", "decay_above", {**ramp, "stop_pressure_Pa": 0.0}),
]
lines, ok = [], True
for what, name, case, parameters in REFUSALS:
    control = control_file(name, case, **parameters)
    message = refusal(lambda: shear.build(control))
    written = (CASES / name / f"{name}_wind.nc").exists()
    ok = ok and message is not None and not written
    lines.append(f"{what}: {'refused' if message else 'NOT REFUSED'}, nothing written {not written}"
                 f"\n    {message}")
record(10, "the named refusals, and only those: unknown case, unknown shape, missing parameter, a "
           "parameter the case does not use (uncertainty_ms with identity included), p_s outside "
           "the wind grid, p_stop on the wrong side of p_s, linear_ln_p with p_stop = 0",
       ok, "\n".join(lines))

hash_after = cio.sha256(SOURCE)
print(f"\ninput {SOURCE} sha256 after the run {hash_after}, unchanged {hash_after == hash_before}")
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results)
    + f"\n\ninput sha256 before {hash_before}, after {hash_after}\n"
    + f"{len(results) - len(failed)} of {len(results)} checks pass\n", encoding="utf-8")
sys.exit(1 if failed else 0)
