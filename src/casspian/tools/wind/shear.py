"""Write a kind W file under a named shear hypothesis. SPEC_05 Step 1.

The input is any kind W file; the output is the same file with its total wind rebuilt by a named
case from `u_s(phi)`, the input's own total wind read at the shear reference pressure `p_s`. The
model reads only `u_total_ms` (SPEC_05 section 1.1), so every shear hypothesis is a file this tool
writes, and no model code changes.

Three layers. The parts know nothing about any case: `u_at_pressure` reads `u_s` through
`lib.windfield.WindField`, the model's own interpolant, so `u_s` is what the model would read at
`p_s`; `position_ln_p` and `position_p` give the position `x(p)` between `p_s` and the stop
pressure; `ramp` gives `F(x)`. The case functions build `u_total(phi, p)` from those parts, and
`CASES` maps each case name to its function and its parameters. `assemble` is the one function
every case passes through: where the input carries a shear part it replaces its values with
`u_total - u_reference`, and where the input carries none it writes none (SPEC_08 section 6
ruling 2).

`construct` touches no file (SPEC_05 section 1.5): the Monte Carlo driver will call it later
without writing. `build` reads the control section and the source file, calls `construct`, and
writes.

Replace, never add (SPEC_05 section 1.7). The output carries the input's attributes and auxiliary
variables as they are. The only values overwritten are those that would otherwise be false about
the new wind: `vertical_structure` (the case name), `value_provenance` (2, parameterized, in every
cell), and, only when `uncertainty_ms` is given, the uncertainty's values and its two descriptive
attributes. `build` adds `input_hashes`; `lib.io.write` stamps what it stamps on every file.
`identity` returns its input untouched.

**`lindal_iris`** (SPEC_11 v0.4 Step 2) is the one case that does not rebuild the total on the
input's own grid: it returns the whole wind of `casspian.tools.lindal.lindal_wind.construct`, the
input being its cloud wind, on that tool's own pressure grid and reference level. Its four inputs,
`temperatures`, `composition`, `gravity` and `rotation`, are paths in the control section; `build`
reads them and passes them in memory, so `construct` still touches no file.

What is refused, and only this (SPEC_05 Step 1 check 10, decision N): an unknown case, an unknown
shape, a missing parameter, a parameter the case does not use, an input without `u_total_ms`
(SPEC_08), `p_s` outside the wind grid, the stop pressure on the wrong side of `p_s`, and
`linear_ln_p` with a stop pressure of zero. A non-finite value reaches the poles, where the run's
loader refuses the file.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from casspian.lib import io as cio
from casspian.lib.control import ControlFileError, load_section
from casspian.lib.windfield import WindField

TOOL = "casspian-wind-shear"

#: The parameter keys any case may carry. Which ones a case requires or accepts is in `CASES`.
PARAMETER_KEYS = (
    "shear_reference_pressure_Pa", "scale", "shape", "stop_pressure_Pa", "stop_fraction",
    "uncertainty_ms", "temperatures", "composition", "gravity", "rotation", "top_level_Pa",
    "above_top",
)

SECTION_KEYS = {"source": True, "output": True, "case": True,
                **{key: False for key in PARAMETER_KEYS}}

#: The `lindal_iris` inputs `build` reads into memory, by key and kind (`None` is the IRIS CSV).
INPUT_KEYS = {"temperatures": None, "composition": "composition", "gravity": "gravity",
              "rotation": "rotation"}

PATH_KEYS = ("source", "output", *INPUT_KEYS)

SHAPES = ("linear_ln_p", "linear_p")

#: `value_provenance` of every cell a case builds: the existing meaning `parameterized`.
PARAMETERIZED = 2

UNCERTAINTY_STATEMENT = "set in the shear control file (uncertainty_ms)"


# ---------------------------------------------------------------------------
# The parts
# ---------------------------------------------------------------------------

def u_at_pressure(source, p_s: float) -> np.ndarray:
    """`u_s(phi)`: the source's `u_total` at `p_s` at every latitude node, in the file's order.

    Read through `WindField`, linear in `ln p` between the bracketing nodes, so the value is the
    one the model would read. `WindField`'s refusal outside the grid is the `p_s` refusal.
    """
    field = WindField(source)
    latitude = np.radians(np.asarray(source["latitude_planetocentric_deg"].values, dtype="float64"))
    try:
        return field.wind_at(latitude, np.full(latitude.shape, float(p_s)))
    except ValueError as exc:
        low, high = field.pressure_bounds_Pa
        raise ControlFileError(
            f"shear_reference_pressure_Pa = {p_s!r} is outside the wind grid of the source, "
            f"{low} to {high} Pa ({exc})."
        ) from exc


def position_ln_p(p, p_s: float, p_stop: float) -> np.ndarray:
    """`x(p)` for `linear_ln_p`: `ln(p / p_s) / ln(p_stop / p_s)`."""
    return np.log(np.asarray(p, dtype="float64") / p_s) / np.log(p_stop / p_s)


def position_p(p, p_s: float, p_stop: float) -> np.ndarray:
    """`x(p)` for `linear_p`: `(p - p_s) / (p_stop - p_s)`."""
    return (np.asarray(p, dtype="float64") - p_s) / (p_stop - p_s)


def ramp(x, f: float) -> np.ndarray:
    """`F(x)`: 1 for `x <= 0`, `1 + (f - 1) x` between, `f` for `x >= 1`."""
    x = np.asarray(x, dtype="float64")
    return np.where(x <= 0.0, 1.0, np.where(x >= 1.0, f, 1.0 + (f - 1.0) * x))


POSITIONS = {"linear_ln_p": position_ln_p, "linear_p": position_p}


# ---------------------------------------------------------------------------
# The cases
# ---------------------------------------------------------------------------

def case_uniform(u_s, pressure, parameters) -> np.ndarray:
    """`c u_s(phi)` at every pressure."""
    c = float(parameters.get("scale", 1.0))
    return np.repeat((c * u_s)[:, None], pressure.size, axis=1)


def case_ramp(u_s, pressure, parameters) -> np.ndarray:
    """`u_s(phi) F(p)`, the one definition of `decay_above` and `increase_below`."""
    p_s = float(parameters["shear_reference_pressure_Pa"])
    p_stop = float(parameters["stop_pressure_Pa"])
    x = POSITIONS[parameters["shape"]](pressure, p_s, p_stop)
    return u_s[:, None] * ramp(x, float(parameters["stop_fraction"]))[None, :]


@dataclass(frozen=True)
class Case:
    """A case's function and the parameters it requires and accepts."""

    function: Callable | None
    required: tuple[str, ...] = ()
    optional: tuple[str, ...] = ()
    #: For the two ramp cases: which side of `p_s` the stop pressure must lie on, +1 or -1.
    stop_side: int = 0
    #: True for a case whose function returns the whole output dataset from the source.
    whole: bool = False


def case_lindal_iris(source, parameters):
    """The Lindal-only wind of SPEC_11 on the cloud wind `source`, its inputs in memory."""
    from casspian.tools.lindal import lindal_wind

    return lindal_wind.construct(
        source, parameters["temperatures"], parameters["composition"], parameters["gravity"],
        parameters["rotation"], float(parameters["top_level_Pa"]), parameters["above_top"]).dataset


RAMP_PARAMETERS = ("shear_reference_pressure_Pa", "shape", "stop_pressure_Pa", "stop_fraction")

CASES = {
    "identity": Case(None),
    "uniform": Case(case_uniform, ("shear_reference_pressure_Pa",), ("scale", "uncertainty_ms")),
    "decay_above": Case(case_ramp, RAMP_PARAMETERS, ("uncertainty_ms",), stop_side=-1),
    "increase_below": Case(case_ramp, RAMP_PARAMETERS, ("uncertainty_ms",), stop_side=+1),
    "lindal_iris": Case(case_lindal_iris, (*INPUT_KEYS, "top_level_Pa", "above_top"), whole=True),
}


def check_parameters(case: str, parameters: dict) -> Case:
    """Refuse the named faults of SPEC_05 Step 1 check 10 that need no grid; return the case."""
    if case not in CASES:
        raise ControlFileError(f"case = {case!r} is not a shear case; the cases are {list(CASES)}.")
    spec = CASES[case]
    missing = [key for key in spec.required if key not in parameters]
    if missing:
        raise ControlFileError(f"case {case!r} requires {missing}, which the section does not give.")
    unused = sorted(set(parameters) - set(spec.required) - set(spec.optional))
    if unused:
        raise ControlFileError(
            f"case {case!r} does not use {unused}; a key no code reads is refused as a typo. "
            f"This case takes {list(spec.required + spec.optional) or 'no parameters'}.")
    if "shape" in spec.required:
        shape = parameters["shape"]
        if shape not in SHAPES:
            raise ControlFileError(f"shape = {shape!r} is not a shape; the shapes are {list(SHAPES)}.")
        p_s = float(parameters["shear_reference_pressure_Pa"])
        p_stop = float(parameters["stop_pressure_Pa"])
        if spec.stop_side * (p_stop - p_s) <= 0.0:
            where = "below p_s (a higher pressure)" if spec.stop_side > 0 else "above p_s (a lower pressure)"
            raise ControlFileError(
                f"case {case!r} changes the wind {where}: stop_pressure_Pa = {p_stop} must lie "
                f"on that side of shear_reference_pressure_Pa = {p_s}.")
        if shape == "linear_ln_p" and p_stop == 0.0:
            raise ControlFileError(
                "shape 'linear_ln_p' with stop_pressure_Pa = 0 puts the stop at an infinite "
                "distance in ln p and returns the input unchanged; a decay to zero pressure is "
                "shape 'linear_p'.")
    return spec


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def assemble(source, u_total, case: str, uncertainty_ms=None):
    """The output dataset: the source with its total replaced, and its shear where it has one.

    Every case passes through here. Where the source carries `u_shear_ms` its values become
    `u_total - u_reference`; where it carries none, none is written (SPEC_08 section 6 ruling 2).
    Only the values that would be false are overwritten (SPEC_05 section 1.7).
    """
    out = source.copy(deep=True)
    out["u_total_ms"].values[...] = u_total
    if "u_shear_ms" in out.variables:
        reference = np.asarray(source["u_reference_ms"].values, dtype="float64")
        out["u_shear_ms"].values[...] = u_total - reference[:, None]
    out["value_provenance"].values[...] = PARAMETERIZED
    out.attrs["vertical_structure"] = case
    if uncertainty_ms is not None:
        out["u_total_uncertainty_ms"].values[...] = float(uncertainty_ms)
        for name in ("long_name", "uncertainty_method"):
            out["u_total_uncertainty_ms"].attrs[name] = UNCERTAINTY_STATEMENT
    return out


def construct(source, case: str, parameters: dict):
    """The complete output dataset for `case`, from an in-memory kind W dataset. No file is read
    or written (SPEC_05 section 1.5)."""
    parameters = dict(parameters)
    spec = check_parameters(case, parameters)
    if "u_total_ms" not in source.variables:
        raise ControlFileError(
            "the source carries no u_total_ms; every case reads u_s from the total, so a source "
            "without it is refused (SPEC_08 section 2 deliverable 5).")
    if spec.function is None:
        return source.copy(deep=True)
    if spec.whole:
        return spec.function(source, parameters)
    u_s = u_at_pressure(source, float(parameters["shear_reference_pressure_Pa"]))
    pressure = np.asarray(source["pressure_Pa"].values, dtype="float64")
    u_total = spec.function(u_s, pressure, parameters)
    return assemble(source, u_total, case, parameters.get("uncertainty_ms"))


def _read(path, kind):
    """A file into memory under its kind, released; the IRIS CSV (`kind` None) as a table."""
    if kind is None:
        return np.genfromtxt(path, delimiter=",", names=True)
    handle = cio.read(path, kind)
    try:
        return handle.load()
    finally:
        handle.close()


def build(control_path, section: str = "shear") -> Path:
    """Read the control section and its source, construct, and write. Returns the path written."""
    control = load_section(control_path, section, SECTION_KEYS, path_keys=PATH_KEYS)
    source_path, output = Path(control["source"]), Path(control["output"])
    parameters = {key: control[key] for key in PARAMETER_KEYS if key in control}
    inputs = [source_path]
    for key, kind in INPUT_KEYS.items():
        if key in parameters:
            inputs.append(Path(parameters[key]))
            parameters[key] = _read(Path(parameters[key]), kind)
    source = _read(source_path, "wind")
    dataset = construct(source, control["case"], parameters)
    dataset.attrs["input_hashes"] = cio.input_hashes(inputs + [Path(control_path)], output)
    return cio.write(output, dataset, "wind", created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Write a kind W file under a named shear case.")
    parser.add_argument("control", help="path to the TOML control file")
    parser.add_argument("--section", default="shear", help="section to read (default shear)")
    args = parser.parse_args(argv)
    print(f"wrote {build(args.control, args.section)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
