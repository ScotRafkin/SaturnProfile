"""`casspian-wind-fit`: the fit's control file, the relinearization, and the kind W file. SPEC_19.

One pass: the wind from the state; the targets' runs (copies of their base runs with only the wind
changed, under `work/`); the weights; the solve. It repeats until the delivered temperature changes by
less than the tolerance at every target, and the wind and its posterior standard deviation are written.
"""

from __future__ import annotations

import argparse
import os
import shutil
import tomllib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import transfer as forward_transfer
from casspian.lib import io as cio
from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.wind.estimation import solver
from casspian.tools.wind.estimation.gradient import TemperatureGradients, layer_pressures, read_levels, wind_term
from casspian.tools.wind.estimation.prior import Prior
from casspian.tools.wind.estimation.radio import RadioTargets, Target
from casspian.tools.wind.estimation.state import Basis

TOOL = "casspian-wind-fit"
#: 66 levels at about 8 s each; the marches share nothing, so they run apart.
WORKERS = max(1, min(10, (os.cpu_count() or 2) - 2))


def _read(path, kind):
    handle = cio.read(path, kind)
    try:
        return handle.load()
    finally:
        handle.close()


def tree(path):
    handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


def _weight_job(args):
    wind, composition, gravity, rotation, pressure, lat_size = args
    return lw.balance_weight(wind, composition, gravity, rotation, pressure, np.ones((pressure.size, lat_size)))


def _term_job(args):
    return wind_term(*args)


def _run_job(namelist):
    return str(forward_transfer.run(namelist).product)


def copy_run(base, name, wind, work, root) -> Path:
    """A run directory for `base`'s target under `wind`, nothing but the wind changed (as SPEC_17)."""
    source = root / "forward" / base
    target = work / name
    if target.exists():
        shutil.rmtree(target)
    (target / "inputs").mkdir(parents=True)
    for kind in ("composition", "gravity", "rotation"):
        shutil.copyfile(source / "inputs" / f"{base}_{kind}.nc", target / "inputs" / f"{name}_{kind}.nc")
    cio.write(target / "inputs" / f"{name}_wind.nc", wind, "wind", created_by=TOOL)
    text = (source / f"{base}.toml").read_text(encoding="utf-8")
    text = text.replace(f'"{base}', f'"{name}').replace(f"inputs/{base}_", f"inputs/{name}_")
    occul = os.path.relpath(root / "occul_data", target).replace(os.sep, "/")
    text = text.replace('path   = "../../occul_data/', f'path   = "{occul}/')
    text = text.replace("figures = true", "figures = false")
    namelist = target / f"{name}.toml"
    namelist.write_bytes(text.encode("utf-8"))
    return namelist


@dataclass
class Problem:
    """The control file, read: everything one pass needs."""

    name: str
    directory: Path
    root: Path
    basis: Basis
    prior: Prior
    radio: RadioTargets
    gradients: TemperatureGradients
    wind: xr.Dataset
    composition: object
    gravity: object
    rotation: object
    tolerance_K: float
    max_passes: int
    output: Path
    control_path: Path
    lat: np.ndarray = field(init=False)
    pressure: np.ndarray = field(init=False)

    def __post_init__(self):
        self.lat = np.asarray(self.wind["latitude_planetocentric_deg"].values, dtype="float64")
        self.pressure = np.asarray(self.wind["pressure_Pa"].values, dtype="float64")


def read_control(control_path) -> Problem:
    control_path = Path(control_path).resolve()
    with open(control_path, "rb") as handle:
        c = tomllib.load(handle)
    here = control_path.parent
    root = (here / c["fit"]["repository_root"]).resolve()
    base = c["base"]["run"]
    inputs = root / "forward" / base / "inputs"
    wind = _read(inputs / f"{base}_wind.nc", "wind")
    composition = _read(inputs / f"{base}_composition.nc", "composition")
    gravity = _read(inputs / f"{base}_gravity.nc", "gravity")
    rotation = _read(inputs / f"{base}_rotation.nc", "rotation")
    s, p = c["state"], c["prior"]
    basis = Basis.build(c["fit"]["epoch"], s["latitude_step_deg"], s["ln_p_step"], s["bottom_Pa"], s["top_Pa"],
                        float(wind["reference_level_pressure_Pa"]))
    prior = Prior.build(basis, p["amplitude_ms"], p["latitude_correlation_deg"], p["ln_p_correlation"])
    r = c["radio"]
    targets = [Target(t["name"], t["run"], t["latitude_planetocentric_deg"], str(root / t["observed"]),
                      tuple(tuple(e) for e in t.get("extra_K", []))) for t in r["target"]]
    radio = RadioTargets(targets, r["anchor_latitude_planetocentric_deg"], str(root / r["sigma_profile"]),
                         r["correlation_ln_p"], r["floor_K"])
    g = c["gradient"]
    table = np.genfromtxt(root / g["data"], delimiter=",", names=True)
    levels = {int(k): float(v) for k, v in g["levels_Pa"].items()}
    lat_c = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    lat_g = np.asarray(wind["latitude_planetographic_deg"].values, dtype="float64")
    gradients = TemperatureGradients(read_levels(table, levels, lat_c, lat_g, g["sample_deg"], g["equatorial_deg"], g["gap_deg"]))
    return Problem(c["fit"]["name"], here, root, basis, prior, radio, gradients, wind, composition, gravity, rotation,
                   c["solve"]["tolerance_K"], c["solve"]["max_passes"], here / c["fit"]["output"], control_path)


def base_shear(problem: Problem) -> np.ndarray:
    """The base wind's shear on its own grid, `(pressure, latitude)`: its `du/dln p`."""
    x = np.log(problem.pressure)
    order = np.argsort(x)
    u = np.asarray(problem.wind["u_total_ms"].values, dtype="float64")
    out = np.empty((x.size, problem.lat.size))
    out[order] = np.gradient(u[:, order], x[order], axis=1).T
    return out


def wind_dataset(problem: Problem, state, sigma=None) -> xr.Dataset:
    out = problem.wind.copy(deep=True)
    u_ref = np.asarray(out["u_reference_ms"].values, dtype="float64")
    u = u_ref[:, None] + problem.basis.wind_change(state, problem.lat, problem.pressure)
    out["u_total_ms"].values[...] = u
    out["u_shear_ms"].values[...] = u - u_ref[:, None]
    if sigma is not None:
        out["u_total_uncertainty_ms"].values[...] = sigma
        out["u_total_uncertainty_ms"].attrs.update({
            "long_name": "posterior standard deviation of the wind", "uncertainty_kind": "1sigma",
            "uncertainty_method": "the fit's posterior covariance (SPEC_19)"})
    out.attrs["vertical_structure"] = (f"the shear fitted on a B-spline basis (SPEC_19), epoch {problem.basis.epoch}; "
                                       "the reference wind held")
    return out


@dataclass
class Pass:
    """One linearization."""

    state: np.ndarray
    products: dict
    blocks: solver.Blocks
    radio_block: solver.Blocks
    wind: xr.Dataset


def linearize(problem: Problem, state, work: Path, pool, data=None) -> Pass:
    """The runs and the weights at `state`. `data`, for a synthetic case, replaces the observed `y`."""
    wind = wind_dataset(problem, state)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    lists = [copy_run(t.run, f"{problem.name}_{i}", wind, work, problem.root) for i, t in enumerate(problem.radio.targets)]
    paths = list(pool.map(_run_job, lists))
    products = {t.name: tree(path) for t, path in zip(problem.radio.targets, paths)}
    first = products[problem.radio.targets[0].name].to_dataset(inherit=False)
    labels = np.asarray(first["pressure_label_Pa"].values, dtype="float64")
    level_p = np.stack([layer_pressures(level.pressure_Pa) for level in problem.gradients.levels])
    all_p = np.concatenate([labels, level_p.ravel()])
    chunks = np.array_split(np.arange(all_p.size), WORKERS)
    jobs = [(wind, problem.composition, problem.gravity, problem.rotation, all_p[c], problem.lat.size) for c in chunks if c.size]
    unit = np.vstack(list(pool.map(_weight_job, jobs)))
    terms = np.vstack(list(pool.map(_term_job, [(wind, problem.gravity, problem.rotation, p) for p in level_p.ravel()])))
    shape = (*level_p.shape, problem.lat.size)
    radio_block = problem.radio.block(problem.basis, products, unit[:labels.size], problem.lat)
    gradient_block = problem.gradients.block(problem.basis, state, unit[labels.size:].reshape(shape), terms.reshape(shape), problem.lat)
    blocks = solver.concatenate([radio_block, gradient_block])
    if data is not None:
        blocks = solver.Blocks(data, blocks.h, blocks.H, blocks.R, blocks.groups)
    return Pass(state, products, blocks, radio_block, wind)


def _delivered(products, radio):
    return {t.name: np.asarray(products[t.name].to_dataset(inherit=False)["temperature_K"].values, dtype="float64")
            for t in radio.targets}


@dataclass
class Fit:
    state: np.ndarray
    solution: solver.Solution
    last: Pass                 #: the linearization the solution and its diagnostics stand on
    final: Pass                #: the runs at the final state
    changes_K: list
    passes: int
    sigma_wind: np.ndarray = None


def iterate(problem: Problem, start, work: Path, pool, data=None) -> Fit:
    """Solve, rebuild, rerun, relinearize, until the delivered temperatures settle. SPEC_19 section 1
    item 5. The first pass stands at `start`.

    Under the v0.3 budget plain passes swing from pass to pass and had not settled by pass 8, so a pass whose
    solve would still move the delivered temperature by the tolerance or more takes an Anderson(1) step, the
    solve mixed with the one before it. A fixed point of the acceleration is a fixed point of the plain pass:
    the route changes, the equations do not. An accelerated step can land near-still by chance, so it cannot end
    the iteration: only a plain step that changes the delivered temperature by under the tolerance does, and
    the fit returned is that plain solve."""
    state = np.asarray(start, dtype="float64")
    current = linearize(problem, state, work, pool, data)
    changes = []
    residual_before, solve_before = None, None
    for count in range(1, problem.max_passes + 1):
        blocks = current.blocks
        solution = solver.linearized(problem.prior, state, blocks)
        residual = solution.state - state
        moved = blocks.H @ residual
        expected = max(float(np.max(np.abs(moved[blocks.groups == t.name] * np.exp(blocks.y[blocks.groups == t.name]))))
                       for t in problem.radio.targets)
        plain = residual_before is None or expected < problem.tolerance_K
        if plain:
            following_state = solution.state
        else:
            step = residual - residual_before
            mix = float(residual @ step / (step @ step))
            following_state = solution.state - mix * (solution.state - solve_before)
        residual_before, solve_before = residual, solution.state
        following = linearize(problem, following_state, work, pool, data)
        now, before = _delivered(following.products, problem.radio), _delivered(current.products, problem.radio)
        change = {k: float(np.nanmax(np.abs(now[k] - before[k]))) for k in now}
        changes.append(change)
        print(f"pass {count} ({'plain' if plain else f'Anderson {mix:+.3f}'}): the solve's linear change in T "
              f"{expected:.3f}; largest change in delivered T {change}", flush=True)
        last = current
        state, current = following_state, following
        if plain and max(change.values()) < problem.tolerance_K:
            return Fit(state, solution, last, current, changes, count)
    raise RuntimeError(f"not converged in {problem.max_passes} passes; the changes {changes}")


def wind_sigma(problem: Problem, solution) -> np.ndarray:
    """The posterior standard deviation of the wind on its grid, `(latitude, pressure)`."""
    lat_part = problem.basis.lat_matrix(problem.lat)
    x_part = problem.basis.x_integral(problem.pressure)
    variance = np.empty((problem.lat.size, problem.pressure.size))
    for j in range(problem.pressure.size):
        matrix = problem.basis.row(lat_part, np.repeat(x_part[j:j + 1], problem.lat.size, axis=0))
        variance[:, j] = solver.posterior_variance(problem.prior, solution, matrix)
    return np.sqrt(np.maximum(variance, 0.0))


def write_wind(problem: Problem, fit: Fit) -> Path:
    fit.sigma_wind = wind_sigma(problem, fit.solution)
    out = wind_dataset(problem, fit.state, fit.sigma_wind)
    out.attrs["title"] = f"{problem.name}: the wind fitted to the radio targets and IRIS gradients (SPEC_19)"
    out.attrs["epoch_of_fit"] = problem.basis.epoch
    out.attrs["input_hashes"] = cio.input_hashes([problem.control_path], problem.output)
    cio.history_append(out, f"{TOOL}: {fit.passes} passes; largest last change in delivered T "
                            f"{max(fit.changes_K[-1].values()):.3f} K")
    return cio.write(problem.output, out, "wind", created_by=TOOL)


def run(control_path, pool=None) -> tuple:
    problem = read_control(control_path)
    own = pool is None
    pool = pool or ProcessPoolExecutor(max_workers=WORKERS)
    try:
        start = problem.basis.project(base_shear(problem), problem.lat, problem.pressure)
        fit = iterate(problem, start, problem.directory / "work", pool)
    finally:
        if own:
            pool.shutdown()
    write_wind(problem, fit)
    return problem, fit


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog=TOOL, description="Fit the wind to its observations; write kind W.")
    parser.add_argument("control", help="path to the fit's TOML control file")
    args = parser.parse_args(argv)
    problem, fit = run(args.control)
    print(f"wrote {problem.output}: {fit.passes} passes, last changes {fit.changes_K[-1]}")
    final = fit.final.blocks
    check = solver.consistency(final, solver.linearized(problem.prior, fit.state, final))
    print("misfit over expected, on the final linearization: " + ", ".join(f"{k} {v['ratio']:.2f}" for k, v in check.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
