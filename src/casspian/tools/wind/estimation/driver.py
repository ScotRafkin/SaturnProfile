"""`casspian-wind-fit`: the fit's control file, the relinearization, and the kind W file. SPEC_19, SPEC_21.

One pass: the wind from the state; the targets' runs (copies of their base runs with only the wind
changed, under `work/`); the weights; the solve, robust by Huber's rule. It repeats until the delivered
temperature changes by less than the tolerance at every target. Then each IRIS group's systematic part is sized
by the consistency check and the fit is made again once from there (SPEC_21 section 1 item 7), and the wind and its
posterior standard deviation are written.

Every amplitude, length, form and threshold is the control file's, each with its basis in a comment (SPEC_21
section 1 item 10); a source's own errors are in its data (item 5).
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
from casspian.lib import heat_capacity as hc
from casspian.lib import io as cio
from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.wind.estimation import solver
from casspian.tools.wind.estimation.gradient import Systematic, TemperatureGradients, layer_pressures, read_levels, wind_term
from casspian.tools.wind.estimation.prior import DensePrior, fine_grid, temperature_prior
from casspian.tools.wind.estimation.radio import Profile, RadioTargets, Representativeness, Target, interp_ln_p
from casspian.tools.wind.estimation.stability import GROUP as STABILITY, Stability
from casspian.tools.wind.estimation.state import Basis

TOOL = "casspian-wind-fit"
#: 66 levels at about 8 s each; the marches share nothing, so they run apart.
WORKERS = max(1, min(10, (os.cpu_count() or 2) - 2))
#: Every wind file the fit writes reaches this many decades above the base wind's top (1 Pa, the forward grid's), the
#: shear held at its top value there, the state's own rule above its highest knot: a column whose top moves a little
#: past 1 Pa must not stop the forward run (REVIEW_21_step1, check 5 item 1).
TOP_EXTENSION_DECADES = 1.0
#: `value_provenance` of the levels added above the base wind's top: extrapolated.
EXTRAPOLATED = 3


def extend_top(wind: xr.Dataset, decades=TOP_EXTENSION_DECADES) -> xr.Dataset:
    """The base wind with levels added above its top at its own top spacing in `ln p`, to `decades` above it. The
    values there are the top level's, to be overwritten from the state; their provenance is extrapolated."""
    p = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    if not np.all(np.diff(p) > 0):
        raise ValueError("the base wind's pressure must increase")
    step = np.log(p[1] / p[0])
    n = int(np.ceil(decades * np.log(10.0) / step - 1e-9))
    added = p[0] * np.exp(-step * np.arange(n, 0, -1))
    top = wind.isel(pressure=np.zeros(n, dtype=int))
    top = top.assign_coords(pressure_Pa=("pressure", added))
    top["value_provenance"].values[...] = EXTRAPOLATED
    out = xr.concat([top, wind], dim="pressure", data_vars="minimal", coords="minimal", compat="override",
                    combine_attrs="override")
    out.attrs["coverage_pressure_Pa"] = np.array([float(added[0]), float(p[-1])])
    return out


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
    prior: DensePrior
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
    huber_threshold: float
    systematic_length_deg: float
    systematic_form: str
    stability: Stability = None
    lat: np.ndarray = field(init=False)
    pressure: np.ndarray = field(init=False)
    wind_out: xr.Dataset = field(init=False)          #: the base wind extended above its top: the files' grid
    pressure_out: np.ndarray = field(init=False)

    def __post_init__(self):
        self.lat = np.asarray(self.wind["latitude_planetocentric_deg"].values, dtype="float64")
        self.pressure = np.asarray(self.wind["pressure_Pa"].values, dtype="float64")
        self.wind_out = extend_top(self.wind)
        self.pressure_out = np.asarray(self.wind_out["pressure_Pa"].values, dtype="float64")


def amplitude(table, prefix):
    """`a (p_0 / p)^e`, capped, from the control file's four keys under `prefix`."""
    a, p0, e, cap = (float(table[f"{prefix}_K"]), float(table[f"{prefix}_reference_Pa"]),
                     float(table[f"{prefix}_exponent"]), float(table[f"{prefix}_cap_K"]))
    return lambda p: np.minimum(a * (p0 / np.asarray(p, dtype="float64")) ** e, cap)


def mean_molar_mass(composition, pressure_Pa, latitude_deg) -> np.ndarray:
    """The kind C mixture's mean molar mass (kg/mol) at `pressure_Pa` and one latitude."""
    root = composition.to_dataset(inherit=False)
    species = composition["species"].to_dataset(inherit=False)
    names = [str(n) for n in species["species_name"].values]
    masses = np.asarray(species["molar_mass_kg_mol"].values, dtype="float64")
    levels = np.asarray(root["pressure_Pa"].values, dtype="float64")
    lat = np.asarray(root["latitude_planetocentric_deg"].values, dtype="float64")
    rising, north = np.argsort(np.log(levels)), np.argsort(lat)
    p = np.atleast_1d(np.asarray(pressure_Pa, dtype="float64"))
    total = np.zeros(p.shape)
    for name, mass in zip(names, masses):
        x = np.asarray(root[f"x_{name}"].values, dtype="float64")[np.ix_(rising, north)]
        column = np.array([np.interp(latitude_deg, lat[north], row) for row in x])
        total = total + np.interp(np.log(p), np.log(levels[rising]), column) * mass
    return total


def gravity_wave_speed(thermo_path, composition, latitude_deg, x, x_length) -> np.ndarray:
    """`c = N h` on the `ln p` values `x` (SPEC_21 v0.4 section 1 item 1a): `N` from the anchor's own profile, `h` the
    prior's vertical length in height, `x_length` scale heights; `a` in `L_R` is the gravity's normalization radius, the
    equatorial radius as the geoid takes it. With `H = R_s T / g`,
    `N^2 = (g^2 / R_s T)(kappa - d ln T / d ln p)`, so `N H = sqrt(R_s T (kappa - d ln T / d ln p))` and
    `c = x_length N H`, which needs no gravity. The slope is taken across the vertical length itself, centered,
    and inside the profile where it ends; `kappa = R / c_p` with `c_p` for normal hydrogen from the composition (SPEC_20), as the
    stability constraint takes it. Zero where the profile is not stable."""
    with xr.open_dataset(thermo_path) as th:
        p = np.asarray(th["pressure_Pa"].values, dtype="float64")
        T = np.asarray(th["temperature_K"].values, dtype="float64")
    good = np.isfinite(T) & np.isfinite(p)
    order = np.argsort(p[good])
    lp, lT = np.log(p[good][order]), np.log(T[good][order])
    half = 0.5 * float(x_length)
    low, high = np.maximum(x - half, lp[0]), np.minimum(x + half, lp[-1])
    slope = (np.interp(high, lp, lT) - np.interp(low, lp, lT)) / (high - low)
    T_x = np.exp(np.interp(x, lp, lT))
    p_x = np.exp(x)
    kappa = MOLAR_GAS_CONSTANT / hc.mixture_cp(composition, p_x, latitude_deg, T_x, "normal")
    R_s = MOLAR_GAS_CONSTANT / mean_molar_mass(composition, p_x, latitude_deg)
    return float(x_length) * np.sqrt(np.maximum(R_s * T_x * (kappa - slope), 0.0))


def build_prior(table, basis, wind, composition, gravity, rotation, root, anchor_latitude_deg) -> DensePrior:
    """The shear prior stated in temperature (SPEC_21 section 1 item 1): `w` held at the base run's wind and the
    anchor's `T_bar(p)`, uniform in latitude (REVIEW_21_step0 ruling 3); the latitude length floored at the
    deformation radius, with `c` from the same profile (v0.4 section 1 item 1a)."""
    lat_f, x_f = fine_grid(basis)
    p_f = np.exp(x_f)
    with xr.open_dataset(root / table["temperature_profile"]) as th:
        T_bar = interp_ln_p(np.asarray(th["pressure_Pa"].values), np.asarray(th["temperature_K"].values), p_f)
    edge = np.isfinite(T_bar)
    T_bar = np.interp(x_f, x_f[edge], T_bar[edge])
    lat_w = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    # balance_weight with T = 1 on each level, spread over the workers as a pass spreads it; w goes as 1 / T.
    chunks = np.array_split(np.arange(p_f.size), WORKERS)
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        unit = np.vstack(list(pool.map(_weight_job, [(wind, composition, gravity, rotation, p_f[c], lat_w.size)
                                                     for c in chunks if c.size])))
    w = np.stack([np.interp(lat_f, lat_w, row) for row in unit / T_bar[:, None]])
    speed = None
    if table.get("deformation_floor", False):
        speed = gravity_wave_speed(root / table["temperature_profile"], composition, anchor_latitude_deg, x_f,
                                   table["ln_p_correlation"])
    return temperature_prior(basis, w, T_bar, amplitude(table, "sigma_T")(p_f), table["latitude_correlation_deg"],
                             table["latitude_form"], table["ln_p_correlation"], table["ln_p_form"], speed=speed,
                             radius_m=float(gravity["normalization_radius_m"]),
                             omega=float(rotation["angular_rate_rad_s"]))


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
    if float(c["reference_wind"]["spread_ms"]) != 0.0:
        raise ValueError(f"{control_path}: the reference wind is held until wind observations enter (SPEC_23)")
    s = c["state"]
    basis = Basis.build(c["fit"]["epoch"], s["latitude_step_deg"], s["ln_p_step"], s["bottom_Pa"], s["top_Pa"],
                        float(wind["reference_level_pressure_Pa"]))
    prior = build_prior(c["prior"], basis, wind, composition, gravity, rotation, root,
                        float(c["radio"]["reference"]["latitude_planetocentric_deg"]))
    r, rep = c["radio"], c["representativeness"]
    path = lambda entry: str(root / entry) if entry else ""
    reference = Profile(r["reference"]["name"], r["reference"]["latitude_planetocentric_deg"], path(r["reference"]["thermo"]),
                        path(r["reference"].get("uncertainty", "")))
    targets = [Target(t["name"], t["latitude_planetocentric_deg"], path(t["thermo"]), path(t.get("uncertainty", "")), t["run"])
               for t in r["target"]]
    representativeness = Representativeness(float(rep["sigma_K"]), float(rep["sigma_reference_Pa"]), float(rep["sigma_exponent"]),
                                            float(rep["sigma_cap_K"]), float(rep["correlation_ln_p"]), rep["correlation_form"])
    radio = RadioTargets(reference, targets, r["correlation_ln_p"], r["floor_K"], representativeness, r["correlation_form"],
                         float(r.get("floor_correlation_ln_p", 0.0)), r.get("floor_correlation_form", "exponential"),
                         bool(r.get("average_layers", False)))
    g = c["gradient"]
    table = np.genfromtxt(root / g["data"], delimiter=",", names=True)
    levels = {int(k): float(v) for k, v in g["levels_Pa"].items()}
    with open(root / g["properties"], "rb") as handle:
        ranges = {int(k): v for k, v in tomllib.load(handle).get("level_range_Pa", {}).items()}
    lat_c = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    lat_g = np.asarray(wind["latitude_planetographic_deg"].values, dtype="float64")
    gradients = TemperatureGradients(read_levels(table, levels, lat_c, lat_g, g["sample_deg"], g["gap_deg"], ranges))
    stability = None
    if "stability" in c:
        ap, aT, _ = radio.profiles[reference.name]
        stability = Stability(float(c["stability"]["sigma"]), float(c["stability"]["width"]), composition,
                              {t.name: float(t.latitude_deg) for t in targets}, (ap, aT), float(reference.latitude_deg),
                              np.asarray(wind["pressure_Pa"].values, dtype="float64"))
    return Problem(c["fit"]["name"], here, root, basis, prior, radio, gradients, wind, composition, gravity, rotation,
                   c["solve"]["tolerance_K"], c["solve"]["max_passes"], here / c["fit"]["output"], control_path,
                   float(c["robust"]["huber_threshold"]), float(g["systematic"]["latitude_correlation_deg"]),
                   g["systematic"]["form"], stability)


def base_shear(problem: Problem) -> np.ndarray:
    """The base wind's shear on its own grid, `(pressure, latitude)`: its `du/dln p`."""
    x = np.log(problem.pressure)
    order = np.argsort(x)
    u = np.asarray(problem.wind["u_total_ms"].values, dtype="float64")
    out = np.empty((x.size, problem.lat.size))
    out[order] = np.gradient(u[:, order], x[order], axis=1).T
    return out


def wind_dataset(problem: Problem, state, sigma=None) -> xr.Dataset:
    out = problem.wind_out.copy(deep=True)
    u_ref = np.asarray(out["u_reference_ms"].values, dtype="float64")
    u = u_ref[:, None] + problem.basis.wind_change(state, problem.lat, problem.pressure_out)
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
    active: str = "none"      #: where the stability constraint is active at this linearization
    shared: np.ndarray = None  #: the radio's shared anchor term, on the first rows of `blocks`


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
    centers = problem.gradients.placements()
    level_p = np.stack([layer_pressures(center) for center in centers])
    all_p = np.concatenate([labels, level_p.ravel()])
    chunks = np.array_split(np.arange(all_p.size), WORKERS)
    jobs = [(wind, problem.composition, problem.gravity, problem.rotation, all_p[c], problem.lat.size) for c in chunks if c.size]
    unit = np.vstack(list(pool.map(_weight_job, jobs)))
    terms = np.vstack(list(pool.map(_term_job, [(wind, problem.gravity, problem.rotation, p) for p in level_p.ravel()])))
    shape = (*level_p.shape, problem.lat.size)
    unit_levels, terms_levels = unit[labels.size:].reshape(shape), terms.reshape(shape)
    at = {center: (unit_levels[k], terms_levels[k], center) for k, center in enumerate(centers)}
    nodes = None if problem.stability is None else problem.pressure
    delivered = problem.radio.delivered(problem.basis, products, unit[:labels.size], problem.lat, nodes)
    radio_block = problem.radio.block(problem.basis, products, unit[:labels.size], problem.lat, delivered,
                                      problem.pressure)
    gradient_block = problem.gradients.block(problem.basis, state, at, problem.lat)
    parts = [radio_block, gradient_block]
    active = "none"
    if problem.stability is not None:
        parts.append(problem.stability.block(delivered, problem.basis.size))
        active = (problem.stability.report() + "; excess already in the anchor, not penalized: "
                  + problem.stability.report_source())
    blocks = solver.concatenate(parts)
    if data is not None:
        # A synthetic case replaces the observations, `data` one value per observation row; the constraint's rows
        # keep their own `y`.
        observed = blocks.groups != STABILITY
        y = blocks.y.copy()
        y[observed] = data
        blocks = solver.Blocks(y, blocks.h, blocks.H, blocks.R, blocks.groups)
    return Pass(state, products, blocks, radio_block, wind, active, problem.radio.shared)


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
    weights: np.ndarray = None   #: Huber's weights on the last linearization, 1 where untouched
    weighted: solver.Blocks = None
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
        solution, weighted, weights = solver.robust(problem.prior, state, blocks, problem.huber_threshold,
                                                    exempt=(STABILITY,))
        residual = solution.state - state
        moved = blocks.H @ residual
        parts = solver.cost_parts(problem.prior, state, weighted, solution, current.shared)
        print(f"pass {count} solve, the cost's parts before -> after: " + "; ".join(
            f"{k} {parts['before'][k]:.4g} -> {parts['after'][k]:.4g}" for k in parts["before"])
            + f"; total {sum(parts['before'].values()):.5g} -> {sum(parts['after'].values()):.5g}", flush=True)
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
              f"{expected:.3f}; largest change in delivered T {change}; stability active at {current.active}", flush=True)
        last = current
        state, current = following_state, following
        if plain and max(change.values()) < problem.tolerance_K:
            return Fit(state, solution, last, current, changes, count, weights, weighted)
    raise RuntimeError(f"not converged in {problem.max_passes} passes; the changes {changes}")


def wind_sigma(problem: Problem, solution) -> np.ndarray:
    """The posterior standard deviation of the wind on the files' grid, `(latitude, pressure)`."""
    lat_part = problem.basis.lat_matrix(problem.lat)
    x_part = problem.basis.x_integral(problem.pressure_out)
    variance = np.empty((problem.lat.size, problem.pressure_out.size))
    for j in range(problem.pressure_out.size):
        matrix = problem.basis.row(lat_part, np.repeat(x_part[j:j + 1], problem.lat.size, axis=0))
        variance[:, j] = solver.posterior_variance(problem.prior, solution, matrix)
    return np.sqrt(np.maximum(variance, 0.0))


def size_systematic(problem: Problem, fit: Fit) -> dict:
    """Each IRIS group's systematic part, sized on the fit's last linearization to bring the group's weighted
    misfit to its expected value (SPEC_21 section 1 item 7). Zero for a group at or below it."""
    blocks = fit.weighted
    innovation = fit.solution.innovation
    sizes = {}
    for level in problem.gradients.levels:
        shape = problem.gradients.shape(level.name, problem.lat, problem.systematic_length_deg, problem.systematic_form)
        sizes[level.name] = solver.systematic_term(problem.prior, blocks, innovation, level.name, shape)
    return sizes


def write_wind(problem: Problem, fit: Fit, record: dict = None) -> Path:
    fit.sigma_wind = wind_sigma(problem, fit.solution)
    out = wind_dataset(problem, fit.state, fit.sigma_wind)
    out.attrs["title"] = f"{problem.name}: the wind fitted to the radio targets and IRIS gradients (SPEC_19, SPEC_21)"
    out.attrs["epoch_of_fit"] = problem.basis.epoch
    out.attrs["input_hashes"] = cio.input_hashes([problem.control_path], problem.output)
    for key, value in (record or {}).items():
        out.attrs[key] = value
    cio.history_append(out, f"{TOOL}: {fit.passes} passes; largest last change in delivered T "
                            f"{max(fit.changes_K[-1].values()):.3f} K")
    return cio.write(problem.output, out, "wind", created_by=TOOL)


def run(control_path, pool=None, start=None, write=True) -> tuple:
    """The fit, its IRIS systematic parts sized, and the fit again once from there. Returns the problem, the final
    fit and a record: the first fit, the sizes, and the consistency of each fit."""
    problem = read_control(control_path)
    own = pool is None
    pool = pool or ProcessPoolExecutor(max_workers=WORKERS)
    try:
        if start is None:
            start = problem.basis.project(base_shear(problem), problem.lat, problem.pressure)
        first = iterate(problem, start, problem.directory / "work", pool)
        before = solver.consistency(first.weighted, first.solution, exclude=(STABILITY,))
        sizes = size_systematic(problem, first)
        print("IRIS systematic parts, per radian: " + ", ".join(f"{k} {v:.4g}" for k, v in sizes.items()), flush=True)
        problem.gradients.systematic = {k: Systematic(v, problem.systematic_length_deg, problem.systematic_form)
                                        for k, v in sizes.items() if np.isfinite(v) and v > 0.0}
        fit = iterate(problem, first.state, problem.directory / "work", pool) if problem.gradients.systematic else first
        after = solver.consistency(fit.weighted, fit.solution, exclude=(STABILITY,))
    finally:
        if own:
            pool.shutdown()
    record = {"first": first, "sizes": sizes, "before": before, "after": after}
    if write:
        write_wind(problem, fit, {
            "iris_systematic_per_rad": ", ".join(f"{k} {v:.6g}" for k, v in sizes.items()),
            "iris_systematic_latitude_correlation": f"{problem.systematic_length_deg:g} deg, {problem.systematic_form}",
            "passes": f"first fit {first.passes}, second fit {fit.passes if fit is not first else 0}"})
    return problem, fit, record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog=TOOL, description="Fit the wind to its observations; write kind W.")
    parser.add_argument("control", help="path to the fit's TOML control file")
    args = parser.parse_args(argv)
    problem, fit, record = run(args.control)
    print(f"wrote {problem.output}: first fit {record['first'].passes} passes, second fit {fit.passes} passes, "
          f"last changes {fit.changes_K[-1]}")
    for name, check in (("before the systematic parts", record["before"]), ("after", record["after"])):
        print(f"misfit over expected, {name}: " + ", ".join(f"{k} {v['ratio']:.2f}" for k, v in check.items()))
    down = fit.weights < 1.0
    print(f"Huber's rule downweighted {int(down.sum())} of {down.size} rows; the smallest weight {float(fit.weights.min()):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
