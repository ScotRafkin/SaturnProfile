"""Acceptance for SPEC_19 v0.4 Step 1: the wind fit, and hypothesis 3 at the Voyager epoch.

The resolution test first; then `casspian-wind-fit` on `fits/voyager_fit_1/voyager_fit_1.toml`, the linear
ensemble on its final linearization, the three runs under its wind, and the tables and figures for the
author. Writes under `reports/step19_1/`, the wind file and the runs. The hypotheses 1 and 2 runs are read
as they are.

Latitudes are given planetographic, with planetocentric in parentheses, where they name a target.

Run from the repository root; it takes about four hours.
"""

import json
import shutil
import subprocess
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import cho_solve

from casspian.forward import transfer as forward_transfer
from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.lindal import radio_profiles as rp
from casspian.tools.run import new_run, run_inputs
from casspian.tools.wind.estimation import driver, solver
from casspian.tools.wind.estimation.gradient import LAYER_THICKNESS, LAYER_WEIGHTS, layer_pressures, wind_term
from casspian.tools.wind.estimation.radio import interp_ln_p, temperature_on_wind_grid

warnings.filterwarnings("ignore", message=".*input_hashes.*")
warnings.filterwarnings("ignore", message=".*solar_longitude.*")
HERE = Path("reports/step19_1")
CONTROL = Path("fits/voyager_fit_1/voyager_fit_1.toml")
TABLE_MBAR = (730.0, 290.0, 110.0, 30.0, 10.0, 4.5, 1.0)
#: Target, the three hypotheses' runs (hypothesis 3's are made here), its observed file.
RUNS = {"31.2 S": ("lindal_iris_v_sl_ii_110", "lindal_iris_v_sl_iii_110", "lindal_voyager_fit_1"),
        "74.2 S": ("lindal_iris_v_sl_ii_110_74s", "lindal_iris_v_sl_iii_110_74s", "lindal_voyager_fit_1_74s"),
        "3.0 S": ("lindal_iris_v_sl_ii_110_03s", "lindal_iris_v_sl_iii_110_03s", "lindal_voyager_fit_1_03s")}
#: Above this the shear is held; the grid the shear is compared on.
SHEAR_RANGE_Pa = (20.0, 1.3e5)
SEED = 19
#: SPEC_19 section 3 acceptance 1(b): the members.
ENSEMBLE = 20
#: Planetocentric, section 3 acceptance 2.
PROFILE_LATITUDES = (50.0, 20.0, 5.0, -15.0, -50.0)
#: The latitude correlations of the candidate systematic term, section 2 item 3: the term, and a sensitivity.
TERM_CORRELATIONS = (20.0, 40.0)
#: The first fit, under v0.2's budget, delivered minus observed at TABLE_MBAR (REPORT_19_step1 section 4).
FIRST_FIT = {"31.2 S": (-1.34, 1.56, 6.97, 0.12, 0.22, -0.43, -1.50),
             "74.2 S": (-1.35, -3.73, -1.68, 0.04, -0.24, -0.24, 1.02),
             "3.0 S": (-0.72, 2.44, 4.89, 1.58, -11.60, 4.09, 7.43)}
AT_110_Pa = 11000.0
#: Planetocentric: the anchor, then segments along the path to 2.4 S and on to 26.1 S (3.0 S and 31.2 S planetographic).
PATH_EDGES = (30.8056, 20.0, 10.0, 5.0, -2.3736, -5.0, -15.0, -26.0610)
results, summary = [], {}


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'REPORTED' if passed is None else 'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")
    sys.stdout.flush()


def shear_sigma(problem, solution, lat, pressure):
    """The posterior standard deviation of the shear, `(pressure, latitude)`."""
    a = problem.basis.lat_matrix(lat)
    b = problem.basis.x_matrix(pressure)
    out = np.empty((pressure.size, lat.size))
    for k in range(pressure.size):
        matrix = problem.basis.row(a, np.repeat(b[k:k + 1], lat.size, axis=0))
        out[k] = np.sqrt(np.maximum(solver.posterior_variance(problem.prior, solution, matrix), 0.0))
    return out


def make_run(name, base, wind_path):
    """A run from `base` with only the wind changed: `[shear]` the identity on the fit's wind."""
    directory = Path("forward") / name
    if directory.exists():
        shutil.rmtree(directory)
    new_run.new_run(name, Path("forward") / base)
    build = directory / f"{name}_build.toml"
    text = build.read_bytes().decode("utf-8")
    start = text.index("\n[shear]") + 1
    stop = text.index("\n[composition]") + 1
    text = (text[:start] + "[shear]\n# SPEC_19 Step 1: hypothesis 3, the fit's wind; the only change from " + base + ".\n"
            f'source = "../../{wind_path.as_posix()}"\noutput = "inputs/{name}_wind.nc"\ncase   = "identity"\n\n' + text[stop:])
    build.write_bytes(text.encode("utf-8"))
    run_inputs.build(build)
    return forward_transfer.run(directory / f"{name}.toml").product


def draw_from_prior(problem, rng):
    """A state drawn from the prior: its covariance is a Kronecker product of two nearly singular
    correlations, so each is square-rooted by its eigenvalues."""
    def root(c):
        values, vectors = np.linalg.eigh(c)
        return vectors * np.sqrt(np.maximum(values, 0.0))
    z = rng.standard_normal((problem.basis.n_lat, problem.basis.n_x))
    return problem.prior.mean + problem.prior.amplitude * (root(problem.prior.c_lat) @ z @ root(problem.prior.c_x).T).ravel()


def at_level(problem, products, blocks, values, at=AT_110_Pa):
    """Delivered minus observed (K) at a pressure at each radio target, from a `ln T` per row."""
    out = {}
    for target in problem.radio.targets:
        table = problem.radio.observed[target.name]
        p = np.asarray(products[target.name].to_dataset(inherit=False)["pressure_Pa"].values, dtype="float64")
        p = p[np.isfinite(interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], p))]
        T = np.exp(interp_ln_p(p, values[blocks.groups == target.name], np.array([at]))[0])
        out[target.name] = float(T - interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], np.array([at]))[0])
    return out


def landing(problem, products, blocks, state, prior, R=None):
    """Where 110 mbar lands, with all the data and with IRIS alone, and the group ratios: one solve on this
    linearization, not a fit. `R` replaces the stated errors."""
    R = blocks.R if R is None else R
    innovation = blocks.y - blocks.h + blocks.H @ (state - prior.mean)
    solved = solver.solve(prior, blocks.H, R, innovation)
    iris = np.array([g.startswith("IRIS") for g in blocks.groups])
    alone = solver.solve(prior, blocks.H[iris], R[np.ix_(iris, iris)], innovation[iris])
    check = solver.consistency(solver.Blocks(blocks.y, blocks.h, blocks.H, R, blocks.groups), solved)
    return {"all_data": at_level(problem, products, blocks, blocks.h + blocks.H @ (solved.state - state)),
            "iris_only": at_level(problem, products, blocks, blocks.h + blocks.H @ (alone.state - state)),
            "ratios": {k: v["ratio"] for k, v in check.items()}}


def size_terms(problem, products, blocks, state, length):
    """The candidate systematic term of section 2 item 3, an error correlated in latitude: its size per IRIS
    level, one group at a time, and what all of them together do; one solve each. If no size brings a level
    to its expected value, a long correlation cannot absorb the misfit: a scan of sizes is given instead."""
    prior = problem.prior
    innovation = blocks.y - blocks.h + blocks.H @ (state - prior.mean)
    shapes, stated, per_deg = {}, {}, {}
    for level in problem.gradients.levels:
        la = problem.lat[level.indices]
        shapes[level.name] = np.exp(-0.5 * ((la[:, None] - la[None, :]) / length) ** 2)
        stated[level.name] = float(np.median(level.sigma))
        per_deg[level.name] = float(np.median(level.temperature_K)) * np.pi / 180.0

    def applied(sizes):
        R = blocks.R.copy()
        for name, size in sizes.items():
            g = np.flatnonzero(blocks.groups == name)
            R[np.ix_(g, g)] += size ** 2 * shapes[name]
        return landing(problem, products, blocks, state, prior, R)

    sizes = {name: solver.systematic_term(prior, blocks, innovation, name, shape) for name, shape in shapes.items()}
    out = {"correlation_deg": length, "sizes_K_per_deg": {k: v * per_deg[k] for k, v in sizes.items()}}
    finite = {k: v for k, v in sizes.items() if np.isfinite(v)}
    out["applied"] = applied(finite)
    out["scan"] = {f"{factor:g} x the median stated error": {"sizes_K_per_deg": {k: factor * stated[k] * per_deg[k] for k in shapes},
                                                             **applied({k: factor * stated[k] for k in shapes})} for factor in (1.0, 2.0, 4.0)}
    return out


def path_segments(problem, fit, start, hyp1, at=AT_110_Pa):
    """The shear differences integrated along the path from the anchor at `at`, K per segment, with the
    final weight; and the three shears on the path."""
    x = np.log(problem.pressure)
    order = np.argsort(x)
    s1 = np.array([np.interp(np.log(at), x[order], hyp1[order, j]) for j in range(problem.lat.size)])
    s0 = problem.basis.shear(start, problem.lat, np.array([at]))[0]
    s3 = problem.basis.shear(fit.state, problem.lat, np.array([at]))[0]
    unit = lw.balance_weight(fit.final.wind, problem.composition, problem.gravity, problem.rotation, np.array([at]), np.ones((1, problem.lat.size)))[0]
    root = fit.final.products["3.0 S"]
    labels, T_wind = temperature_on_wind_grid(root.to_dataset(inherit=False), root["isobars"].to_dataset(inherit=False), problem.lat)
    row = T_wind[np.argmin(np.abs(np.log(labels / at)))]
    T = np.interp(problem.lat, problem.lat[np.isfinite(row)], row[np.isfinite(row)])
    out = {}
    for name, ds in (("fit minus hypothesis 1", s3 - s1), ("projection minus hypothesis 1", s0 - s1)):
        parts = []
        for a, b in zip(PATH_EDGES[:-1], PATH_EDGES[1:]):
            nodes = np.concatenate([[b], problem.lat[(problem.lat > b) & (problem.lat < a)], [a]])
            parts.append(-float(np.trapezoid(np.interp(nodes, problem.lat, unit / T * ds), np.radians(nodes)) * np.interp(nodes.mean(), problem.lat, T)))
        out[name] = parts
    return out, {"lat": problem.lat, "hypothesis 1": s1, "projection": s0, "fit": s3}


def radio_modes(problem, blocks, solution):
    """Each radio group's misfit in the eigenmodes of its `R`, against expected: which scales carry it."""
    r = blocks.R @ cho_solve(solution.factor, solution.innovation)
    expected_cov = blocks.R @ cho_solve(solution.factor, blocks.R)
    out = {}
    for target in problem.radio.targets:
        g = blocks.groups == target.name
        values, vectors = np.linalg.eigh(blocks.R[np.ix_(g, g)])
        values, vectors = values[::-1], vectors[:, ::-1]
        misfit = (vectors.T @ r[g]) ** 2 / values
        expected = np.einsum("ij,ij->j", vectors, expected_cov[np.ix_(g, g)] @ vectors) / values
        out[target.name] = [(lo + 1, hi, float(misfit[lo:hi].sum()), float(expected[lo:hi].sum()))
                            for lo, hi in ((0, 3), (3, 8), (8, 20), (20, int(g.sum())))]
    return out


def synthetic(problem, pool, truth, label, rng, work):
    """SPEC_19 section 3 acceptance 1(a). The truth is `truth`, the data its model equivalents with noise
    drawn from section 2's errors, and the fit starts from the prior mean."""
    at_truth = driver.linearize(problem, truth, work, pool)
    data = at_truth.blocks.h + np.linalg.cholesky(at_truth.blocks.R) @ rng.standard_normal(at_truth.blocks.h.size)
    fit = driver.iterate(problem, problem.prior.mean, work, pool, data)
    solution = solver.linearized(problem.prior, fit.state, fit.final.blocks)
    levels = np.geomspace(*SHEAR_RANGE_Pa, 40)
    lat = np.arange(-85.0, 85.1, 2.5)
    fitted, true = (problem.basis.shear(x, lat, levels) for x in (fit.state, truth))
    sigma = shear_sigma(problem, solution, lat, levels)
    prior_sd = float(np.sqrt(np.mean(np.einsum(
        "mi,mi->m", *(lambda m: (problem.prior.apply(m), m))(problem.basis.row(problem.basis.lat_matrix(lat), np.repeat(problem.basis.x_matrix([1000.0]), lat.size, axis=0)))))))
    z = (fitted - true) / sigma
    check = solver.consistency(fit.final.blocks, solution)
    constrained = sigma < 0.5 * prior_sd
    out = {"passes": fit.passes, "changes": fit.changes_K, "rms_error": float(np.sqrt(np.mean((fitted - true) ** 2))),
           "rms_true": float(np.sqrt(np.mean(true ** 2))), "within_1": float(np.mean(np.abs(z) <= 1.0)),
           "within_2": float(np.mean(np.abs(z) <= 2.0)), "within_2_constrained": float(np.mean(np.abs(z[constrained]) <= 2.0)),
           "constrained_fraction": float(np.mean(constrained)), "prior_sd": prior_sd,
           "rms_error_constrained": float(np.sqrt(np.mean((fitted - true)[constrained] ** 2))),
           "total_ratio": check["all"]["ratio"], "groups": {k: v["ratio"] for k, v in check.items()},
           "count_shear": int(z.size), "count_within_2": int(np.sum(np.abs(z) <= 2.0)),
           "count_within_1": int(np.sum(np.abs(z) <= 1.0))}
    picks = (15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 0.0, -2.5, -5.0, -7.5, -10.0, -15.0)
    at = np.array([AT_110_Pa])
    out["jet_at_110_mbar"] = {"latitudes": picks,
                              "truth": np.interp(picks, problem.lat, problem.basis.shear(truth, problem.lat, at)[0]).tolist(),
                              "recovered": np.interp(picks, problem.lat, problem.basis.shear(fit.state, problem.lat, at)[0]).tolist()}
    shutil.rmtree(work, ignore_errors=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    for ax, field, title, limit in ((axes[0], true, f"truth ({label})", 40), (axes[1], fitted, "recovered", 40),
                                    (axes[2], z, "(recovered - truth) / posterior sd", 3)):
        mesh = ax.pcolormesh(lat, levels / 100, field, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="nearest")
        fig.colorbar(mesh, ax=ax)
        ax.set_yscale("log")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("planetocentric latitude (deg)")
    axes[0].invert_yaxis()
    axes[0].set_ylabel("pressure (mbar)")
    fig.suptitle(f"The synthetic test, {label}: shear, m/s per unit ln p")
    fig.tight_layout()
    fig.savefig(HERE / f"synthetic_{label}.png", dpi=110)
    plt.close(fig)
    return out


def linear_ensemble(problem, state, blocks):
    """Section 3 acceptance 1(b), SPEC_19 v0.4: truths drawn from the prior, their data made with the operator
    at the final linearization, `y = h + H (truth - state) + noise`, and one solve each. No truth is marched:
    those of this prior go unphysical at 74.2 S (tests/step19_1/survey_truths.py)."""
    levels = np.geomspace(*SHEAR_RANGE_Pa, 40)
    lat = np.arange(-85.0, 85.1, 2.5)
    root = np.linalg.cholesky(blocks.R)
    # The posterior covariance does not depend on the data, so one sd serves every member.
    sigma = shear_sigma(problem, solver.linearized(problem.prior, state, blocks), lat, levels)
    members = {}
    for i in range(ENSEMBLE):
        rng = np.random.default_rng([SEED, i])
        truth = draw_from_prior(problem, rng)
        data = blocks.h + blocks.H @ (truth - state) + root @ rng.standard_normal(blocks.h.size)
        drawn = solver.Blocks(data, blocks.h, blocks.H, blocks.R, blocks.groups)
        solution = solver.linearized(problem.prior, state, drawn)
        z = (problem.basis.shear(solution.state, lat, levels) - problem.basis.shear(truth, lat, levels)) / sigma
        members[i] = {"count_shear": int(z.size), "count_within_1": int(np.sum(np.abs(z) <= 1.0)), "count_within_2": int(np.sum(np.abs(z) <= 2.0)),
                      "within_2": float(np.mean(np.abs(z) <= 2.0)),
                      "groups": {k: v["ratio"] for k, v in solver.consistency(drawn, solution).items()}}
    (HERE / "ensemble.json").write_text(json.dumps(members, indent=1, default=float), encoding="utf-8")
    return members


def report_linear_ensemble(problem, state, blocks):
    members = linear_ensemble(problem, state, blocks)
    total = sum(m["count_shear"] for m in members.values())
    pooled_2 = sum(m["count_within_2"] for m in members.values()) / total
    pooled_1 = sum(m["count_within_1"] for m in members.values()) / total
    mean_ratio = {k: float(np.mean([m["groups"][k] for m in members.values()])) for k in members[0]["groups"]}
    spread_ratio = {k: float(np.std([m["groups"][k] for m in members.values()], ddof=1)) for k in mean_ratio}
    per_group = {k: k == "all" or 0.8 <= v <= 1.2 for k, v in mean_ratio.items()}
    prior_T = {t.name: float(np.max(np.sqrt(np.einsum("mi,mi->m", problem.prior.apply(blocks.H[blocks.groups == t.name]), blocks.H[blocks.groups == t.name])) * np.exp(blocks.h[blocks.groups == t.name])))
               for t in problem.radio.targets}
    summary["synthetic_ensemble"] = {"members": len(members), "pooled_within_2": pooled_2, "pooled_within_1": pooled_1, "mean_ratio": mean_ratio,
                                     "spread_ratio": spread_ratio, "prior_largest_sd_in_T_K": prior_T}
    record("1b", f"the linear ensemble, {ENSEMBLE} truths drawn from the prior, their data made with the operator at the final linearization: "
                 "pooled at least 90% of the shear within 2 posterior sd, each group's mean ratio 0.8 to 1.2",
           pooled_2 >= 0.90 and all(per_group.values()),
           f"pooled within 1 posterior sd {pooled_1:.1%}, within 2 {pooled_2:.1%}\n"
           "each group's mean ratio (member spread): " + ", ".join(f"{k} {v:.2f} ({spread_ratio[k]:.2f})" for k, v in mean_ratio.items()) + "\n"
           "within 2 sd per member: " + " ".join(f"{m['within_2']:.0%}" for m in members.values()) + "\n"
           "the prior's largest sd in delivered T on this linearization (K): " + ", ".join(f"{k} {v:.0f}" for k, v in prior_T.items()) + "\n"
           f"members {HERE / 'ensemble.json'}")


def report_one_solve(problem, at_truth):
    """Item 1(a)'s draw of noise in one linear solve at its truth, beside the iterated fit's ratios: which of those
    the draw explains and which the iteration adds."""
    b = at_truth.blocks
    data = b.h + np.linalg.cholesky(b.R) @ np.random.default_rng(SEED).standard_normal(b.h.size)
    drawn = solver.Blocks(data, b.h, b.H, b.R, b.groups)
    solution = solver.linearized(problem.prior, at_truth.state, drawn)
    ratios = {k: v["ratio"] for k, v in solver.consistency(drawn, solution).items()}
    picks = np.array((15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 0.0, -2.5, -5.0, -7.5, -10.0, -15.0))
    shear = problem.basis.shear(solution.state, picks, np.array([AT_110_Pa]))[0]
    summary["synthetic_resolution_one_solve"] = {"groups": ratios, "jet_at_110_mbar": shear.tolist()}
    record("1a, continued", "the same noise in one linear solve at the truth, no iteration; no bound", None,
           "misfit over expected: " + ", ".join(f"{k} {v:.2f}" for k, v in ratios.items()) + "\n"
           "the jet at 110 mbar, planetocentric latitude: " + " ".join(f"{v:g}" for v in picks) + "\n"
           "  one solve " + " ".join(f"{v:+6.1f}" for v in shear))


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    pool = ProcessPoolExecutor(max_workers=driver.WORKERS)
    problem = driver.read_control(CONTROL)

    # ---- 1. the synthetic tests ------------------------------------------------------------------
    truth = problem.basis.project(driver.base_shear(problem), problem.lat, problem.pressure)
    syn = synthetic(problem, pool, truth, "hypothesis_1", np.random.default_rng(SEED), problem.directory / "work_synthetic")
    # Check 1(b) stands on the final linearization, so it follows the fit.
    summary["synthetic_resolution"] = syn
    jet = syn["jet_at_110_mbar"]
    record("1a", "the resolution test, hypothesis 1's shear on the basis as the truth; no bound", None,
           f"{syn['passes']} passes from the prior mean, the changes {syn['changes']}\n"
           f"RMS shear error {syn['rms_error']:.2f} against a true RMS of {syn['rms_true']:.2f} m/s per unit ln p; the prior's point sd {syn['prior_sd']:.1f}\n"
           f"within 1 posterior sd {syn['within_1']:.1%}, within 2 {syn['within_2']:.1%}\n"
           f"where the data constrain it (posterior sd under half the prior's, {syn['constrained_fraction']:.0%} of the grid): "
           f"within 2 sd {syn['within_2_constrained']:.1%}, RMS error {syn['rms_error_constrained']:.2f}\n"
           "misfit over expected: " + ", ".join(f"{k} {v:.2f}" for k, v in syn["groups"].items()) + "\n"
           "the jet at 110 mbar, planetocentric latitude: " + " ".join(f"{v:g}" for v in jet["latitudes"]) + "\n"
           "  truth     " + " ".join(f"{v:+6.1f}" for v in jet["truth"]) + "\n"
           "  recovered " + " ".join(f"{v:+6.1f}" for v in jet["recovered"]) + "\n"
           f"figure {HERE / 'synthetic_hypothesis_1.png'}")

    # ---- the fit --------------------------------------------------------------------------------
    problem, fit = driver.run(CONTROL, pool)
    wind_path = problem.output.relative_to(Path.cwd())
    from_run = driver._read(problem.output, "wind")["u_total_ms"].values.copy()
    subprocess.run(["casspian-wind-fit", str(CONTROL)], check=True)
    from_command = driver._read(problem.output, "wind")["u_total_ms"].values
    command_difference = float(np.max(np.abs(from_command - from_run)))
    for target, runs in RUNS.items():
        make_run(runs[2], runs[0], wind_path)
    blocks = fit.final.blocks
    solution = solver.linearized(problem.prior, fit.state, blocks)
    groups = list(dict.fromkeys(blocks.groups.tolist()))

    report_linear_ensemble(problem, fit.state, blocks)

    # Hypothesis 1 on the basis, the prior not involved: what the basis itself can hold.
    start = truth
    projected = driver.linearize(problem, start, problem.directory / "work_projection", pool)
    shutil.rmtree(problem.directory / "work_projection", ignore_errors=True)
    report_one_solve(problem, projected)

    # ---- 2. hypothesis 3, for the author --------------------------------------------------------
    products = {name: driver.tree(Path("forward") / name / "output" / f"{name}_profile.nc") for runs in RUNS.values() for name in runs}
    roots = {k: v.to_dataset(inherit=False) for k, v in products.items()}
    lines, delivered = [f"delivered minus observed (K) at {', '.join(f'{p:g}' for p in TABLE_MBAR)} mbar; sigma_T the observed error:"], {}
    observed = {t.name: rp.read_csv(t.observed) for t in problem.radio.targets}
    for target, runs in RUNS.items():
        table = observed[target]
        at = np.array(TABLE_MBAR) * 100.0
        obs = interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], at)
        delivered[target] = {}
        lines.append(f"{target}:")
        for label, name in zip(("hypothesis 1", "hypothesis 2", "hypothesis 3"), runs):
            T = interp_ln_p(roots[name]["pressure_Pa"].values, roots[name]["temperature_K"].values, at)
            delivered[target][label] = (T - obs).tolist()
            lines.append(f"  {label:13s}" + " ".join(f"{v:+7.2f}" for v in T - obs))
        lines.append(f"  {'first fit':13s}" + " ".join(f"{v:+7.2f}" for v in FIRST_FIT[target]))
        lines.append(f"  {'sigma_T':13s}" + " ".join(f"{v:7.2f}" for v in problem.radio.sigma(at)))
    summary["delivered_minus_observed"] = delivered

    levels = np.geomspace(0.3 * 100, 1.3e5, 60)
    lat = problem.lat
    hyp1 = driver.base_shear(problem)
    wind2 = driver._read(Path("forward") / RUNS["31.2 S"][1] / "inputs" / f"{RUNS['31.2 S'][1]}_wind.nc", "wind")
    x = np.log(problem.pressure)
    order = np.argsort(x)
    hyp2 = np.empty(hyp1.shape)
    hyp2[order] = np.gradient(np.asarray(wind2["u_total_ms"].values)[:, order], x[order], axis=1).T
    on_levels = lambda s: np.stack([np.interp(np.log(levels), x[order], s[order, j]) for j in range(lat.size)], axis=1)
    s3 = problem.basis.shear(fit.state, lat, levels)
    sd3 = shear_sigma(problem, solution, lat, levels)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True, sharey=True)
    for ax, field, title, cmap, limit in ((axes[0, 0], on_levels(hyp1), "hypothesis 1", "RdBu_r", 40), (axes[0, 1], on_levels(hyp2), "hypothesis 2", "RdBu_r", 40),
                                          (axes[1, 0], s3, "hypothesis 3, voyager_fit_1", "RdBu_r", 40), (axes[1, 1], sd3, "hypothesis 3, posterior standard deviation", "viridis", None)):
        mesh = ax.pcolormesh(lat, levels / 100, field, cmap=cmap, vmin=-limit if limit else 0, vmax=limit if limit else 20, shading="nearest")
        fig.colorbar(mesh, ax=ax, label="m/s per unit ln p")
        ax.set_yscale("log")
        ax.set_title(title, fontsize=10)
        for target in problem.radio.targets:
            ax.axvline(target.latitude_deg, color="0.3", linewidth=0.6, linestyle=":")
    axes[0, 0].invert_yaxis()
    for ax in axes[:, 0]:
        ax.set_ylabel("pressure (mbar)")
    for ax in axes[1]:
        ax.set_xlabel("planetocentric latitude (deg)")
    fig.tight_layout()
    fig.savefig(HERE / "shear.png", dpi=110)
    plt.close(fig)

    # The shear and its posterior sd against pressure, through IRIS's top level.
    span = np.geomspace(SHEAR_RANGE_Pa[0], SHEAR_RANGE_Pa[1], 160)
    lat_g = np.asarray(problem.wind["latitude_planetographic_deg"].values, dtype="float64")
    hyp1_span = np.stack([np.interp(np.log(span), x[order], hyp1[order, j]) for j in range(lat.size)], axis=1)
    named = lambda c: f"{abs(np.interp(c, lat, lat_g)):.1f} {'N' if c > 0 else 'S'} ({abs(c):g} {'N' if c > 0 else 'S'} planetocentric)"
    fig, axes = plt.subplots(2, len(PROFILE_LATITUDES), figsize=(17, 9), sharey=True, sharex="row")
    summary["profiles"] = {}
    table_p = np.array([1000.0, 730.0, 500.0, 290.0, 200.0, 150.0, 110.0, 90.0, 70.0, 50.0, 30.0, 10.0, 1.0]) * 100.0
    lines.append("the posterior sd of the shear (m/s per unit ln p) against pressure, mbar: " + " ".join(f"{p / 100:g}" for p in table_p))
    for k, c in enumerate(PROFILE_LATITUDES):
        s_fit = problem.basis.shear(fit.state, np.array([c]), span)[:, 0]
        sd_fit = shear_sigma(problem, solution, np.array([c]), span)[:, 0]
        s_h1 = np.array([np.interp(c, lat, row) for row in hyp1_span])
        top, bottom = axes[0, k], axes[1, k]
        for level in problem.gradients.levels:
            for ax in (top, bottom):
                ax.axhspan(level.pressure_Pa * np.exp(-0.5 * LAYER_THICKNESS) / 100, level.pressure_Pa * np.exp(0.5 * LAYER_THICKNESS) / 100, color="0.92", linewidth=0)
        top.fill_betweenx(span / 100, s_fit - sd_fit, s_fit + sd_fit, color="C0", alpha=0.25, linewidth=0, label="fit, 1 sd")
        top.plot(s_fit, span / 100, color="C0", linewidth=1.5, label="hypothesis 3")
        top.plot(s_h1, span / 100, color="k", linewidth=1.2, linestyle="--", label="hypothesis 1")
        top.axvline(0, color="0.6", linewidth=0.6)
        top.set_title(named(c), fontsize=10)
        top.set_xlabel("shear (m/s per unit ln p)")
        bottom.plot(sd_fit, span / 100, color="C0", linewidth=1.5)
        bottom.set_xlabel("posterior sd (m/s per unit ln p)")
        bottom.set_xscale("log")
        at_table = np.interp(np.log(table_p), np.log(span), sd_fit)
        summary["profiles"][str(c)] = {"sd_at": dict(zip((table_p / 100).tolist(), at_table.tolist()))}
        lines.append(f"  {named(c)}: " + " ".join(f"{v:.1f}" for v in at_table))
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylim(SHEAR_RANGE_Pa[1] / 100, SHEAR_RANGE_Pa[0] / 100)
    axes[0, 0].set_ylabel("pressure (mbar); gray, IRIS's layers")
    axes[1, 0].set_ylabel("pressure (mbar); gray, IRIS's layers")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "profiles.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 7), sharey=True)
    for ax, (target, runs) in zip(axes, RUNS.items()):
        table = observed[target]
        for name, label, color in zip(runs, ("hypothesis 1", "hypothesis 2", "hypothesis 3"), ("C2", "C4", "C3")):
            p = roots[name]["pressure_Pa"].values
            d = roots[name]["temperature_K"].values - interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], p)
            ax.plot(d, p / 100, color=color, label=label, linewidth=1.3)
        p = np.geomspace(30.0, 1.3e5, 100)
        s = problem.radio.sigma(p)
        ax.fill_betweenx(p / 100, -s, s, color="0.85", label="sigma_T")
        ax.axvline(0, color="0.5", linewidth=0.7)
        ax.set_yscale("log")
        ax.set_xlim(-25, 25)
        ax.set_title(target)
        ax.set_xlabel("delivered minus observed (K)")
        ax.legend(fontsize=8)
    axes[0].set_ylim(1300, 0.3)
    axes[0].set_ylabel("pressure (mbar)")
    fig.tight_layout()
    fig.savefig(HERE / "targets.png", dpi=110)
    plt.close(fig)

    # IRIS residuals, K per degree.
    lines.append("IRIS residuals at the fit, data minus the model's equivalent, K per planetocentric degree:")
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    summary["iris"] = {}
    for ax, level in zip(axes, problem.gradients.levels):
        g = blocks.groups == level.name
        scale = level.temperature_K * np.pi / 180.0
        r = (blocks.y - blocks.h)[g] * scale
        s = level.sigma * scale
        la = lat[level.indices]
        before = (level.y - driver_hyp1_equivalent(problem, level, hyp1, x, order)) * scale
        ax.errorbar(la, r, yerr=s, fmt="o", markersize=3, color="C3", label="hypothesis 3")
        ax.plot(la, before, "s", markersize=3, color="C2", label="hypothesis 1")
        ax.axhline(0, color="0.5", linewidth=0.7)
        ax.set_ylabel(f"{level.name} (K/deg)")
        ax.legend(fontsize=8)
        rms, rms_h1 = float(np.sqrt(np.mean(r ** 2))), float(np.sqrt(np.mean(before ** 2)))
        worst = np.argsort(-np.abs(r / s))[:3]
        summary["iris"][level.name] = {"rms_K_per_deg": rms, "median_sigma": float(np.median(s)), "rms_hyp1": rms_h1,
                                       "worst": [(float(la[i]), float(r[i])) for i in worst]}
        lines.append(f"  {level.name}: RMS {rms:.3f} (hypothesis 1 {rms_h1:.3f}), median stated error {np.median(s):.3f}; "
                     "largest at " + ", ".join(f"{la[i]:+.0f} deg {r[i]:+.3f}" for i in worst))
    axes[-1].set_xlabel("planetocentric latitude (deg)")
    fig.tight_layout()
    fig.savefig(HERE / "iris_residuals.png", dpi=110)
    plt.close(fig)

    # 110 mbar: where the fit lands. The radio's value is 0 here, all three being delivered minus observed.
    innovation = blocks.y - blocks.h + blocks.H @ (fit.state - problem.prior.mean)
    products_final = fit.final.products
    delivered_final = at_level(problem, products_final, blocks, blocks.h)
    without = {t.name: solver.leave_out(problem.prior, blocks, innovation, t.name) for t in problem.radio.targets}
    summary["at_110_mbar"] = {
        "sigma_T": float(problem.radio.sigma(np.array([AT_110_Pa]))[0]), "hypothesis_3": delivered_final,
        "hypothesis_1": {t: delivered[t]["hypothesis 1"][TABLE_MBAR.index(110.0)] for t in RUNS},
        "without_own_radio": {k: at_level(problem, products_final, blocks, blocks.h + blocks.H @ (v["state"] - fit.state))[k] for k, v in without.items()},
        "projection": at_level(problem, projected.products, projected.blocks, projected.blocks.h)}
    at110 = summary["at_110_mbar"]
    ruled = summary["landing"] = landing(problem, products_final, blocks, fit.state, problem.prior)
    lines.append(f"at 110 mbar, delivered minus observed (K; the radio is 0, sigma_T {at110['sigma_T']:.2f}), hypothesis 1, projection of hypothesis 1, IRIS only, without the target's own radio, hypothesis 3, first fit:")
    lines += [f"  {t:7s} {at110['hypothesis_1'][t]:+6.2f} {at110['projection'][t]:+6.2f} {ruled['iris_only'][t]:+6.2f} {at110['without_own_radio'][t]:+6.2f} {at110['hypothesis_3'][t]:+6.2f} {FIRST_FIT[t][TABLE_MBAR.index(110.0)]:+6.2f}" for t in RUNS]
    segments, shears = path_segments(problem, fit, start, hyp1)
    summary["path_at_110_mbar"] = {"edges_planetocentric": PATH_EDGES, "segments_K": segments}
    lines.append("hypotheses along the path at 110 mbar, K per segment, planetocentric latitude from the anchor (the delivered change to compare at the end):")
    for name, parts in segments.items():
        lines.append(f"  {name}: " + "  ".join(f"{a:.1f} to {b:.1f} {v:+.2f}" for a, b, v in zip(PATH_EDGES[:-1], PATH_EDGES[1:], parts))
                     + f"; to 3.0 S {sum(parts[:4]):+.2f}, to 31.2 S {sum(parts):+.2f}")
    lines.append(f"  delivered, fit minus hypothesis 1: 3.0 S {delivered_final['3.0 S'] - at110['hypothesis_1']['3.0 S']:+.2f}, 31.2 S {delivered_final['31.2 S'] - at110['hypothesis_1']['31.2 S']:+.2f}")
    picks = (15.0, 12.5, 10.0, 7.5, 5.0, 2.5, 0.0, -2.5, -5.0, -7.5, -10.0, -15.0)
    lines.append("  shear at 110 mbar, m/s per unit ln p, planetocentric latitude: " + " ".join(f"{v:g}" for v in picks))
    lines += [f"    {name:13s}" + " ".join(f"{np.interp(v, shears['lat'], shears[name]):+6.1f}" for v in picks) for name in ("hypothesis 1", "projection", "fit")]
    summary["shear_at_110_mbar"] = {name: [float(np.interp(v, shears["lat"], shears[name])) for v in picks] for name in ("hypothesis 1", "projection", "fit")}
    lines.append(f"figures {HERE / 'shear.png'}, {HERE / 'profiles.png'}, {HERE / 'targets.png'}, {HERE / 'iris_residuals.png'}")
    record(2, "for the author: hypothesis 3 beside hypotheses 1 and 2 and the first fit", None, "\n".join(lines))

    # ---- 3. the consistency check ---------------------------------------------------------------
    check = solver.consistency(blocks, solution)
    summary["consistency"] = check
    lines = ["group, rows, weighted misfit, expected, ratio:"]
    lines += [f"  {g:14s} {check[g]['count']:4d} {check[g]['misfit']:9.1f} {check[g]['expected']:9.1f} {check[g]['ratio']:7.2f}" for g in groups]
    lines.append(f"  all (innovations against S): {check['all']['misfit']:.1f} of {check['all']['expected']:.0f}, ratio {check['all']['ratio']:.2f}")
    summary["basis_limit"] = {}
    lines.append("hypothesis 1 projected on the basis, the prior not involved: weighted misfit over rows, the basis's own limit:")
    for group in groups:
        g = projected.blocks.groups == group
        r = (projected.blocks.y - projected.blocks.h)[g]
        summary["basis_limit"][group] = float(r @ np.linalg.solve(projected.blocks.R[np.ix_(g, g)], r) / g.sum())
        lines.append(f"  {group:14s} {summary['basis_limit'][group]:6.3f}")
    summary["radio_modes"] = radio_modes(problem, blocks, solution)
    lines.append("the radio misfit in the eigenmodes of its R, largest first: modes, misfit, expected, ratio:")
    for target, rows in summary["radio_modes"].items():
        lines.append(f"  {target}: " + "; ".join(f"{a} to {b}: {m:.2f}, {e:.2f}, {m / e:.2f}" for a, b, m, e in rows))
    # Whether the term of section 2 item 3 is sized: an IRIS group above 1.5, or a radio target missed by more
    # than 2 stated sd at an IRIS level when withheld.
    above = {g: check[g]["ratio"] for g in groups if g.startswith("IRIS") and check[g]["ratio"] > 1.5}
    misses = {}
    at_iris = np.array([level.pressure_Pa for level in problem.gradients.levels])
    for target in problem.radio.targets:
        out = solver.leave_out(problem.prior, blocks, innovation, target.name)
        table = observed[target.name]
        p_rows = np.asarray(products_final[target.name].to_dataset(inherit=False)["pressure_Pa"].values)
        p_rows = p_rows[np.isfinite(interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], p_rows))]
        z = interp_ln_p(p_rows, out["error"] / np.sqrt(out["variance"]), at_iris)
        for p, v in zip(at_iris, z):
            if abs(v) > 2.0:
                misses[f"{target.name} at {p / 100:g} mbar"] = float(v)
    summary["term_trigger"] = {"iris_above_1.5": above, "radio_missed_over_2_sd": misses}
    lines.append(f"the systematic term is asked for by: IRIS above 1.5 {({k: round(v, 2) for k, v in above.items()}) or 'none'}; "
                 f"leave-one-out over 2 sd at an IRIS level {({k: round(v, 2) for k, v in misses.items()}) or 'none'}")
    summary["terms"] = []
    if above or misses:
        lines.append("the candidate term, correlated in latitude, sized per level by one solve; reported, not used until the author rules:")
        for length in TERM_CORRELATIONS:
            term = size_terms(problem, products_final, blocks, fit.state, length)
            summary["terms"].append(term)
            lines.append(f"  correlation {length:g} deg: size " + ", ".join(f"{k} {v:.3f} K/deg" for k, v in term["sizes_K_per_deg"].items()) + " (nan: no size reaches the expected value)")
            lines.append("    ratios " + ", ".join(f"{k} {v:.2f}" for k, v in term["applied"]["ratios"].items()))
            lines.append("    at 110 mbar, all data " + ", ".join(f"{k} {v:+.2f}" for k, v in term["applied"]["all_data"].items()) + "; IRIS only " + ", ".join(f"{k} {v:+.2f}" for k, v in term["applied"]["iris_only"].items()))
            for name, row in term["scan"].items():
                lines.append(f"    scan, {name}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["ratios"].items() if k.startswith("IRIS"))
                             + f"; at 110 mbar, 31.2 S {row['all_data']['31.2 S']:+.2f}, 3.0 S {row['all_data']['3.0 S']:+.2f}")
    record(3, "the consistency check, per group", None, "\n".join(lines))

    # ---- 4. leave-one-out -----------------------------------------------------------------------
    lines = ["each group withheld, predicted by the fit to the rest on the final linearization: weighted error over rows (1 if "
             "the stated errors and the prior are right); RMS prediction error; RMS stated (data and prediction):"]
    summary["leave_one_out"] = {}
    for group in groups:
        out = solver.leave_out(problem.prior, blocks, innovation, group)
        g = blocks.groups == group
        radio_group = group in {t.name for t in problem.radio.targets}
        if radio_group:
            scale = np.exp(blocks.y[g])
            unit = "K"
        else:
            level = next(v for v in problem.gradients.levels if v.name == group)
            scale = level.temperature_K * np.pi / 180.0
            unit = "K/deg"
        rms_e, rms_s = float(np.sqrt(np.mean((out["error"] * scale) ** 2))), float(np.sqrt(np.mean(out["variance"] * scale ** 2)))
        summary["leave_one_out"][group] = {"weighted_over_rows": out["weighted"] / out["count"], "rms_error": rms_e, "rms_stated": rms_s, "unit": unit}
        lines.append(f"  {group:14s} {out['weighted'] / out['count']:6.2f}   {rms_e:7.3f} {unit}   {rms_s:7.3f} {unit}")
        if radio_group:
            p_rows = np.asarray(fit.final.products[group].to_dataset(inherit=False)["pressure_Pa"].values)
            table = observed[group]
            keep = np.isfinite(interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], p_rows))
            e = interp_ln_p(p_rows[keep], out["error"] * scale, np.array(TABLE_MBAR) * 100.0)
            s = interp_ln_p(p_rows[keep], np.sqrt(out["variance"]) * scale, np.array(TABLE_MBAR) * 100.0)
            lines.append("      observed minus predicted at table levels: " + " ".join(f"{v:+.1f}({u:.1f})" for v, u in zip(e, s)))
            summary["leave_one_out"][group]["at_table"] = e.tolist()
    record(4, "leave-one-out", None, "\n".join(lines))

    # ---- 5. the solve ---------------------------------------------------------------------------
    lines = [f"{fit.passes} passes; the largest change in delivered T at each (K): "]
    lines += [f"  pass {i + 1}: " + ", ".join(f"{k} {v:.3f}" for k, v in c.items()) for i, c in enumerate(fit.changes_K)]
    identity = {}
    for target, runs in RUNS.items():
        residual = np.asarray(roots[runs[2]]["pressure_identity_residual"].values, dtype="float64")
        identity[target] = float(np.nanmax(np.abs(residual)))
        lines.append(f"  pressure identity, {target}: largest {identity[target]:.3e}")
    reproduced = {}
    for t, runs in RUNS.items():
        reproduced[t] = float(np.nanmax(np.abs(roots[runs[2]]["temperature_K"].values - fit.final.products[t].to_dataset(inherit=False)["temperature_K"].values)))
    lines.append("  registered runs against the fit's last runs, largest T difference: " + ", ".join(f"{k} {v:.2e} K" for k, v in reproduced.items()))
    from_mean = driver.iterate(problem, problem.prior.mean, problem.directory / "work_prior_mean", pool)
    shutil.rmtree(problem.directory / "work_prior_mean", ignore_errors=True)
    final_T, mean_T = driver._delivered(fit.final.products, problem.radio), driver._delivered(from_mean.final.products, problem.radio)
    starts = {k: float(np.nanmax(np.abs(final_T[k] - mean_T[k]))) for k in final_T}
    lines.append(f"  started from the prior mean: {from_mean.passes} passes, the changes " + "; ".join(", ".join(f"{v:.3f}" for v in c.values()) for c in from_mean.changes_K))
    lines.append("  largest difference in delivered T from the fit's: " + ", ".join(f"{k} {v:.3f} K" for k, v in starts.items()) + " (bound 0.1 K)")
    lines.append(f"  `casspian-wind-fit` itself against the same fit in this script, largest wind difference {command_difference:.2e} m/s")
    lines.append(f"  condition number of the radio R: {np.linalg.cond(fit.final.radio_block.R):.0f}")
    summary["solve"] = {"passes": fit.passes, "changes": fit.changes_K, "identity": identity, "reproduced": reproduced,
                        "from_prior_mean": {"passes": from_mean.passes, "changes": from_mean.changes_K, "largest_difference_K": starts},
                        "command_difference_ms": command_difference}
    converged = max(fit.changes_K[-1].values()) < problem.tolerance_K
    record(5, "the solve: the change at each relinearization, the pressure identity of each run, and the start from the prior mean within 0.1 K",
           converged and max(identity.values()) < 1e-2 and max(starts.values()) < 0.1 and command_difference == 0.0, "\n".join(lines))

    # ---- 6. the wind aloft -----------------------------------------------------------------------
    written = driver._read(problem.output, "wind")
    u, sig = written["u_total_ms"].values, written["u_total_uncertainty_ms"].values
    reference = float(written["reference_level_pressure_Pa"])
    ref_col = int(np.flatnonzero(problem.pressure == reference)[0])
    base_u = problem.wind["u_total_ms"].values
    in_range = (problem.pressure > SHEAR_RANGE_Pa[0]) & (problem.pressure < reference)
    change = np.abs(u - base_u) * in_range[None, :]
    i, j = np.unravel_index(np.argmax(change), change.shape)
    summary["wind"] = {"largest_change": float(change[i, j]), "at_planetocentric": float(lat[i]), "at_mbar": float(problem.pressure[j] / 100),
                       "sigma_there": float(sig[i, j]), "sigma_max": float(np.max(sig)), "by_level": {}}
    lines = [f"largest change from the base wind above the reference level {change[i, j]:.1f} m/s at {lat[i]:.1f} deg planetocentric, "
             f"{problem.pressure[j] / 100:.2f} mbar, against a posterior sd there of {sig[i, j]:.1f} m/s",
             "by level, mbar: largest change (m/s), at planetocentric, its posterior sd, change over sd; RMS change; median sd:"]
    for mbar in TABLE_MBAR:
        k = int(np.argmin(np.abs(np.log(problem.pressure / (mbar * 100.0)))))
        d = np.abs(u[:, k] - base_u[:, k])
        n = int(np.argmax(d))
        summary["wind"]["by_level"][mbar] = {"largest": float(d[n]), "at_planetocentric": float(lat[n]), "sigma_there": float(sig[n, k]),
                                             "rms": float(np.sqrt(np.mean(d ** 2))), "sigma_median": float(np.median(sig[:, k]))}
        lines.append(f"  {mbar:g}: {d[n]:.1f}, {lat[n]:.1f}, {sig[n, k]:.1f}, {d[n] / max(sig[n, k], 1e-9):.1f}; {np.sqrt(np.mean(d ** 2)):.1f}; {np.median(sig[:, k]):.1f}")
    ok = np.abs(u[[0, -1]]).max() == 0.0 and np.max(np.abs(u[:, ref_col] - problem.wind["u_reference_ms"].values)) == 0.0
    lines.append(f"poles zero and the reference level unchanged: {bool(ok)}")
    lines.append("the layer operator at hypothesis 1's wind and shear against the data, which hypothesis 1 solved at a point, per radian:")
    for level in problem.gradients.levels:
        r = level.y - driver_hyp1_equivalent(problem, level, hyp1, x, order)
        lines.append(f"  {level.name}: RMS {np.sqrt(np.mean(r ** 2)):.4f} against a median error {np.median(level.sigma):.4f}")
    record(6, "the wind aloft: the largest change from the base wind, and where, against its posterior sd; the wind file holds", bool(ok), "\n".join(lines))

    pool.shutdown()
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    failed = [r for r in results if r[2] is False]
    bounded = [r for r in results if r[2] is not None]
    print(f"{len(bounded) - len(failed)} of {len(bounded)} bounded checks pass; {len(results) - len(bounded)} reported without a bound")
    return 1 if failed else 0


def driver_hyp1_equivalent(problem, level, hyp1, x, order):
    """The layer operator at hypothesis 1's wind and shear, at `level`'s points."""
    nodes = layer_pressures(level.pressure_Pa)
    unit = lw.balance_weight(problem.wind, problem.composition, problem.gravity, problem.rotation, nodes, np.ones((nodes.size, problem.lat.size)))
    out = np.zeros(level.indices.size)
    for q, quadrature in enumerate(LAYER_WEIGHTS):
        term = wind_term(problem.wind, problem.gravity, problem.rotation, nodes[q])
        s_level = np.array([np.interp(np.log(nodes[q]), x[order], hyp1[order, j]) for j in range(problem.lat.size)])
        out += quadrature * (unit[q, level.indices] / level.temperature_K * s_level[level.indices] + term[level.indices])
    return out


if __name__ == "__main__":
    sys.exit(main())
