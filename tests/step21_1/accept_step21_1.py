"""SPEC_21 Step 1 acceptance: the general error budget, and hypothesis 4.

The synthetic tests first (1a resolution, 1b a balanced equatorial anomaly, 1c the linear statistics, 1d the
prior's draws); then `casspian-wind-fit`'s fit, `voyager_fit_2`, with its runs, its anchor matrix and its sections;
the checks; the solve; and the fit again with the Voyager 2 egress as the reference. Writes under
`reports/step21_1/`, `reports/figures/step21_1_*.png`, the wind file, the runs and the sections.

Figures in the standing style: absolute profiles first, both hydrogen adiabats on every temperature panel; sections
with the run's zonal wind shaded under temperature contours, no hatching.

SPEC_21 v0.4 and v0.5: the prior's latitude length floored at the deformation radius and no superadiabatic delivered
profile (REVIEW_21_step1_preliminary); 1(d) also reports the equatorial deformation radius used and where the
stability constraint was active.

REVIEW_21_step1: 1(c)'s bound is judged on the mean ratios over RATIO_DRAWS draws, the 20-draw values reported; the
fit's wind files reach a decade above 1 Pa (`driver.extend_top`); check 2 adds the equatorial wind aloft with its
posterior sd beside hypothesis 3's.

Run from the repository root; about six hours alone. `--resume` starts after the matrix runs from what the last
run left in `reports/step21_1/resume.pkl` (REVIEW_21_step1_preliminary, the departures).
"""

import json
import pickle
import shutil
import subprocess
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "step20_1"))
import build_matrix as bm  # noqa: E402
from casspian.forward import transfer as forward_transfer  # noqa: E402
from casspian.lib import heat_capacity as hc  # noqa: E402
from casspian.lib import io as cio  # noqa: E402
from casspian.tools.run import new_run, run_inputs  # noqa: E402
from casspian.tools.wind.estimation import driver, solver  # noqa: E402
from casspian.tools.wind.estimation.gradient import Systematic  # noqa: E402
from casspian.tools.wind.estimation.stability import GROUP as STABILITY  # noqa: E402
from casspian.tools.wind.estimation.prior import fine_grid  # noqa: E402
from casspian.tools.wind.estimation.radio import interp_ln_p  # noqa: E402

warnings.filterwarnings("ignore")
HERE = Path("reports/step21_1")
FIG = Path("reports/figures")
LINDAL = Path("occul_data/lindal")
CONTROL = Path("fits/voyager_fit_2/voyager_fit_2.toml")
CONTROL_31S = Path("fits/voyager_fit_2/voyager_fit_2_reference_31s.toml")
HYP3_WIND = Path("fits/voyager_fit_1/voyager_fit_1_wind.nc")
SEED = 21
DRAWS = 20
#: 1(c), reported beside the bound: the ratios alone over more draws, since a group with few effective degrees of
#: freedom has a 20-draw mean too wide for the bound (REVIEW_21_step1_preliminary, diagnosis 4.2).
RATIO_DRAWS = 500
TABLE_MBAR = (730.0, 290.0, 110.0, 30.0, 10.0, 1.0)
ANCHORS = ("lindal", "voyager2_egress", "voyager1_ingress", "voyager1_egress")
NAME = {"lindal": "36.3 N", "voyager2_egress": "31.2 S", "voyager1_ingress": "74.2 S", "voyager1_egress": "3.0 S"}
LABEL = {"lindal": "36.3 N (30.8 N), Lindal", "voyager2_egress": "31.2 S (26.1 S), V2 egress",
         "voyager1_ingress": "74.2 S (71.1 S), V1 ingress", "voyager1_egress": "3.0 S (2.4 S), V1 egress"}
COLOR = {"lindal": "C0", "voyager2_egress": "C3", "voyager1_ingress": "C2", "voyager1_egress": "C4"}
#: Hypothesis 4's runs from the Lindal anchor, made here from SPEC_19's base runs with only the wind changed.
H4_RUNS = {"voyager2_egress": ("lindal_voyager_fit_2", "lindal_iris_v_sl_ii_110"),
           "voyager1_ingress": ("lindal_voyager_fit_2_74s", "lindal_iris_v_sl_ii_110_74s"),
           "voyager1_egress": ("lindal_voyager_fit_2_03s", "lindal_iris_v_sl_ii_110_03s")}
H3_RUNS = bm.LINDAL_RUNS
SUFFIX_H4 = "_h4"
#: 1(b): a 20 K warm layer, half-width 10 deg, about one scale height thick near 15 mbar (SPEC_21 section 4).
ANOMALY_K, ANOMALY_HALF_WIDTH_DEG, ANOMALY_CENTER_PA, ANOMALY_HALF_THICKNESS = 20.0, 10.0, 1500.0, 0.5
RESUME = HERE / "resume.pkl"
results, summary = [], {}


def record(number, description, passed, detail):
    results.append(passed)
    print(f"[{'REPORTED' if passed is None else 'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")
    sys.stdout.flush()


def tree(path):
    return driver.tree(path)


def composition(slug):
    handle = cio.read(LINDAL / f"{slug}_composition.nc", "composition")
    try:
        return handle.load()
    finally:
        handle.close()


def latitude(slug):
    return float(tree(LINDAL / f"{slug}_refractivity.nc").to_dataset(inherit=False)["latitude_planetocentric_deg"])


def ratios(check):
    return ", ".join(f"{k} {v['ratio']:.2f}" for k, v in check.items())


# ---- the synthetic fits -----------------------------------------------------------------------------------

def synthetic_fit(problem, truth, rng, pool, label):
    """Data from `truth` at the real observation points with noise from `R` there, fitted as the tool fits: from
    the prior mean, IRIS's systematic parts sized, and again once. Returns the fit and the truth's linearization."""
    work = problem.directory / f"work_{label}"
    at_truth = driver.linearize(problem, truth, work, pool)
    b = observations(at_truth.blocks)
    data = b.h + np.linalg.cholesky(b.R) @ rng.standard_normal(b.h.size)
    problem.gradients.systematic = {}
    try:
        first = driver.iterate(problem, problem.prior.mean, work, pool, data)
        sizes = driver.size_systematic(problem, first)
        problem.gradients.systematic = {k: Systematic(v, problem.systematic_length_deg, problem.systematic_form)
                                        for k, v in sizes.items() if np.isfinite(v) and v > 0}
        fit = driver.iterate(problem, first.state, work, pool, data) if problem.gradients.systematic else first
    except (RuntimeError, ValueError) as error:
        # Not converged, or a march refused: recorded as the result, not a crash of the whole acceptance.
        print(f"{label}: the synthetic fit did not finish: {str(error)[:300]}", flush=True)
        fit, first, sizes = None, None, {}
    finally:
        problem.gradients.systematic = {}
        shutil.rmtree(work, ignore_errors=True)
    return fit, first, sizes


def observations(blocks):
    """The observation rows of a linearization, the stability constraint's left out."""
    keep = blocks.groups != STABILITY
    return solver.Blocks(blocks.y[keep], blocks.h[keep], blocks.H[keep], blocks.R[np.ix_(keep, keep)], blocks.groups[keep])


def shear_sd(problem, solution, lat, pressure):
    a, b = problem.basis.lat_matrix(lat), problem.basis.x_matrix(pressure)
    out = np.empty((pressure.size, lat.size))
    for k in range(pressure.size):
        m = problem.basis.row(a, np.repeat(b[k:k + 1], lat.size, axis=0))
        out[k] = np.sqrt(np.maximum(solver.posterior_variance(problem.prior, solution, m), 0.0))
    return out


def anomaly_state(problem):
    """1(b)'s truth: the shear of a balanced warm layer by the prior's own construction, `s = (1/T_bar) dT'/dphi / w`,
    with `w` and `T_bar` as the prior holds them, projected on the basis."""
    lat_f, x_f = fine_grid(problem.basis)
    p_f = np.exp(x_f)
    with xr.open_dataset("occul_data/lindal/lindal_thermo.nc") as th:
        T_bar = interp_ln_p(th["pressure_Pa"].values, th["temperature_K"].values, p_f)
    T_bar = np.interp(x_f, x_f[np.isfinite(T_bar)], T_bar[np.isfinite(T_bar)])
    width = ANOMALY_HALF_WIDTH_DEG / np.sqrt(2.0 * np.log(2.0))
    thick = ANOMALY_HALF_THICKNESS / np.sqrt(2.0 * np.log(2.0))
    phi = np.radians(lat_f)
    vertical = np.exp(-0.5 * ((x_f - np.log(ANOMALY_CENTER_PA)) / thick) ** 2)
    dT_dphi = -ANOMALY_K * (lat_f / width ** 2) * (180.0 / np.pi) * np.exp(-0.5 * (lat_f / width) ** 2)
    g = (vertical / T_bar)[:, None] * dT_dphi[None, :]
    lat_w = problem.lat
    with ProcessPoolExecutor(max_workers=driver.WORKERS) as pool:
        chunks = np.array_split(np.arange(p_f.size), driver.WORKERS)
        unit = np.vstack(list(pool.map(driver._weight_job, [(problem.wind, problem.composition, problem.gravity,
                                                             problem.rotation, p_f[c], lat_w.size) for c in chunks if c.size])))
    w = np.stack([np.interp(lat_f, lat_w, row) for row in unit / T_bar[:, None]])
    shear = g / w
    return problem.basis.project(shear, lat_f, p_f)


def report_resolution(problem, label, truth, fit, sizes, lat_s, p_s, number, description):
    if fit is None:
        summary[label] = {"finished": False}
        record(number, description, None, "the synthetic fit did not finish; see the log above")
        tru = problem.basis.shear(truth, lat_s, p_s)
        return np.full_like(tru, np.nan), tru, np.full_like(tru, np.nan)
    est = problem.basis.shear(fit.state, lat_s, p_s)
    tru = problem.basis.shear(truth, lat_s, p_s)
    sd = shear_sd(problem, fit.solution, lat_s, p_s)
    z = (est - tru) / sd
    check = solver.consistency(fit.weighted, fit.solution, exclude=(STABILITY,))
    j = int(np.argmin(np.abs(p_s - 1.1e4)))
    k = int(np.argmin(np.abs(p_s - ANOMALY_CENTER_PA)))
    picks = np.array([15.0, 10.0, 5.0, 2.5, 0.0, -2.5, -5.0, -10.0, -15.0])
    idx = [int(np.argmin(np.abs(lat_s - v))) for v in picks]
    rows = {"truth 110": tru[j, idx], "fit 110": est[j, idx], "sd 110": sd[j, idx],
            "truth 15": tru[k, idx], "fit 15": est[k, idx], "sd 15": sd[k, idx]}
    summary[label] = {"passes": fit.passes, "rms_error": float(np.sqrt(np.mean((est - tru) ** 2))), "rms_true": float(np.sqrt(np.mean(tru ** 2))),
                      "within_1": float(np.mean(np.abs(z) <= 1)), "within_2": float(np.mean(np.abs(z) <= 2)),
                      "ratios": {g: v["ratio"] for g, v in check.items()}, "sizes": sizes,
                      "rows": {k2: v.tolist() for k2, v in rows.items()}, "stability": fit.final.active}
    lines = [f"{fit.passes} passes in its last fit; IRIS systematic parts sized at " + ", ".join(f"{k2} {v:.4g}" for k2, v in sizes.items()),
             f"RMS shear error {summary[label]['rms_error']:.2f} against a true RMS {summary[label]['rms_true']:.2f} m/s per unit ln p; "
             f"within 1 posterior sd {100 * summary[label]['within_1']:.1f}%, within 2 {100 * summary[label]['within_2']:.1f}%",
             "misfit over expected: " + ratios(check),
             f"stability constraint active at the fit: {fit.final.active}",
             "shear, m/s per unit ln p, planetocentric " + " ".join(f"{v:g}" for v in picks) + ":"]
    lines += [f"  {name:9s}" + "".join(f"{v:+8.1f}" for v in vals) for name, vals in rows.items()]
    record(number, description, None, "\n".join(lines))
    return est, tru, sd


# ---- the figures ------------------------------------------------------------------------------------------

def draw_adiabats(ax, slug, p_mbar, T, latitude_deg):
    deepest = int(np.argmax(p_mbar))
    for mode, style in (("normal", "--"), ("equilibrium", ":")):
        p, t = hc.adiabat(float(T[deepest]), float(p_mbar[deepest]) * 100.0, 1.0e4, composition(slug), latitude_deg, mode)
        ax.plot(t, p / 100.0, color="0.35", linestyle=style, linewidth=1.1, label=f"dry adiabat, {mode} H2")


def section_figure(sections, lats, path, title):
    grid = np.geomspace(1.3e5, 20.0, 160)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
    levels_wind = np.concatenate([np.arange(-25, 0, 5), np.arange(0, 551, 25)])
    norm = matplotlib.colors.TwoSlopeNorm(vmin=-25, vcenter=0, vmax=550)
    for ax, slug in zip(axes.ravel(), ANCHORS):
        ds = sections[slug]
        lat = ds["latitude_planetocentric_deg"].values
        T = np.array([interp_ln_p(ds["pressure_Pa"].values[i], ds["temperature_K"].values[i], grid) for i in range(lat.size)])
        u = np.array([interp_ln_p(ds["pressure_Pa"].values[i], ds["u_ms"].values[i], grid) for i in range(lat.size)])
        cs = ax.contourf(lat, grid / 100, u.T, levels=levels_wind, cmap="RdBu_r", norm=norm, extend="both")
        ct = ax.contour(lat, grid / 100, T.T, levels=np.arange(60, 181, 5), colors="k", linewidths=0.6)
        ax.clabel(ct, levels=np.arange(60, 181, 10), fmt="%d K", fontsize=7, inline_spacing=2)
        for s in ANCHORS:
            ax.axvline(lats[s], color="k", linewidth=0.7, linestyle="--")
        ax.axvline(-9.5, color="k", linewidth=0.7, linestyle=":")
        ax.set_yscale("log")
        ax.set_ylim(1300, 0.2)
        ax.set_title(f"{LABEL[slug]} anchor", fontsize=10)
    for ax in axes[1]:
        ax.set_xlabel("planetocentric latitude (deg)")
    for ax in axes[:, 0]:
        ax.set_ylabel("pressure (mbar)")
    fig.colorbar(cs, ax=axes, label="zonal wind of the run (m/s); blue and red on separate scales", shrink=0.8,
                 ticks=[-25, -20, -15, -10, -5, 0, 100, 200, 300, 400, 500])
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=95)
    plt.close(fig)


def make_run(name, base, wind_path):
    """A run from `base` with only the wind changed, the fit's wind by the identity shear (as SPEC_19)."""
    directory = Path("forward") / name
    if directory.exists():
        shutil.rmtree(directory)
    new_run.new_run(name, Path("forward") / base)
    build = directory / f"{name}_build.toml"
    text = build.read_bytes().decode("utf-8")
    start, stop = text.index("\n[shear]") + 1, text.index("\n[composition]") + 1
    text = (text[:start] + "[shear]\n# SPEC_21 Step 1: hypothesis 4, the fit's wind; the only change from " + base + ".\n"
            f'source = "../../{wind_path.as_posix()}"\noutput = "inputs/{name}_wind.nc"\ncase   = "identity"\n\n' + text[stop:])
    build.write_bytes(text.encode("utf-8"))
    run_inputs.build(build)
    return forward_transfer.run(directory / f"{name}.toml").product


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    lats = {s: latitude(s) for s in ANCHORS}
    rng = np.random.default_rng(SEED)
    problem = driver.read_control(CONTROL)
    lat_s = np.arange(-87.5, 87.6, 2.5)
    p_s = np.geomspace(1.3e5, 20.0, 60)
    pool = ProcessPoolExecutor(max_workers=driver.WORKERS)

    # ---- 1(a), 1(b) ----------------------------------------------------------------------------------------
    truth_a = problem.basis.project(driver.base_shear(problem), problem.lat, problem.pressure)
    fit_a, _, sizes_a = synthetic_fit(problem, truth_a, rng, pool, "synthetic_a")
    est_a, tru_a, sd_a = report_resolution(problem, "1a", truth_a, fit_a, sizes_a, lat_s, p_s, "1a",
                                           "resolution, hypothesis 1 as the truth; no bound")
    truth_b = truth_a + anomaly_state(problem)
    fit_b, _, sizes_b = synthetic_fit(problem, truth_b, rng, pool, "synthetic_b")
    est_b, tru_b, sd_b = report_resolution(problem, "1b", truth_b, fit_b, sizes_b, lat_s, p_s, "1b",
                                           "a balanced 20 K warm layer at the equator near 15 mbar, on hypothesis 1; no bound")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
    for row, (label, tru, est, sd) in enumerate((("1(a)", tru_a, est_a, sd_a), ("1(b)", tru_b, est_b, sd_b))):
        for ax, field, title, cmap, levels in ((axes[row, 0], tru, "truth", "RdBu_r", np.arange(-150, 151, 10)),
                                               (axes[row, 1], est, "recovered", "RdBu_r", np.arange(-150, 151, 10)),
                                               (axes[row, 2], (est - tru) / sd, "(recovered - truth) / posterior sd", "RdBu_r", np.arange(-3, 3.1, 0.5))):
            cs = ax.contourf(lat_s, p_s / 100, field, levels=levels, cmap=cmap, extend="both")
            fig.colorbar(cs, ax=ax)
            ax.set_yscale("log")
            ax.set_ylim(1300, 0.2)
            ax.set_title(f"{label}: shear {title}", fontsize=10)
    for ax in axes[1]:
        ax.set_xlabel("planetocentric latitude (deg)")
    fig.savefig(FIG / "step21_1_synthetic.png", dpi=90)
    plt.close(fig)

    # ---- hypothesis 4 ------------------------------------------------------------------------------------------
    problem, fit, rec = driver.run(CONTROL, pool)
    first = rec["first"]
    wind_path = problem.output.relative_to(Path.cwd())
    products = {}
    for target, (name, base) in H4_RUNS.items():
        products[target] = make_run(name, base, wind_path)
    final = fit.final.blocks
    solution, weighted, weights = solver.robust(problem.prior, fit.state, final, problem.huber_threshold,
                                                exempt=(STABILITY,))

    # ---- 1(c), 1(d): the prior's draws -------------------------------------------------------------------------
    draws = problem.prior.draw(np.random.default_rng([SEED, 3]), DRAWS)
    linear = observations(final)
    noise_L = np.linalg.cholesky(linear.R)
    pooled, count, group_ratios = 0, 0, {}
    for i, truth in enumerate(draws):
        r = np.random.default_rng([SEED, 4, i])
        y = linear.h + linear.H @ (truth - fit.state) + noise_L @ r.standard_normal(linear.h.size)
        blocks = solver.Blocks(y, linear.h, linear.H, linear.R, linear.groups)
        sol = solver.linearized(problem.prior, fit.state, blocks)
        est, tru = problem.basis.shear(sol.state, lat_s, p_s), problem.basis.shear(truth, lat_s, p_s)
        sd = shear_sd(problem, sol, lat_s, p_s)
        pooled += int(np.sum(np.abs(est - tru) <= 2 * sd))
        count += est.size
        for g, v in solver.consistency(blocks, sol).items():
            group_ratios.setdefault(g, []).append(v["ratio"])
    mean_ratio = {g: float(np.mean(v)) for g, v in group_ratios.items()}
    spread = {g: float(np.std(v, ddof=1)) for g, v in group_ratios.items()}
    # The ratios' bound is judged on the mean over RATIO_DRAWS draws, ratios only (REVIEW_21_step1); the 20-draw
    # means are reported with each group's effective degrees of freedom and its 20-draw mean's sd under the model.
    from scipy.linalg import cho_factor, cho_solve, eigh
    S = linear.H @ problem.prior.apply(linear.H).T + linear.R
    RSR = linear.R @ cho_solve(cho_factor(S, lower=True), linear.R)
    dof = {}
    for g in dict.fromkeys(linear.groups.tolist()):
        m = linear.groups == g
        lam = eigh(RSR[np.ix_(m, m)], linear.R[np.ix_(m, m)], eigvals_only=True)
        lam = lam[lam > 0]
        dof[g] = (float(lam.sum() ** 2 / np.sum(lam ** 2)), float(np.sqrt(2 * np.sum(lam ** 2)) / lam.sum() / np.sqrt(DRAWS)))
    many = problem.prior.draw(np.random.default_rng([SEED, 5]), RATIO_DRAWS)
    many_ratios = {}
    for i, truth in enumerate(many):
        r = np.random.default_rng([SEED, 6, i])
        y = linear.h + linear.H @ (truth - fit.state) + noise_L @ r.standard_normal(linear.h.size)
        blocks = solver.Blocks(y, linear.h, linear.H, linear.R, linear.groups)
        for g, v in solver.consistency(blocks, solver.linearized(problem.prior, fit.state, blocks)).items():
            many_ratios.setdefault(g, []).append(v["ratio"])
    many_mean = {g: (float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(len(v)))) for g, v in many_ratios.items()}
    passed = pooled / count >= 0.9 and all(0.8 <= v[0] <= 1.2 for v in many_mean.values())
    summary["1c"] = {"within_2": pooled / count, "mean_ratio": mean_ratio, "spread": spread, "dof": dof,
                     f"mean_ratio_{RATIO_DRAWS}": many_mean, "condition_S": float(np.linalg.cond(S))}
    record("1c", f"statistics, linear: {DRAWS} truths from the prior, data by the operator at the final linearization; at least 90% "
           f"of the shear within 2 posterior sd; each group's mean ratio over {RATIO_DRAWS} draws 0.8 to 1.2", passed,
           f"pooled within 2 sd {100 * pooled / count:.1f}%\neach group's {DRAWS}-draw mean ratio (spread), reported: "
           + ", ".join(f"{g} {mean_ratio[g]:.2f} ({spread[g]:.2f})" for g in mean_ratio)
           + "\nreported beside the bound: each group's effective degrees of freedom, and its 20-draw mean's sd under the model: "
           + ", ".join(f"{g} {v[0]:.1f} ({v[1]:.2f})" for g, v in dof.items())
           + f"\nthe mean ratio over {RATIO_DRAWS} draws, ratios only (its standard error), judged: "
           + ", ".join(f"{g} {v[0]:.2f} ({v[1]:.2f})" for g, v in many_mean.items())
           + f"\nthe condition number of S {summary['1c']['condition_S']:.3g}")

    # 1(d): the draws' temperature, recovered from their shear by the balance, and their wind at the equator.
    lat_f, x_f = fine_grid(problem.basis)
    p_d = np.array([1.0e5, 1.1e4, 1.0e3, 1.0e2])
    with xr.open_dataset("occul_data/lindal/lindal_thermo.nc") as th:
        T_bar = interp_ln_p(th["pressure_Pa"].values, th["temperature_K"].values, p_d)
    chunks = np.array_split(np.arange(p_d.size), min(driver.WORKERS, p_d.size))
    unit = np.vstack(list(pool.map(driver._weight_job, [(problem.wind, problem.composition, problem.gravity, problem.rotation,
                                                         p_d[c], problem.lat.size) for c in chunks])))
    w = unit / T_bar[:, None]
    lines = ["the largest |T'| of the draws (K), on an isobar, the field's mean along it removed, by the balance "
             "T' = T_bar int w s dphi; sigma_T last:"]
    at = {**{NAME[s]: lats[s] for s in ANCHORS}, "87.5 N": 87.5, "87.5 S": -87.5}
    sigma_T = driver.amplitude(problem_control(CONTROL)["prior"], "sigma_T")(p_d)
    phi = np.radians(problem.lat)
    summary["1d"] = {"T": {}, "equator": {}}
    for k, p in enumerate(p_d):
        Tp = []
        for truth in draws:
            s = problem.basis.shear(truth, problem.lat, np.array([p]))[0]
            g = w[k] * s
            integral = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(phi))]) * T_bar[k]
            Tp.append(integral - np.mean(integral))
        Tp = np.abs(np.array(Tp))
        row = {name: float(np.max(Tp[:, int(np.argmin(np.abs(problem.lat - v)))])) for name, v in at.items()}
        summary["1d"]["T"][f"{p / 100:g}"] = row
        lines.append(f"  {p / 100:7g} mbar: " + ", ".join(f"{n} {v:.1f}" for n, v in row.items()) + f"; sigma_T {sigma_T[k]:.1f}")
    s_eq = np.array([problem.basis.shear(truth, np.array([0.0]), p_s)[:, 0] for truth in draws])
    du = np.array([problem.basis.wind_change(truth, np.array([0.0]), np.array([1.0e3, 1.0e5]))[0] for truth in draws])
    change = du[:, 0] - du[:, 1]
    summary["1d"]["equator"] = {"largest_shear": float(np.max(np.abs(s_eq))), "wind_change_10_to_1000": change.tolist()}
    lines.append(f"at the equator: the largest |shear| of the draws {np.max(np.abs(s_eq)):.0f} m/s per unit ln p; the wind change "
                 f"from 1000 to 10 mbar, by draw (m/s): " + " ".join(f"{v:+.0f}" for v in change))
    # The equatorial deformation radius the prior used, and where the stability constraint was active.
    lengths = problem.prior.lengths_deg
    if lengths is not None:
        eq = int(np.argmin(np.abs(lat_f)))
        used = lengths[:, eq]
        summary["1d"]["equator_length_deg"] = dict(zip([f"{np.exp(v) / 100:.4g}" for v in x_f], used.tolist()))
        floor = lengths > float(problem_control(CONTROL)["prior"]["latitude_correlation_deg"]) + 1e-9
        reach = [float(np.max(np.abs(lat_f[row]))) if row.any() else 0.0 for row in floor]
        lines.append("the latitude length used at the equator (deg), by level (mbar): "
                     + ", ".join(f"{np.exp(v) / 100:.3g} {u:.2f}" for v, u in zip(x_f[::3], used[::3])))
        lines.append(f"  largest {used.max():.2f} deg at {np.exp(x_f[int(np.argmax(used))]) / 100:.3g} mbar, smallest "
                     f"{used.min():.2f} deg at {np.exp(x_f[int(np.argmin(used))]) / 100:.3g} mbar; the floor binds within "
                     f"{min(reach):.1f} to {max(reach):.1f} deg of the equator")
    lines.append(f"the stability constraint active at the hypothesis 4 fit: {fit.final.active}; at its first fit: {first.final.active}")
    summary["1d"]["stability"] = {"fit": fit.final.active, "first": first.final.active}
    record("1d", "the prior's draws are physical: temperature at the Voyager latitudes and 87.5 N and S, and the wind at the "
           "equator; the equatorial deformation radius used and where the stability constraint was active; no bound",
           None, "\n".join(lines))

    # ---- 2. hypothesis 4, absolute first -----------------------------------------------------------------------
    for (anchor, target), run in bm.runs({t: v[0] for t, v in H4_RUNS.items()}, SUFFIX_H4).items():
        if anchor != "lindal":
            bm.make((anchor, target), Path("forward/lindal_voyager_fit_2"), SUFFIX_H4, "hypothesis 4")
    matrix_h4 = bm.runs({t: v[0] for t, v in H4_RUNS.items()}, SUFFIX_H4)
    matrix_h3 = bm.runs()
    observed = problem.radio.observed
    reference = problem.radio.reference.name

    def obs(slug):
        o = observed[NAME[slug]]
        return o["pressure_mbar"], o["temperature_K"], o["temperature_uncertainty_K"]

    fig, axes = plt.subplots(1, 4, figsize=(22, 8.5), sharey=True)
    lines = [f"delivered minus observed (K) at {', '.join(f'{v:g}' for v in TABLE_MBAR)} mbar, hypothesis 4 then 3; the "
             "profile's own stated sigma last:"]
    summary["matrix"] = {}
    at_p = np.array(TABLE_MBAR) * 100.0
    for ax, target in zip(axes, ANCHORS):
        p_o, T_o, s_o = obs(target)
        ax.fill_betweenx(p_o, T_o - s_o, T_o + s_o, color="0.85", label="observed, its stated sigma")
        ax.plot(T_o, p_o, color="k", linewidth=1.6, label=f"observed, {LABEL[target]}")
        draw_adiabats(ax, target, p_o, T_o, lats[target])
        lines.append(f"  at {LABEL[target]}:")
        for anchor in ANCHORS:
            if anchor == target:
                continue
            for hyp, matrix, style in (("4", matrix_h4, "-"), ("3", matrix_h3, "--")):
                run = matrix[(anchor, target)]
                r = tree(f"forward/{run}/output/{run}_profile.nc").to_dataset(inherit=False)
                ax.plot(r["temperature_K"].values, r["pressure_Pa"].values / 100.0, color=COLOR[anchor], linestyle=style,
                        linewidth=1.2 if hyp == "4" else 0.8, label=f"from {LABEL[anchor]}, hypothesis {hyp}")
                d = interp_ln_p(r["pressure_Pa"].values, r["temperature_K"].values, at_p) - np.interp(np.log(at_p), np.log(p_o * 100), T_o)
                summary["matrix"][f"{anchor} to {target}, hypothesis {hyp}"] = d.tolist()
                lines.append(f"    from {anchor:17s} h{hyp}" + "".join(f"{v:+8.2f}" for v in d))
        lines.append(f"    {'sigma':24s}" + "".join(f"{v:8.2f}" for v in np.interp(np.log(at_p), np.log(p_o * 100), s_o)))
        ax.set_yscale("log")
        ax.set_ylim(1300, 0.2)
        ax.set_xlabel("temperature (K)")
        ax.set_title(f"at {LABEL[target]}", fontsize=10)
        ax.grid(True, color="0.92", linewidth=0.6)
        ax.legend(fontsize=6, loc="lower left")
    axes[0].set_ylabel("pressure (mbar)")
    fig.suptitle("The anchor matrix: hypothesis 4 solid, hypothesis 3 dashed; absolute profiles", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "step21_1_matrix_absolute.png", dpi=95)
    plt.close(fig)
    lines.append("  Voyager 1 egress to Voyager 1 ingress, hypothesis 4: "
                 + " ".join(f"{v:+.2f}" for v in summary["matrix"]["voyager1_egress to voyager1_ingress, hypothesis 4"]))
    context = {"problem": problem, "fit": fit, "rec": rec, "first": first, "products": products, "final": final,
               "solution": solution, "weighted": weighted, "weights": weights, "lines": lines, "lats": lats,
               "summary": summary, "results": results}
    try:
        with open(RESUME, "wb") as handle:
            pickle.dump(context, handle)
    except Exception as error:   # the resume point is a convenience; its absence only means a full rerun
        print(f"the resume point was not written: {error}", flush=True)
    try:
        return rest(context, pool)
    finally:
        pool.shutdown()


def rest(context, pool):
    """From the sections to the end: what follows the resume point."""
    problem, fit, rec, first = context["problem"], context["fit"], context["rec"], context["first"]
    products, final, solution = context["products"], context["final"], context["solution"]
    weighted, weights, lines, lats = context["weighted"], context["weights"], context["lines"], context["lats"]
    summary.update(context["summary"])
    results[:] = context["results"]
    lat_s = np.arange(-87.5, 87.6, 2.5)
    p_s = np.geomspace(1.3e5, 20.0, 60)
    at_p = np.array(TABLE_MBAR) * 100.0

    # The sections under hypothesis 4.
    sections = {}
    for slug in ANCHORS:
        name = f"{slug}_voyager_fit_2"
        directory = Path("sections") / name
        directory.mkdir(parents=True, exist_ok=True)
        text = Path(f"sections/{slug}_voyager_fit_1/{slug}_voyager_fit_1.toml").read_text(encoding="utf-8")
        text = text.replace("voyager_fit_1", "voyager_fit_2").replace("hypothesis 3", "hypothesis 4")
        (directory / f"{name}.toml").write_text(text, encoding="utf-8", newline="\n")
        subprocess.run(["casspian-section", str(directory / f"{name}.toml")], check=True, capture_output=True)
        sections[slug] = xr.open_dataset(directory / f"{name}_section.nc")
    section_figure(sections, lats, FIG / "step21_1_sections_temperature.png",
                   "Hypothesis 4 sections: temperature in black contours every 5 K, the run's zonal wind shaded. "
                   "Vertical lines: the four Voyager anchors (dashed), Pioneer 11 (dotted).")

    # The shear and the wind with their posterior sd, beside hypotheses 1 and 3.
    p_w = problem.pressure
    keep = (p_w <= 1.3e5) & (p_w >= 20.0)
    def on_grid(wind):
        """A wind on the base grid's levels: the fit's files reach above it (`driver.extend_top`)."""
        return wind.isel(pressure=np.searchsorted(np.asarray(wind["pressure_Pa"].values), p_w))

    files = {"hypothesis 1": driver._read(Path("forward/lindal_iris_v_sl_ii_110/inputs/lindal_iris_v_sl_ii_110_wind.nc"), "wind"),
             "hypothesis 3": driver._read(HYP3_WIND, "wind"), "hypothesis 4": driver._read(problem.output, "wind")}
    winds = {name: on_grid(wind) for name, wind in files.items()}

    def shear_of(wind):
        u = np.asarray(wind["u_total_ms"].values, dtype="float64")
        x = np.log(p_w)
        o = np.argsort(x)
        s = np.empty_like(u)
        s[:, o] = np.gradient(u[:, o], x[o], axis=1)
        return s

    sd4 = shear_sd(problem, solution, lat_s, p_s)
    fig, axes = plt.subplots(2, 4, figsize=(24, 11), sharey=True)
    for col, (name, wind) in enumerate(winds.items()):
        cs = axes[0, col].contourf(problem.lat, p_w[keep] / 100, shear_of(wind)[:, keep].T, levels=np.arange(-150, 151, 10), cmap="RdBu_r", extend="both")
        axes[0, col].set_title(f"{name}: shear (m/s per unit ln p)", fontsize=10)
        cw = axes[1, col].contourf(problem.lat, p_w[keep] / 100, wind["u_total_ms"].values[:, keep].T, levels=np.arange(-200, 601, 25), cmap="RdBu_r",
                                   norm=matplotlib.colors.TwoSlopeNorm(vmin=-200, vcenter=0, vmax=600))
        axes[1, col].set_title(f"{name}: wind (m/s)", fontsize=10)
    fig.colorbar(cs, ax=axes[0, :3])
    fig.colorbar(cw, ax=axes[1, :3])
    c1 = axes[0, 3].contourf(lat_s, p_s / 100, sd4, levels=np.geomspace(1, 1500, 20), cmap="magma_r", norm=matplotlib.colors.LogNorm())
    fig.colorbar(c1, ax=axes[0, 3])
    axes[0, 3].set_title("hypothesis 4: shear's posterior sd", fontsize=10)
    c2 = axes[1, 3].contourf(problem.lat, p_w[keep] / 100, winds["hypothesis 4"]["u_total_uncertainty_ms"].values[:, keep].T,
                             levels=np.geomspace(1, 3000, 20), cmap="magma_r", norm=matplotlib.colors.LogNorm())
    fig.colorbar(c2, ax=axes[1, 3])
    axes[1, 3].set_title("hypothesis 4: wind's posterior sd (m/s)", fontsize=10)
    for ax in axes.ravel():
        ax.set_yscale("log")
        ax.set_ylim(1300, 0.2)
    fig.savefig(FIG / "step21_1_shear_wind.png", dpi=85)
    plt.close(fig)
    j = int(np.argmin(np.abs(p_s - 1.1e4)))
    picks = np.array([15.0, 10.0, 5.0, 2.5, 0.0, -2.5, -5.0, -10.0, -15.0])
    idx = [int(np.argmin(np.abs(lat_s - v))) for v in picks]
    s4 = problem.basis.shear(fit.state, lat_s, np.array([1.1e4]))[0][idx]
    lines.append("shear at 110 mbar, planetocentric " + " ".join(f"{v:g}" for v in picks) + ":")
    for name, wind in winds.items():
        s = shear_of(wind)
        k = int(np.argmin(np.abs(p_w - 1.1e4)))
        lines.append(f"  {name:13s}" + "".join(f"{np.interp(v, problem.lat, s[:, k]):+8.1f}" for v in picks))
    lines.append(f"  {'h4 sd':13s}" + "".join(f"{v:8.1f}" for v in sd4[j, idx]))
    # The equatorial wind aloft with its posterior sd, hypothesis 4 beside 3: it carries the 3.0 S structure
    # (REVIEW_21_step1, check 5 item 3).
    aloft = np.array([30.0, 10.0, 3.0, 1.0])
    lines.append("the equatorial wind (m/s, its posterior sd in parentheses) at " + ", ".join(f"{v:g}" for v in aloft) + " mbar:")
    for name in ("hypothesis 3", "hypothesis 4"):
        wind = files[name]
        x = np.log(np.asarray(wind["pressure_Pa"].values))
        u0 = np.array([np.interp(0.0, problem.lat, col) for col in wind["u_total_ms"].values.T])
        sd0 = np.array([np.interp(0.0, problem.lat, col) for col in wind["u_total_uncertainty_ms"].values.T])
        cells = [f"{np.interp(np.log(v * 100.0), x, u0):+7.0f} ({np.interp(np.log(v * 100.0), x, sd0):.0f})" for v in aloft]
        lines.append(f"  {name:13s}" + "  ".join(cells))
        summary.setdefault("equatorial_wind", {})[name] = cells

    # IRIS residuals.
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for ax, level in zip(axes, problem.gradients.levels):
        g = final.groups == level.name
        scale = level.temperature_K * np.pi / 180.0
        r = (final.y - final.h)[g] * scale
        ax.errorbar(problem.lat[level.indices], r, yerr=level.sigma * scale, fmt="o", markersize=3, color="C3")
        ax.axhline(0, color="0.5", linewidth=0.7)
        ax.set_ylabel(f"{level.name} (K/deg)")
        lines.append(f"IRIS residuals at the fit, {level.name}: RMS {np.sqrt(np.mean(r ** 2)):.3f} K/deg, median stated "
                     f"{np.median(level.sigma * scale):.3f}")
    axes[-1].set_xlabel("planetocentric latitude (deg)")
    fig.tight_layout()
    fig.savefig(FIG / "step21_1_iris_residuals.png", dpi=100)
    plt.close(fig)
    record(2, "hypothesis 4, absolute first: the matrix beside hypothesis 3, the sections, the shear and wind, IRIS", None, "\n".join(lines))

    # ---- 3. the checks -----------------------------------------------------------------------------------------
    lines = ["misfit over expected, before IRIS's systematic parts: " + ratios(rec["before"]),
             "and after: " + ratios(rec["after"]),
             "the parts, per radian, latitude correlation "
             f"{problem.systematic_length_deg:g} deg {problem.systematic_form}: " + ", ".join(f"{k} {v:.4g}" for k, v in rec["sizes"].items())]
    innovation = solution.innovation
    summary["loo"] = {}
    lines.append("leave-one-out, weighted error over rows (1 if the stated errors and prior are right):")
    for group in dict.fromkeys(final.groups.tolist()):
        if group == STABILITY:
            continue
        out = solver.leave_out(problem.prior, weighted, innovation, group)
        summary["loo"][group] = out["weighted"] / out["count"]
        lines.append(f"  {group:14s} {out['weighted'] / out['count']:6.2f}")
    lines.append(f"the stability constraint active at the fit: {fit.final.active}")
    down = weights < 1.0
    lines.append(f"Huber's rule downweighted {int(down.sum())} of {down.size} rows" + (
        ": " + ", ".join(f"{final.groups[i]} row {i} weight {weights[i]:.2f}" for i in np.flatnonzero(down)) if down.any() else ""))
    summary["consistency"] = {"before": rec["before"], "after": rec["after"], "sizes": rec["sizes"], "downweighted": int(down.sum())}
    record(3, "the checks: consistency before and after IRIS's systematic parts, leave-one-out, the robust weighting", None, "\n".join(lines))

    # ---- 4. the solve ------------------------------------------------------------------------------------------
    _, from_mean, _ = driver.run(CONTROL, pool, start=problem.prior.mean, write=False)
    delivered = driver._delivered(fit.final.products, problem.radio)
    from_mean_T = driver._delivered(from_mean.final.products, problem.radio)
    starts = {k: float(np.nanmax(np.abs(delivered[k] - from_mean_T[k]))) for k in delivered}
    identity = {t: float(np.nanmax(np.abs(tree(p).to_dataset(inherit=False)["pressure_identity_residual"].values))) for t, p in products.items()}
    lines = [f"first fit {first.passes} passes, second fit {fit.passes} passes; their changes in delivered T (K):"]
    lines += [f"  first, pass {i + 1}: " + ", ".join(f"{k} {v:.3f}" for k, v in c.items()) for i, c in enumerate(first.changes_K)]
    lines += [f"  second, pass {i + 1}: " + ", ".join(f"{k} {v:.3f}" for k, v in c.items()) for i, c in enumerate(fit.changes_K)]
    lines.append("from the prior mean: largest difference in delivered T " + ", ".join(f"{k} {v:.3f} K" for k, v in starts.items()) + " (bound 0.1 K)")
    lines.append("pressure identity of each run: " + ", ".join(f"{k} {v:.2e}" for k, v in identity.items()))
    summary["solve"] = {"first": first.passes, "second": fit.passes, "starts": starts, "identity": identity}
    record(4, "the solve: passes, the start from the prior mean within 0.1 K, the pressure identity", max(starts.values()) < 0.1 and max(identity.values()) < 1e-2, "\n".join(lines))

    check_reference(problem, pool)
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    bounded = [r for r in results if r is not None]
    print(f"{sum(bounded)} of {len(bounded)} bounded checks pass; {results.count(None)} reported without a bound")
    return 0 if all(bounded) else 1


def check_reference(problem, pool):
    """Check 5: the reference does not matter. The fit with the Voyager 2 egress as reference against the
    Lindal-referenced one (`problem.output`), delivered at the four profiles, bound 0.5 K."""
    at_p = np.array(TABLE_MBAR) * 100.0
    problem_31, fit_31, rec_31 = driver.run(CONTROL_31S, pool)
    lines, worst = [], 0.0
    pairs = [("lindal", t) for t in H4_RUNS] + [("voyager2_egress", "lindal")]
    winds_ab = {"Lindal-referenced": problem.output, "31.2 S-referenced": problem_31.output}
    delivered_ab = {}
    for key, wind_file in winds_ab.items():
        for anchor, target in pairs:
            base = H4_RUNS[target][1] if anchor == "lindal" else "voyager2_egress_to_36n"
            name = f"reference_check_{anchor}_to_{target}"
            work = Path("fits/voyager_fit_2/work_reference") / name
            namelist = driver.copy_run(base, name, driver._read(wind_file, "wind"), work.parent, Path.cwd())
            delivered_ab[(key, anchor, target)] = tree(str(forward_transfer.run(namelist).product)).to_dataset(inherit=False)
    shutil.rmtree("fits/voyager_fit_2/work_reference", ignore_errors=True)
    for anchor, target in pairs:
        a = delivered_ab[("Lindal-referenced", anchor, target)]
        b = delivered_ab[("31.2 S-referenced", anchor, target)]
        d = b["temperature_K"].values - interp_ln_p(a["pressure_Pa"].values, a["temperature_K"].values, b["pressure_Pa"].values)
        worst = max(worst, float(np.nanmax(np.abs(d))))
        lines.append(f"  at {LABEL[target]} from the {anchor} anchor: largest |difference| {np.nanmax(np.abs(d)):.3f} K, RMS {np.sqrt(np.nanmean(d ** 2)):.3f}; at "
                     + ", ".join(f"{v:g}" for v in TABLE_MBAR) + " mbar " + " ".join(f"{v:+.2f}" for v in interp_ln_p(b["pressure_Pa"].values, d, at_p)))
    lines.insert(0, f"the 31.2 S-referenced fit: first {rec_31['first'].passes} and second {fit_31.passes} passes; ratios after "
                 + ratios(rec_31["after"]))
    summary["reference"] = {"worst": worst}
    record(5, "the reference does not matter: with the Voyager 2 egress as reference, delivered within 0.5 K at the four profiles", worst <= 0.5, "\n".join(lines))


def resume():
    with open(RESUME, "rb") as handle:
        context = pickle.load(handle)
    problem = context["problem"]
    if not hasattr(problem, "wind_out"):
        # A resume point from before the wind files reached above 1 Pa: the problem's files' grid, and the
        # Lindal-referenced wind written on it from the same state (REVIEW_21_step1, check 5 item 1).
        problem.__post_init__()
        rec = context["rec"]
        driver.write_wind(problem, context["fit"], {
            "iris_systematic_per_rad": ", ".join(f"{k} {v:.6g}" for k, v in rec["sizes"].items()),
            "iris_systematic_latitude_correlation": f"{problem.systematic_length_deg:g} deg, {problem.systematic_form}",
            "passes": f"first fit {context['first'].passes}, second fit {context['fit'].passes}"})
    with ProcessPoolExecutor(max_workers=driver.WORKERS) as pool:
        return rest(context, pool)


def refit():
    """Hypothesis 4 again under the current control, after the check 5 diagnosis: the radio profiles averaged into
    the forward wind grid's layers. Reports the change in delivered T at each target from the resume point's fit,
    bound 1 K, writes the wind file, then check 5."""
    with open(RESUME, "rb") as handle:
        context = pickle.load(handle)
    old_problem, old = context["problem"], context["fit"]
    at_p = np.array(TABLE_MBAR) * 100.0
    with ProcessPoolExecutor(max_workers=driver.WORKERS) as pool:
        problem, fit, rec = driver.run(CONTROL, pool)
        lines = [f"radio profiles averaged into the wind grid's layers: {problem.radio.average_layers}; the floor "
                 f"{problem.radio.floor_K:g} K once per layer; first fit {rec['first'].passes} and second {fit.passes} passes; "
                 "ratios after " + ratios(rec["after"]),
                 "the IRIS parts, per radian: " + ", ".join(f"{k} {v:.4g}" for k, v in rec["sizes"].items()),
                 "delivered T, the refit minus the resume point's fit (K), at " + ", ".join(f"{v:g}" for v in TABLE_MBAR) + " mbar:"]
        worst = 0.0
        for t in problem.radio.targets:
            a = old.final.products[t.name].to_dataset(inherit=False)
            b = fit.final.products[t.name].to_dataset(inherit=False)
            d = b["temperature_K"].values - interp_ln_p(a["pressure_Pa"].values, a["temperature_K"].values, b["pressure_Pa"].values)
            keep = (b["pressure_Pa"].values >= 100.0) & (b["pressure_Pa"].values <= 1.3e5)
            worst = max(worst, float(np.nanmax(np.abs(d[keep]))))
            lines.append(f"  {t.name}: " + " ".join(f"{v:+.2f}" for v in interp_ln_p(b["pressure_Pa"].values, d, at_p))
                         + f"; largest |difference| 1 mbar to 1.3 bar {np.nanmax(np.abs(d[keep])):.2f} K")
        summary["refit"] = {"worst": worst}
        record("2r", "hypothesis 4 refit with the radio profiles averaged into the wind grid's layers: its change from the resume point's fit, bound 1 K",
               worst < 1.0, "\n".join(lines))
        check_reference(problem, pool)
    (HERE / "summary_refit.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    bounded = [r for r in results if r is not None]
    print(f"{sum(bounded)} of {len(bounded)} bounded checks pass")
    return 0 if all(bounded) else 1


def problem_control(path):
    import tomllib
    with open(path, "rb") as handle:
        return tomllib.load(handle)


if __name__ == "__main__":
    raise SystemExit(refit() if "--refit" in sys.argv[1:] else resume() if "--resume" in sys.argv[1:] else main())
