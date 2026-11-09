"""Hypothesis 1's projection as the truth, in the linear system at that truth: the noise-free solve and 20 noise draws.
Does the truth's own structure give IRIS 110's ratio and the jet's overshoot, or was it 1(a)'s draw? And the operator's
inexactness: actual h at the solve's state against the linear prediction, weighted by R. Prints; writes no product."""
import shutil
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, "tests/step19_1")

if __name__ == "__main__":
    warnings.simplefilter("ignore")
    import accept_step19_1 as m
    from casspian.tools.wind.estimation import driver as d, solver
    problem = d.read_control(m.CONTROL)
    truth = problem.basis.project(d.base_shear(problem), problem.lat, problem.pressure)
    work = problem.directory / "work_linear"
    lat = np.array([15.0, 10.0, 5.0, 2.5, 0.0, -2.5, -5.0, -10.0])
    jet = lambda x: problem.basis.shear(x, lat, np.array([11000.0]))[0]
    names = None
    with ProcessPoolExecutor(max_workers=d.WORKERS) as pool:
        at = d.linearize(problem, truth, work, pool)
        b = at.blocks
        names = list(dict.fromkeys(b.groups))
        L = np.linalg.cholesky(b.R)
        ratios = lambda bl, s: {k: v["ratio"] for k, v in solver.consistency(bl, s).items()}
        clean = solver.Blocks(b.h, b.h, b.H, b.R, b.groups)
        s0 = solver.linearized(problem.prior, truth, clean)
        print("no noise, ratios", {k: round(v, 3) for k, v in ratios(clean, s0).items()}, flush=True)
        print("latitudes", lat.tolist(), flush=True)
        print("truth   ", np.round(jet(truth), 1).tolist(), flush=True)
        print("no noise", np.round(jet(s0.state), 1).tolist(), flush=True)

        rng = np.random.default_rng(m.SEED)
        first = solver.Blocks(b.h + L @ rng.standard_normal(b.h.size), b.h, b.H, b.R, b.groups)
        s1 = solver.linearized(problem.prior, truth, first)
        print("1(a)'s draw, ratios", {k: round(v, 3) for k, v in ratios(first, s1).items()}, flush=True)
        print("1(a)'s draw, jet", np.round(jet(s1.state), 1).tolist(), flush=True)

        rows, jets = [], []
        for i in range(20):
            rng = np.random.default_rng([m.SEED, 2000 + i])
            drawn = solver.Blocks(b.h + L @ rng.standard_normal(b.h.size), b.h, b.H, b.R, b.groups)
            s = solver.linearized(problem.prior, truth, drawn)
            rows.append(ratios(drawn, s))
            jets.append(jet(s.state))
        keys = list(rows[0])
        print("20 draws at this truth, mean ratio (spread):", {k: (round(float(np.mean([r[k] for r in rows])), 2), round(float(np.std([r[k] for r in rows], ddof=1)), 2)) for k in keys}, flush=True)
        print("IRIS 110 ratio, each draw", [round(r["IRIS 110 mbar"], 2) for r in rows], flush=True)
        jets = np.array(jets)
        print("jet, mean of 20", np.round(jets.mean(axis=0), 1).tolist(), flush=True)
        print("jet, spread of 20", np.round(jets.std(axis=0, ddof=1), 1).tolist(), flush=True)

        # The operator: march at 1(a)'s one-solve state, and compare the actual h with the linear prediction there.
        delta = s1.state - truth
        there = d.linearize(problem, s1.state, work, pool)
        print("rows at truth", b.h.size, "at the solve's state", there.blocks.h.size, flush=True)
        for name in names:
            i0, i1 = b.groups == name, there.blocks.groups == name
            if i0.sum() != i1.sum():
                print(name, "rows differ", int(i0.sum()), int(i1.sum()), flush=True)
                continue
            actual = there.blocks.h[i1]
            data = first.y[i0]
            predicted = b.h[i0] + (b.H @ delta)[i0]
            Ri = np.linalg.inv(b.R[np.ix_(i0, i0)])
            eps = actual - predicted
            lin, non = data - predicted, data - actual
            w, v = np.linalg.eigh(b.R[np.ix_(i0, i0)])
            order = np.argsort(w)[::-1]
            modes = v[:, order].T @ eps
            scale = np.exp(data) if name[0].isdigit() else 1.0
            tail = float(np.sum(modes[20:] ** 2 / w[order][20:]))
            print(f"{name}: misfit with the linear prediction {lin @ Ri @ lin:.1f}, with the actual h {non @ Ri @ non:.1f}; "
                  f"operator error alone {eps @ Ri @ eps:.1f}, of it modes 21 on {tail:.1f}; RMS error {float(np.sqrt(np.mean((eps * scale) ** 2))):.4f}, "
                  f"smallest mode sd {float(np.sqrt(w.min()) * (np.mean(scale))):.4f}", flush=True)
    shutil.rmtree(work, ignore_errors=True)
