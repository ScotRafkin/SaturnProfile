"""SPEC_19 v0.4: why check 1(b) is linear. The first eight truths of the linear ensemble, each marched once: does
it march, how cold and how hot it gets, and how far it lies from hypothesis 1's projection. Prints; writes no product.

A truth of this prior has so large a geopotential spread across latitude that it asks for the wind beyond the
file's grid in the columns where the top isobar stands low, so the wind grid is carried here to 1e-4 and 1e8 Pa at
the same spacing. No delivered level uses those nodes, and the state defines the wind there (the shear held at its
edge value); the grids are compared on hypothesis 1's projection first.

Run from the repository root; it takes about an hour.
"""

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import accept_step19_1 as accept  # noqa: E402
from casspian.tools.wind.estimation import driver  # noqa: E402

TRUTHS = 8


def wider_grid(problem):
    p = np.asarray(problem.wind["pressure_Pa"].values, dtype="float64")
    above, below = 10.0 ** np.arange(-4.0, 0.0, 0.05), 10.0 ** np.arange(6.05, 8.01, 0.05)
    index = np.concatenate([np.zeros(above.size, dtype=int), np.arange(p.size), np.full(below.size, p.size - 1)])
    wide = problem.wind.isel(pressure=index).load()
    wide["pressure_Pa"] = ("pressure", np.concatenate([above, p, below]), problem.wind["pressure_Pa"].attrs)
    return replace(problem, wind=wide)


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    on_file = driver.read_control(accept.CONTROL)
    problem = wider_grid(on_file)
    work = on_file.directory / "work_truths"
    start = on_file.basis.project(driver.base_shear(on_file), on_file.lat, on_file.pressure)
    with ProcessPoolExecutor(max_workers=driver.WORKERS) as pool:
        narrow = driver.linearize(on_file, start, work, pool)
        base = driver.linearize(problem, start, work, pool)
        T0, T_narrow = driver._delivered(base.products, problem.radio), driver._delivered(narrow.products, problem.radio)
        print("the wider grid against the file's, hypothesis 1's projection: largest difference in T "
              f"{max(float(np.nanmax(np.abs(T0[k] - T_narrow[k]))) for k in T0):.3e} K, h {np.max(np.abs(base.blocks.h - narrow.blocks.h)):.3e}, "
              f"H {np.max(np.abs(base.blocks.H - narrow.blocks.H)):.3e}", flush=True)
        print("hypothesis 1, coldest T (K):", {k: round(float(np.nanmin(v)), 1) for k, v in T0.items()}, "rows", base.blocks.h.size, flush=True)
        for i in range(TRUTHS):
            truth = accept.draw_from_prior(problem, np.random.default_rng([accept.SEED, i]))
            try:
                at = driver.linearize(problem, truth, work, pool)
            except ValueError as error:
                print(i, "refused:", str(error)[:160], flush=True)
                continue
            T = driver._delivered(at.products, problem.radio)
            line = {k: (round(float(np.nanmin(T[k])), 1), round(float(np.nanmax(T[k])), 1),
                        round(float(np.nanmax(np.abs(T[k] - T0[k]))), 1) if T[k].shape == T0[k].shape else None) for k in T}
            print(i, "(coldest, hottest, largest change from hypothesis 1), K:", line, "rows", at.blocks.h.size, flush=True)
