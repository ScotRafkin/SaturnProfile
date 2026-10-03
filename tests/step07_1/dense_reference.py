"""SPEC_07 v0.4 §8 ruling 5: the dense-node reference, before the author decides between L2 and L.

Runs 5, 7 and 7f rebuilt with the same hypothesis on ten times the pressure nodes (the source wind
built on the denser grid, the shear case unchanged: 100 nodes per decade for runs 5 and 7, 200 for
run 7f, every coarse node among them), at 5e4, figures off, in copies under
`reports/step07_1/dense/<reading>/<run>/`.

`--reading L2` runs them with the working tree (decision L2); `--reading L` with the code that
imports, which the caller points at a worktree at `a238306` (decision L) by putting its `src` first on
`PYTHONPATH`. Either way the `-dirty` refusal is relaxed in this process only. `--compare` then
reads both readings' dense products beside the coarse ones (coarse L: SPEC_06's, set aside under
`reports/step07_1/spec06/`; coarse L2: the experiments as rerun by this step's acceptance) and writes
`reports/step07_1/dense/comparison.json` and `comparison.txt`. Run from the repository root.

Every comparison is of the delivered temperature change from run 2 under the same reading and grid.
Run 2's wind has no vertical structure, so its pressure grid does not enter its answer; the coarse
run 2 of each reading serves the dense runs of that reading.
"""

import json
import re
import shutil
import sys
import time
import warnings
from pathlib import Path

import numpy as np

HERE = Path("reports/step07_1/dense")
RUNS = {"shear_r5_decay20": 100, "shear_r7_decay_linp": 100, "shear_r7f_decay_linp_fine": 200}
#: The slope changes of each hypothesis, in Pa: `p_s` and, for run 5, the stop pressure.
CHANGES = {"shear_r5_decay20": (1.0e5, 700.0), "shear_r7_decay_linp": (1.0e5,),
           "shear_r7f_decay_linp_fine": (1.0e5,)}


def dense_grid(per_decade):
    """`10^(k / per_decade)` from 1 Pa to 1 MPa; every power of ten exact, so 1e5 is a node."""
    k = np.arange(6 * per_decade + 1)
    grid = 10.0 ** (k / per_decade)
    grid[k % per_decade == 0] = 10.0 ** (k[k % per_decade == 0] // per_decade)
    return grid


def build_and_run(reading):
    from casspian.forward import transfer as forward_transfer
    from casspian.lib import control as ctl
    from casspian.tools.run import run_inputs

    original = ctl._refuse_dirty_commit
    ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
        attrs.get("casspian_git_commit", ""))
    print(f"RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this process only; "
          f"reading {reading}, casspian from {Path(ctl.__file__).parents[1]}", flush=True)
    out = {}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for run, per_decade in RUNS.items():
                source = Path("forward") / run
                copy = HERE / reading / run
                if copy.exists():
                    shutil.rmtree(copy)
                copy.mkdir(parents=True)
                build = (source / f"{run}_build.toml").read_text(encoding="utf-8")
                build = build.replace('"../../', '"../../../../../')
                block = re.search(r"pressure_grid_Pa = \[(.*?)\]", build, flags=re.S)
                coarse = np.asarray([float(v) for v in block.group(1).replace(chr(10), " ").split(",")
                                     if v.strip()])
                grid = dense_grid(per_decade)
                missing = [float(p) for p in coarse if not np.any(np.isclose(grid, p, rtol=1e-5))]
                if missing:
                    raise RuntimeError(f"{run}: coarse nodes {missing} are not in the dense grid")
                build = build.replace(block.group(0), "pressure_grid_Pa = ["
                                      + ", ".join(repr(float(v)) for v in grid) + "]")
                (copy / f"{run}_build.toml").write_bytes(build.encode("utf-8"))
                namelist = (source / f"{run}.toml").read_text(encoding="utf-8")
                namelist = namelist.replace('"../../', '"../../../../../').replace(
                    "figures = true", "figures = false")
                (copy / f"{run}.toml").write_bytes(namelist.encode("utf-8"))
                started = time.time()
                run_inputs.build(copy / f"{run}_build.toml")
                product = forward_transfer.run(copy / f"{run}.toml").product
                out[run] = {"product": Path(product).as_posix(), "nodes": int(grid.size),
                            "wall_s": round(time.time() - started, 1)}
                print(f"{reading} {run}: {grid.size} pressure nodes, {out[run]['wall_s']} s", flush=True)
    finally:
        ctl._refuse_dirty_commit = original
    (HERE / reading / "runs.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def compare():
    import xarray as xr

    def root(path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tree = xr.open_datatree(path, engine="netcdf4").load()
        tree.close()
        return tree

    def delivered(path):
        tree = root(path)
        d = tree.to_dataset(inherit=False)
        return (np.asarray(d["temperature_K"].values, dtype="float64"),
                np.asarray(d["pressure_label_Pa"].values, dtype="float64"),
                float(tree["transfer_record"].attrs["pressure_identity_largest_abs"]),
                int(tree["transfer_record"].attrs["pressure_identity_largest_level"]))

    now = json.loads(Path("reports/step05_4/results.json").read_text(encoding="utf-8"))
    spec06 = Path("reports/step07_1/spec06")
    dense = {r: json.loads((HERE / r / "runs.json").read_text(encoding="utf-8")) for r in ("L", "L2")}
    run2 = {"L": delivered(spec06 / "shear_r2_uniform_profile.nc")[0],
            "L2": delivered(now["shear_r2_uniform"]["5e4"]["product"])[0]}
    report, lines = {}, []
    for run in RUNS:
        sets = {"coarse L": (delivered(spec06 / f"{run}_profile.nc"), "L"),
                "coarse L2": (delivered(now[run]["5e4"]["product"]), "L2"),
                "dense L": (delivered(dense["L"][run]["product"]), "L"),
                "dense L2": (delivered(dense["L2"][run]["product"]), "L2")}
        change = {name: d[0] - run2[reading] for name, (d, reading) in sets.items()}
        labels = sets["dense L2"][0][1]
        # The corner intervals: the coarse wind intervals on either side of the coarse nodes that
        # bracket each slope change (one interval either side of a change at a node).
        coarse_grid = np.sort(np.asarray(xr.open_datatree(
            Path("forward") / run / "inputs" / f"{run}_wind.nc", engine="netcdf4").to_dataset(
            inherit=False)["pressure_Pa"].values, dtype="float64"))
        corner = np.zeros(labels.shape, dtype=bool)
        for p in CHANGES[run]:
            i = int(np.searchsorted(coarse_grid, p, side="right") - 1)
            on_node = np.isclose(coarse_grid[i], p, rtol=1e-9)
            low = coarse_grid[max(i - 1, 0)]
            high = coarse_grid[min(i + 1, coarse_grid.size - 1) if on_node else min(i + 2, coarse_grid.size - 1)]
            corner |= (labels >= low) & (labels <= high)
        reference = change["dense L2"]
        rows = {}
        for name in ("coarse L", "coarse L2", "dense L"):
            diff = change[name] - reference
            rows[name] = {"away_largest_K": float(np.max(np.abs(diff[~corner]))),
                          "away_mean_K": float(np.mean(np.abs(diff[~corner]))),
                          "inside_K": {int(k): float(diff[k]) for k in np.flatnonzero(corner)}}
        identities = {name: {"largest": d[2], "level": d[3]} for name, (d, _) in sets.items()}
        report[run] = {"corner_levels": [int(k) for k in np.flatnonzero(corner)], "against_dense_L2": rows,
                       "identities": identities,
                       "dense_change_inside_K": {int(k): float(reference[k]) for k in np.flatnonzero(corner)}}
        lines.append(f"{run}: corner intervals hold levels {[int(k) for k in np.flatnonzero(corner)]}")
        for name, row in rows.items():
            lines.append(f"    {name} minus dense L2, away from the corners: largest {row['away_largest_K']:.3f} K, "
                         f"mean {row['away_mean_K']:.3f} K")
        lines.append("    inside the corners (level: dense L2's change; coarse L, coarse L2, dense L minus it), K:")
        for k in np.flatnonzero(corner):
            lines.append(f"        {k} ({labels[k] / 100:.1f} hPa): {reference[k]:+.2f}; "
                         + ", ".join(f"{rows[n]['inside_K'][int(k)]:+.2f}" for n in ("coarse L", "coarse L2", "dense L")))
        lines.append("    identities: " + "; ".join(f"{n} {v['largest']:.2e} (level {v['level']})"
                                                     for n, v in identities.items()))
    (HERE / "comparison.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    (HERE / "comparison.txt").write_text(chr(10).join(lines) + chr(10), encoding="utf-8")
    print(chr(10).join(lines))


if __name__ == "__main__":
    if "--compare" in sys.argv:
        compare()
    else:
        build_and_run(sys.argv[sys.argv.index("--reading") + 1])
