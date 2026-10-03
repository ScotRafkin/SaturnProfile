"""SPEC_06 v0.4 check 7: what sets the sign and size of the identity at the stop pressure.

Run 5's case (`decay_above`, 20 percent per scale height, stop fraction 0) at 5e4, with the stop
pressure placed at five positions `s` across one anchor layer, the layer from level 15 (632.3 Pa)
to level 16 (795.4 Pa): `p_stop = p_15 (p_16 / p_15)^s`, `s` in 0.02, 0.25, 0.5, 0.75, 0.98; and
run 6's case (40 percent, stop at 8000 Pa) the same way across the layer from level 27 (7243 Pa) to
level 28 (8728 Pa). The
stop pressure is added to the wind file's pressure grid as a node, so the wind file's change of
shear sits exactly at it (between nodes the file is linear in `ln p` and the ramp is linear in
`ln p`, so the rest of the decay is represented exactly). Each run is a copy under
`reports/step06_1/stop_position/`, the `-dirty` refusal relaxed in this process only.

Reported for each: the absolute pressure offset `p - p_label` (Pa) at levels 14 to 18, its
increment across each layer, and the delivered temperature change from run 2 at those levels.
Writes `reports/step06_1/stop_position.json`. Run from the repository root, after
`tests/step06_1/accept_step06_1.py` (it reads run 2's product from `reports/step05_4/results.json`).
"""

import itertools
import json
import re
import shutil
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.tools.run import run_inputs

HERE = Path("reports/step06_1/stop_position")
POSITIONS = (0.02, 0.25, 0.5, 0.75, 0.98)
#: (run, the stop pressure as its build file writes it, the layer's upper level)
CASES = (("shear_r5_decay20", "700.0", 15), ("shear_r6_decay40", "8000.0", 27))


def root(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(path, engine="netcdf4").load()
    tree.close()
    return tree.to_dataset(inherit=False)


def grid_of(text):
    block = re.search(r"pressure_grid_Pa = \[(.*?)\]", text, flags=re.S)
    return block, [float(v) for v in block.group(1).replace("\n", " ").split(",") if v.strip()]


def copy_with_stop(SOURCE, stop_text, p_stop, tag):
    RUN = SOURCE.name
    copy = HERE / tag
    if copy.exists():
        shutil.rmtree(copy)
    copy.mkdir(parents=True)
    build = (SOURCE / f"{RUN}_build.toml").read_text(encoding="utf-8").replace('"../../', '"../../../../')
    block, grid = grid_of(build)
    grid = sorted(set(grid) | {p_stop})
    build = build.replace(block.group(0), "pressure_grid_Pa = [" + ", ".join(f"{v!r}" for v in grid) + "]")
    old = f"stop_pressure_Pa = {stop_text}"
    if build.count(old) != 1:
        raise RuntimeError(f"{old!r} matched {build.count(old)} times")
    build = build.replace(old, f"stop_pressure_Pa = {p_stop!r}")
    (copy / f"{RUN}_build.toml").write_bytes(build.encode("utf-8"))
    namelist = (SOURCE / f"{RUN}.toml").read_text(encoding="utf-8").replace('"../../', '"../../../../')
    namelist = namelist.replace("figures = true", "figures = false")
    (copy / f"{RUN}.toml").write_bytes(namelist.encode("utf-8"))
    run_inputs.build(copy / f"{RUN}_build.toml")
    return copy / f"{RUN}.toml"


results = json.loads(Path("reports/step05_4/results.json").read_text(encoding="utf-8"))
run2 = root(results["shear_r2_uniform"]["5e4"]["product"])
labels = np.asarray(run2["pressure_label_Pa"].values, dtype="float64")

original = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
    attrs.get("casspian_git_commit", ""))
print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this process only\n")
out = {}
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for (run, stop_text, top), s in itertools.product(CASES, POSITIONS):
            SOURCE, LEVELS = Path("forward") / run, range(top - 2, top + 4)
            p_stop = float(labels[top] * (labels[top + 1] / labels[top]) ** s)
            product = forward_transfer.run(copy_with_stop(
                SOURCE, stop_text, p_stop, f"{run}_s{s:.2f}")).product
            d = root(product)
            p, lab = d["pressure_Pa"].values, d["pressure_label_Pa"].values
            offset = p - lab
            dT = d["temperature_K"].values - run2["temperature_K"].values
            out[f"{run} {s:.2f}"] = {
                "p_stop_Pa": p_stop,
                "offset_Pa": {int(k): float(offset[k]) for k in LEVELS},
                "dT_K": {int(k): float(dT[k]) for k in LEVELS},
                "identity_largest": float(np.max(np.abs(d["pressure_identity_residual"].values))),
                "identity_level": int(np.argmax(np.abs(d["pressure_identity_residual"].values))),
            }
            row = out[f"{run} {s:.2f}"]
            print(f"{run} s = {s:.2f}, p_stop {p_stop:.2f} Pa: offset p - p_label (Pa) at levels "
                  + ", ".join(f"{k}: {row['offset_Pa'][k]:+.2f}" for k in LEVELS)
                  + "; layer increments " + ", ".join(
                      f"{k}-{k + 1}: {row['offset_Pa'][k + 1] - row['offset_Pa'][k]:+.2f}"
                      for k in list(LEVELS)[:-1])
                  + "; dT " + ", ".join(f"{k}: {row['dT_K'][k]:+.2f}" for k in LEVELS)
                  + f"; identity {row['identity_largest']:.2e} at level {row['identity_level']}",
                  flush=True)
finally:
    ctl._refuse_dirty_commit = original
(Path("reports/step06_1") / "stop_position.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
