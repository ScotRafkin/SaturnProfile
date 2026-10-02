"""SPEC_05 Step 4, the two additions of REVIEW_05_step4.

1. The outer loop's residual, pass by pass: reruns the two runs that stopped (run 6 at 5e4, run 7
   at 2.5e4); `forward/transfer.py` now prints every pass's largest |d ln p|, and the named
   failure's message carries the whole history.
2. The kink at `p_s`: run 8's case with `shear_reference_pressure_Pa` at 1.05e5 Pa (the review's
   run) and at 1.3e5 Pa (past one mesh cell and the wind grid's 1.259e5 node, at the author's
   request), at 5e4, in copies under `reports/step05_4/r8_ps105/` and `r8_ps130/` (no new run
   directory); reports the delivered temperature change from run 2 at the 998.7 mbar level and
   below, beside run 8's own.

The `-dirty` refusal is relaxed in this process only, as in `run_experiments.py`. Writes
`reports/step05_4/review_additions.json`. Run from the repository root, after
`run_experiments.py` (it reads run 2's and run 8's products).
"""

import json
import shutil
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.tools.run import run_inputs

import run_experiments as rx

HERE = Path("reports/step05_4")
SOURCE = Path("forward/shear_r8_increase25")


def kink_copy(p_s: str) -> Path:
    """Run 8's control files under `reports/step05_4/r8_ps<tag>/`, with p_s at `p_s` Pa, paths for
    the copy's depth, figures off."""
    copy = HERE / f"r8_ps{round(float(p_s) / 1e3)}"
    if copy.exists():
        shutil.rmtree(copy)
    copy.mkdir(parents=True)
    edits = {"_build.toml": (('"../../', '"../../../', 8),
                             ("shear_reference_pressure_Pa = 1.0e5",
                              f"shear_reference_pressure_Pa = {p_s}", 1)),
             ".toml": (('"../../', '"../../../', 1), ("figures = true", "figures = false", 1))}
    for suffix, pairs in edits.items():
        name = f"{SOURCE.name}{suffix}"
        text = (SOURCE / name).read_text(encoding="utf-8")
        for old, new, count in pairs:
            if text.count(old) != count:
                raise RuntimeError(f"{name}: {old!r} matched {text.count(old)} times, not {count}")
            text = text.replace(old, new)
        (copy / name).write_bytes(text.encode("utf-8"))
    run_inputs.build(copy / f"{SOURCE.name}_build.toml")
    return forward_transfer.run(copy / f"{SOURCE.name}.toml").product


def delivered(product):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(product, engine="netcdf4").load()
    tree.close()
    return tree.to_dataset(inherit=False)


def main():
    original = ctl._refuse_dirty_commit
    ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
        attrs.get("casspian_git_commit", ""))
    print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this process only\n")
    out = {}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out["shear_r6_decay40 at 5e4"] = rx.attempt("shear_r6_decay40 at 5e4", lambda: forward_transfer.run(
                Path("forward/shear_r6_decay40/shear_r6_decay40.toml")).product)
            out["shear_r7_decay_linp at 2.5e4"] = rx.attempt("shear_r7_decay_linp at 2.5e4", lambda: forward_transfer.run(
                rx.half_spacing("shear_r7_decay_linp")).product)
            for p_s in ("1.05e5", "1.3e5"):
                out[f"r8 p_s {p_s} at 5e4"] = rx.attempt(f"r8 p_s {p_s} at 5e4", lambda: kink_copy(p_s))
    finally:
        ctl._refuse_dirty_commit = original

    runs = json.loads((HERE / "results.json").read_text(encoding="utf-8"))
    run2 = delivered(runs["shear_r2_uniform"]["5e4"]["product"])
    labels = run2["pressure_label_Pa"].values
    k = int(np.argmin(np.abs(labels - 1.0e5)))
    T2 = run2["temperature_K"].values
    kink = {"level": k, "label_Pa": float(labels[k])}
    for name, product in (("run 8, p_s 1e5", runs["shear_r8_increase25"]["5e4"]["product"]),
                          ("run 8, p_s 1.05e5", out["r8 p_s 1.05e5 at 5e4"].get("product")),
                          ("run 8, p_s 1.3e5", out["r8 p_s 1.3e5 at 5e4"].get("product"))):
        if product is None:
            continue
        dT = delivered(product)["temperature_K"].values - T2
        kink[name] = {"dT_at_level_K": float(dT[k]),
                      "dT_below": [round(float(v), 3) for v in dT[k + 1:]]}
    out["kink"] = kink
    print(json.dumps(kink, indent=1))
    (HERE / "review_additions.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")


if __name__ == "__main__":
    main()
