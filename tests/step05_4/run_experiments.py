"""SPEC_05 Step 4: build and run every named experiment, at both geopotential spacings.

The run of record is each run's own directory under `forward/`, built with `casspian-run-inputs`
and run with `casspian-forward`'s transfer driver at the namelist's 5e4 m2/s2, figures F5 to F9
included. The 2.5e4 run is made in a copy, `reports/step05_4/spacing_2p5e4/<run>/`, with the
inputs of the run of record, the namelist's `[grid] geopotential_spacing_m2s2` halved, the anchor
path rewritten for the copy's depth (as `step04_5` check 11 reruns the closure namelist) and
figures off; run 7f is made at 5e4 only. No committed file is edited.

**The relaxation, named here.** Until Step 4 is committed the working tree is not clean, every
input built here carries `-dirty`, and `forward` refuses a `-dirty` input. The one refusal point,
`casspian.lib.control._refuse_dirty_commit`, is relaxed in this process only; each product keeps
the stamp it earns, which `results.json` records. Run again after the acceptance commit and the
products carry the clean commit.

A run that stops is a result (SPEC_05 Step 4, run 7): its named failure is recorded, with its type
and message, and nothing is tuned. Writes `reports/step05_4/results.json`. Run from the repository
root.
"""

import json
import shutil
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.tools.run import run_inputs

HERE = Path("reports/step05_4")
HALF = HERE / "spacing_2p5e4"
RUNS = ("shear_r2_uniform", "shear_r3a_nowind", "shear_r3b_nowind_anchor", "shear_r3c_half",
        "shear_r4_decay12", "shear_r5_decay20", "shear_r6_decay40", "shear_r7_decay_linp",
        "shear_r7f_decay_linp_fine", "shear_r8_increase25", "shear_r9_increase50")
ONE_SPACING = ("shear_r7f_decay_linp_fine",)


def summary(product: Path) -> dict:
    """The scalars a run is reported by, read back from its product."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(product, engine="netcdf4").load()
    tree.close()
    root = tree.to_dataset(inherit=False)
    record = tree["transfer_record"].attrs
    return {
        "product": product.as_posix(),
        "casspian_git_commit": str(root.attrs.get("casspian_git_commit", "")),
        "pressure_identity_largest_abs": float(record["pressure_identity_largest_abs"]),
        "pressure_identity_rms": float(record["pressure_identity_rms"]),
        "pressure_identity_largest_level": int(record["pressure_identity_largest_level"]),
        "outer_loop_passes": int(record["outer_loop_passes"]),
        "outer_loop_residuals_ln_p": str(record["outer_loop_residuals_ln_p"]),
        "reference_surface_radius_m": float(root["reference_surface_radius_m"].values),
        "target_latitude_deg": float(root["latitude_planetocentric_deg"].values),
    }


def attempt(label, action, read=True):
    """Run `action`, timed; a failure is recorded with its type and message, not raised."""
    started = time.perf_counter()
    try:
        written = Path(action())
        out = {"status": "completed", **(summary(written) if read else {"last_written": written.as_posix()})}
    except Exception as exc:  # a named failure is a result (SPEC_05 Step 4)
        out = {"status": "failed", "failure_type": type(exc).__name__, "failure": str(exc),
               "traceback_tail": traceback.format_exc().strip().splitlines()[-3:]}
    out["wall_s"] = round(time.perf_counter() - started, 1)
    if out["status"] != "completed":
        detail = f", {out['failure_type']}: {out['failure'][:160]}"
    elif read:
        detail = (f", identity {out['pressure_identity_largest_abs']:.3e}, "
                  f"passes {out['outer_loop_passes']}")
    else:
        detail = ""
    print(f"{label}: {out['status']} in {out['wall_s']} s{detail}")
    return out


def half_spacing(run: str) -> Path:
    """The 2.5e4 copy of a run: the namelist with its spacing halved, the record's inputs."""
    source = Path("forward") / run
    copy = HALF / run
    if copy.exists():
        shutil.rmtree(copy)
    (copy / "inputs").mkdir(parents=True)
    for path in (source / "inputs").glob("*.nc"):
        shutil.copy2(path, copy / "inputs" / path.name)
    text = (source / f"{run}.toml").read_text(encoding="utf-8")
    for old, new in (('"../../occul_data/', '"../../../../occul_data/'),
                     ("geopotential_spacing_m2s2 = 5.0e4", "geopotential_spacing_m2s2 = 2.5e4"),
                     ("figures = true", "figures = false")):
        if text.count(old) != 1:
            raise RuntimeError(f"{run}: the namelist edit {old!r} matched {text.count(old)} times")
        text = text.replace(old, new)
    (copy / f"{run}.toml").write_bytes(text.encode("utf-8"))
    return copy / f"{run}.toml"


def main(runs=RUNS) -> dict:
    HERE.mkdir(parents=True, exist_ok=True)
    original = ctl._refuse_dirty_commit
    ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
        attrs.get("casspian_git_commit", ""))
    print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this process only\n")
    results = {}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for run in runs:
                directory = Path("forward") / run
                built = attempt(f"{run} inputs",
                                lambda: run_inputs.build(directory / f"{run}_build.toml")[-1],
                                read=False)
                entry = {"inputs": built}
                if built["status"] == "completed":
                    entry["5e4"] = attempt(f"{run} at 5e4", lambda: forward_transfer.run(
                        directory / f"{run}.toml").product)
                    if run not in ONE_SPACING:
                        entry["2.5e4"] = attempt(f"{run} at 2.5e4", lambda: forward_transfer.run(
                            half_spacing(run)).product)
                results[run] = entry
    finally:
        ctl._refuse_dirty_commit = original
    # Runs named on the command line replace their own entries; the others are kept.
    path = HERE / "results.json"
    merged = json.loads(path.read_text(encoding="utf-8")) if path.exists() and runs != RUNS else {}
    merged.update(results)
    merged = {run: merged[run] for run in RUNS if run in merged}
    path.write_text(json.dumps(merged, indent=1, default=float), encoding="utf-8")
    return merged


if __name__ == "__main__":
    chosen = tuple(sys.argv[1:]) or RUNS
    main(chosen)
