"""SPEC_20 section 2 deliverable 4: the anchor matrix under hypothesis 3. Each of the three new anchors
transferred to the latitude of each of the other three anchors. The Lindal anchor's three are SPEC_19's
runs, `lindal_voyager_fit_1`, `_74s` and `_03s`, read as they are.

Each run is made by `casspian-new-run` from `lindal_voyager_fit_1`, whose inputs are hypothesis 3's, with
only the anchor and the target changed; its inputs are built by `casspian-run-inputs` and it runs by
`casspian-forward`. A target is the other anchor's planetocentric latitude as its kind N file states it,
which is where SPEC_19's runs went for the Voyager anchors.

Run from the repository root; about half an hour.
"""

import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import xarray as xr

from casspian.tools.run import new_run

BASE = Path("forward/lindal_voyager_fit_1")
#: Anchor slug, and its short name in a run's name.
ANCHORS = {"lindal": "36n", "voyager2_egress": "31s", "voyager1_ingress": "74s", "voyager1_egress": "03s"}
#: SPEC_19's runs, the Lindal anchor's row of the matrix.
LINDAL_RUNS = {"voyager2_egress": "lindal_voyager_fit_1", "voyager1_ingress": "lindal_voyager_fit_1_74s",
               "voyager1_egress": "lindal_voyager_fit_1_03s"}


def latitude(slug) -> float:
    root = xr.open_datatree(f"occul_data/lindal/{slug}_refractivity.nc").to_dataset(inherit=False)
    return float(root["latitude_planetocentric_deg"])


def name(anchor, target, suffix="") -> str:
    return f"{anchor}_to_{ANCHORS[target]}{suffix}"


def runs(lindal_runs=None, suffix="") -> dict:
    """`{(anchor, target): run name}` for every pair of the matrix. SPEC_21 passes hypothesis 4's Lindal runs and
    a suffix for the others' names."""
    lindal_runs = lindal_runs or LINDAL_RUNS
    out = {}
    for anchor in ANCHORS:
        for target in ANCHORS:
            if anchor != target:
                out[(anchor, target)] = lindal_runs[target] if anchor == "lindal" else name(anchor, target, suffix)
    return out


def make(pair, base=BASE, suffix="", hypothesis="hypothesis 3"):
    anchor, target = pair
    run = name(anchor, target, suffix)
    directory = Path("forward") / run
    if directory.exists():
        shutil.rmtree(directory)
    new_run.new_run(run, base)
    namelist = directory / f"{run}.toml"
    text = namelist.read_text(encoding="utf-8")
    text = re.sub(r'(?m)^slug\s*=.*$', f'slug   = "{anchor}"', text)
    text = re.sub(r'(?m)^path\s*=.*$', f'path   = "../../occul_data/lindal/{anchor}_refractivity.nc"', text)
    text = re.sub(r'(?m)^latitude_planetocentric_deg\s*=.*$', f"latitude_planetocentric_deg = {latitude(target)!r}", text)
    if anchor == "voyager1_ingress":
        # The trace stops at 1000 mbar, so 1 bar is not a produced level; the datum places only the altitude.
        text = re.sub(r'(?m)^datum_isobar_Pa\s*=.*$',
                      "datum_isobar_Pa = 1.0e4   # SPEC_20: the anchor reaches 1000 mbar only; the datum places only the altitude", text)
    text = text.replace("# CASSPIAN run namelist, transfer mode",
                        f"# SPEC_20 anchor matrix: the {anchor} anchor to the {target} latitude under hypothesis 3, made\n"
                        f"# from lindal_voyager_fit_1 by casspian-new-run with only the anchor and the target changed.\n"
                        "# CASSPIAN run namelist, transfer mode", 1)
    namelist.write_bytes(text.encode("utf-8"))
    # The run's composition on the anchor's levels, which the transfer places the anchor on: the same
    # mixture, 94/6 with Lindal's ammonia (REVIEW_20_step0 ruling 2), from the anchor's own raw bundle.
    build = directory / f"{run}_build.toml"
    text = build.read_text(encoding="utf-8")
    head, tail = text.split("[composition]", 1)
    tail = tail.replace('raw_bundle     = "../../occul_data/lindal/raw/lindal_raw.nc"',
                        f'raw_bundle     = "../../occul_data/lindal/raw/{anchor}_raw.nc"', 1)
    build.write_bytes((head + "[composition]" + tail).encode("utf-8"))
    subprocess.run(["casspian-run-inputs", str(directory / f"{run}_build.toml")], check=True, capture_output=True)
    done = subprocess.run(["casspian-forward", str(namelist)], capture_output=True, text=True)
    if done.returncode:
        raise RuntimeError(f"{run}: {done.stderr[-2000:]}")
    return run


if __name__ == "__main__":
    pairs = [p for p, r in runs().items() if p[0] != "lindal"
             and not (Path("forward") / r / "output" / f"{r}_profile.nc").exists()]
    with ProcessPoolExecutor(max_workers=int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool:
        for run in pool.map(make, pairs):
            print("ran", run, flush=True)
