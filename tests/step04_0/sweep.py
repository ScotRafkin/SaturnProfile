# The full-chain sweep: every registered product rebuilt clean; run after the acceptance commit of a step that changes one.
"""The Step 0 sweep: every product rebuilt on the clean tree at the acceptance commit.

Handoff of 17 September section 8. The tree is clean, so nothing is relaxed: `refrac` and
`forward` accept their inputs as they stand and every product carries the acceptance commit
rather than `-dirty`. The script rebuilds the chain with the SPEC_01 tools and
`casspian-refrac`, the closure run's inputs with `casspian-run-inputs` and its product with
`casspian-forward`, and the transfer run's inputs; then reads `casspian_git_commit` back from
every file and refuses if any ends in `-dirty`; then prints the SHA-256 of each in the row
format the `step02_1` suite reads.

Run from the repository root on a tree whose porcelain output is empty.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import xarray as xr

ROOT = Path(__file__).resolve().parents[2]

from casspian.lib import io as cio  # noqa: E402
from casspian.refrac import product as refrac_product  # noqa: E402
from casspian.forward import production as forward_production  # noqa: E402
from casspian.tools.composition import build_composition  # noqa: E402
from casspian.tools.gravity import build_gravity, build_rotation  # noqa: E402
from casspian.tools.lindal import build_inputs, build_raw  # noqa: E402
from casspian.tools.run import run_inputs  # noqa: E402
from casspian.tools.wind import build_wind  # noqa: E402

LINDAL = ROOT / "occul_data" / "lindal"
CLOSURE = ROOT / "forward" / "lindal_closure"
TRANSFER = ROOT / "forward" / "lindal_transfer"

KIND_OF = {
    "lindal_raw.nc": "raw",
    "lindal_thermo.nc": "thermo",
    "lindal_geodesy.nc": "geodesy",
    "lindal_reduction.toml": "manifest",
    "lindal_refractivity.nc": "refractivity",
    "lindal_closure_profile.nc": "profile",
    "lindal_transfer_profile.nc": "profile",
}


def kind_of(path: Path) -> str:
    if path.name in KIND_OF:
        return KIND_OF[path.name]
    for suffix, kind in (("_gravity.nc", "gravity"), ("_rotation.nc", "rotation"),
                         ("_wind.nc", "wind"), ("_composition.nc", "composition")):
        if path.name.endswith(suffix):
            return kind
    raise RuntimeError(f"no kind for {path.name}")


def registered() -> list[str]:
    """The committed netCDF products the sweep rewrites (SPEC_05 decision D1), fixtures excluded."""
    listed = subprocess.run(["git", "ls-files", "*.nc"], cwd=ROOT, capture_output=True,
                            text=True, check=True).stdout.split()
    return [name for name in listed if not name.startswith("tests/")]


def main() -> int:
    porcelain = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                               capture_output=True, text=True).stdout.strip()
    if porcelain:
        print("REFUSED: the tree is not clean; the sweep needs an empty porcelain output.")
        print(porcelain)
        return 1
    print(f"clean tree at {cio.git_commit()}")
    # The registered products are committed (D1), so rewriting the first one would make the tree
    # dirty and stamp every later file -dirty. They are marked skip-worktree for the sweep only, so
    # the stamp checks every code and control file and ignores only the outputs being rebuilt.
    products = registered()
    subprocess.run(["git", "update-index", "--skip-worktree", *products], cwd=ROOT, check=True)
    print(f"skip-worktree set on the {len(products)} registered products for the sweep")
    try:
        return sweep()
    finally:
        subprocess.run(["git", "update-index", "--no-skip-worktree", *products], cwd=ROOT,
                       check=True)
        print(f"skip-worktree cleared on the {len(products)} registered products")


def sweep() -> int:

    written = []
    build = LINDAL / "lindal_build.toml"
    written.append(build_raw.build(LINDAL / "raw", LINDAL / "raw" / "lindal_raw.nc"))
    written.append(build_gravity.build(build))
    written.append(build_rotation.build(build))
    written.append(build_wind.build(build))
    written.append(build_composition.build(build))
    written.extend(build_inputs.build(build))
    written.append(refrac_product.build_product(LINDAL / "lindal_reduction.toml"))
    written.extend(run_inputs.build(CLOSURE / "lindal_closure_build.toml"))
    written.append(forward_production.run(CLOSURE / "lindal_closure.toml").product)
    written.extend(run_inputs.build(TRANSFER / "lindal_transfer_build.toml"))
    written.append(forward_production.run(TRANSFER / "lindal_transfer.toml").product)

    print("\n--- casspian_git_commit read back from every rebuilt file ---")
    dirty = []
    for path in written:
        path = Path(path)
        if path.suffix != ".nc":
            continue
        with xr.open_dataset(path, engine="netcdf4") as handle:
            commit = str(handle.attrs.get("casspian_git_commit", ""))
        if commit.endswith("-dirty"):
            dirty.append(path)
        print(f"{path.relative_to(ROOT).as_posix():58s} {commit}")
    if dirty:
        print(f"REFUSED: {len(dirty)} file(s) carry -dirty: {[p.name for p in dirty]}")
        return 1

    print("\n--- the sweep hash table, step02_1 row format ---")
    for path in written:
        path = Path(path)
        print(f"| `{path.relative_to(ROOT).as_posix()}` | {kind_of(path)} | "
              f"`{cio.sha256(path)}` |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
