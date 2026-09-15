"""`casspian-run-inputs`: a forward run's inputs, made by one command. SPEC_03 Step 3 deliverable 1.

A thin driver over the existing tools. It reads `<run>_build.toml` beside the namelist and runs its
four sections in order, gravity, rotation, wind and composition, each with the tool that writes
that kind. It contains no physics and no value; every choice is in the build file, and each tool
refuses what it refuses on its own.

What the driver itself checks, because only it knows the file is a run's build file: the file is
named `<run>_build.toml`; it carries exactly the four sections (no `[stage_two]`, which has no
forward counterpart); every section declares `role = "forward"` and the run's prefix; and every
output is a file under the run's `inputs/` directory (SPEC_00 section 2.3: the tools write
`inputs/`, `forward` writes `output/`).
"""

from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

from casspian.lib.control import ControlFileError
from casspian.tools.composition import build_composition
from casspian.tools.gravity import build_gravity, build_rotation
from casspian.tools.wind import build_wind

TOOL = "casspian-run-inputs"

#: The sections in the order they run: wind reads the gravity and rotation files the first two
#: write.
SECTIONS = (
    ("gravity", build_gravity.build),
    ("rotation", build_rotation.build),
    ("wind", build_wind.build),
    ("composition", build_composition.build),
)


def build(control_path) -> list[Path]:
    """Run the four sections of a run's build file. Returns the paths written, in order."""
    path = Path(control_path).resolve()
    if not path.exists():
        raise ControlFileError(f"{path}: build control file does not exist")
    if not path.name.endswith("_build.toml"):
        raise ControlFileError(
            f"{path}: a run's build control file is <run>_build.toml beside the namelist (SPEC_00 "
            "section 2.3)."
        )
    run = path.name[: -len("_build.toml")]
    with open(path, "rb") as handle:
        document = tomllib.load(handle)

    names = [name for name, _ in SECTIONS]
    if "stage_two" in document:
        raise ControlFileError(
            f"{path}: [stage_two] has no forward counterpart; a run's build file carries "
            f"{names} (SPEC_03 Step 3 deliverable 1)."
        )
    unknown = sorted(set(document) - set(names))
    missing = [name for name in names if name not in document]
    if unknown or missing:
        raise ControlFileError(
            f"{path}: a run's build file carries exactly the sections {names}; unknown "
            f"{unknown or 'none'}, missing {missing or 'none'}."
        )
    inputs_directory = path.parent / "inputs"
    for name in names:
        table = document[name]
        if table.get("role") != "forward":
            raise ControlFileError(
                f"{path}: [{name}] role = {table.get('role')!r}; a run's inputs are role 'forward' "
                "(SPEC_03 Step 3 deliverable 1)."
            )
        if table.get("prefix") != run:
            raise ControlFileError(
                f"{path}: [{name}] prefix = {table.get('prefix')!r}; a run's inputs carry the run "
                f"prefix {run!r} (SPEC_00 section 8)."
            )
        output = (path.parent / str(table.get("output", ""))).resolve()
        if output.parent != inputs_directory or not output.name.startswith(f"{run}_"):
            raise ControlFileError(
                f"{path}: [{name}] output = {table.get('output')!r} is not a file under the run's "
                f"inputs/ directory carrying the prefix '{run}_' (SPEC_00 section 2.3)."
            )

    inputs_directory.mkdir(parents=True, exist_ok=True)
    return [builder(path, name) for name, builder in SECTIONS]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Build a forward run's four inputs from its build control file.")
    parser.add_argument("control", help="path to <run>_build.toml")
    args = parser.parse_args(argv)
    for written in build(args.control):
        print(f"wrote {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
