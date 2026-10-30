"""`casspian-run-inputs`: a forward run's inputs, made by one command. SPEC_03 Step 3 deliverable 1.

A thin driver over the existing tools. It reads `<run>_build.toml` beside the namelist and runs its
five sections in order, gravity, rotation, wind, composition and shear, each with the tool that
writes that kind (SPEC_05 Step 3: `[wind]` writes the source wind, `[shear]` the wind the run
reads, `case = "identity"` when the run wants the source as it is; SPEC_11 v0.4: `[composition]`
comes before `[shear]`, whose `lindal_iris` case reads the run's composition). It contains no
physics and no value; every choice is in the build file, and each tool refuses what it refuses on
its own.

What the driver itself checks, because only it knows the file is a run's build file: it carries
exactly the five sections (no `[stage_two]`, which has no forward counterpart); and every output
is a file inside the run's directory, since a run writes only inside its own directory (SPEC_R1
section 2). Names, prefixes and roles are not checked (SPEC_R1).
"""

from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

from casspian.lib.control import ControlFileError
from casspian.tools.composition import build_composition
from casspian.tools.gravity import build_gravity, build_rotation
from casspian.tools.wind import build_wind, shear

TOOL = "casspian-run-inputs"

#: The sections in the order they run: wind reads the gravity and rotation files the first two
#: write, and shear reads the source wind that wind writes and the composition before it.
SECTIONS = (
    ("gravity", build_gravity.build),
    ("rotation", build_rotation.build),
    ("wind", build_wind.build),
    ("composition", build_composition.build),
    ("shear", shear.build),
)


def build(control_path) -> list[Path]:
    """Run the five sections of a run's build file. Returns the paths written, in order."""
    path = Path(control_path).resolve()
    if not path.exists():
        raise ControlFileError(f"{path}: build control file does not exist")
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
    run_directory = path.parent
    for name in names:
        table = document[name]
        output = (run_directory / str(table.get("output", ""))).resolve()
        if not output.is_relative_to(run_directory) or output == run_directory:
            raise ControlFileError(
                f"{path}: [{name}] output = {table.get('output')!r} is not a file inside the run's "
                "directory; a run writes only inside its own directory (SPEC_R1 section 2)."
            )

    (run_directory / "inputs").mkdir(parents=True, exist_ok=True)
    return [builder(path, name) for name, builder in SECTIONS]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Build a forward run's inputs from its build control file.")
    parser.add_argument("control", help="path to <run>_build.toml")
    args = parser.parse_args(argv)
    for written in build(args.control):
        print(f"wrote {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
