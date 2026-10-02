"""`casspian-new-run`: a new named run, copied from an existing one. SPEC_05 Step 3 deliverable 4.

`casspian-new-run <name> --from <run directory>` makes `<name>/` beside the run it copies, holding
`<name>_build.toml` and `<name>.toml`: the source run's two control files with every whole-token
occurrence of the source run's name, comments included, replaced by the new name, and nothing
else changed. A whole token is the name not preceded by a letter, digit or underscore and not
followed by a letter or digit; an underscore may follow, since `<run>_` is how a run's files carry
its prefix. So the prefixes, output names, `[inputs]` paths, `[output] product`, `[run] name` and
titles follow, and a path that merely contains the name inside a longer word does not.

It makes no inputs. The person edits `[shear]`, `[target]` or whatever the experiment changes,
then runs `casspian-run-inputs` and `casspian-forward` as the runbook says. It refuses if the new
directory exists, and it refuses a source directory without its two control files.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from casspian.lib.control import ControlFileError

TOOL = "casspian-new-run"


def renamed(text: str, old: str, new: str) -> str:
    """`text` with every whole-token occurrence of `old` replaced by `new`."""
    return re.sub(rf"(?<![A-Za-z0-9_]){re.escape(old)}(?![A-Za-z0-9])", new, text)


def new_run(name: str, source_directory) -> Path:
    """Make the run `name` beside `source_directory`, from its two control files. Returns it."""
    source = Path(source_directory).resolve()
    old = source.name
    files = [source / f"{old}_build.toml", source / f"{old}.toml"]
    missing = [path.name for path in files if not path.exists()]
    if missing:
        raise ControlFileError(f"{source}: a run directory holds {old}_build.toml and {old}.toml; "
                               f"missing {missing}.")
    target = source.parent / name
    if target.exists():
        raise ControlFileError(f"{target} exists; casspian-new-run makes a new run and never "
                               "writes over one.")
    target.mkdir()
    for path in files:
        text = path.read_bytes().decode("utf-8")
        (target / path.name.replace(old, name, 1)).write_bytes(
            renamed(text, old, name).encode("utf-8"))
    return target


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Make a new forward run from an existing run's control files.")
    parser.add_argument("name", help="the new run's name, which becomes its prefix")
    parser.add_argument("--from", dest="source", required=True,
                        help="the run directory to copy, for example forward/lindal_transfer")
    args = parser.parse_args(argv)
    print(f"made {new_run(args.name, args.source)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
