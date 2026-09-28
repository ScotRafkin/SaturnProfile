"""The `-dirty` input refusal relaxed in this process only, for the SPEC_03 Step 3 acceptance.

SPEC_03 v0.6 section 0, for a step that changes an input file: the rebuilt inputs cannot be
clean before the acceptance commit, so the acceptance builds its candidates in memory with the
refusal of `lib.control` relaxed here, never in the code. Each input's commit is kept, with the
suffix renamed so the check passes and the relaxation stays visible. `RELAXED` lists the files.
`disabled()` restores the original loader for a case that must see the refusal itself.

Usage as a runner: `python tests/step03_3/relaxed.py <suite.py>` runs a suite with the
relaxation applied.
"""

import runpy
import sys
from contextlib import contextmanager

from casspian.lib import control as ctl

MARK = "-dirty (relaxed in the Step 03_3 acceptance)"
RELAXED = []
_original = ctl._load_into_memory


def _relaxed(path, kind):
    dataset = _original(path, kind)
    commit = str(dataset.attrs.get("casspian_git_commit", ""))
    if commit.endswith("-dirty"):
        dataset.attrs["casspian_git_commit"] = commit[: -len("-dirty")] + MARK
        RELAXED.append(str(path))
    return dataset


@contextmanager
def disabled():
    """The original loader, for the duration of a case that tests the refusal."""
    ctl._load_into_memory = _original
    try:
        yield
    finally:
        ctl._load_into_memory = _relaxed


ctl._load_into_memory = _relaxed

if __name__ == "__main__":
    suite = sys.argv[1]
    sys.argv = [suite]
    print(f"[runner] -dirty input refusal relaxed in memory for {suite}")
    runpy.run_path(suite, run_name="__main__")
