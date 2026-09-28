"""The `-dirty` input refusal relaxed in this process only, for the SPEC_03 Step 0 acceptance.

Author ruling of 14 September 2026 (option 1): the Step 0 inputs cannot be clean before the
acceptance commit, so the acceptance builds its candidate kind N in memory with the refusal of
`lib.control.load_reduction_inputs` relaxed here, never in the code. Each input's commit is
kept, with the suffix renamed so the check passes and the relaxation stays visible.

Usage as a runner: `python tests/step03_0/relaxed.py <suite.py>` runs a suite with the
relaxation applied.
"""

import runpy
import sys

from casspian.lib import control as ctl

MARK = "-dirty (relaxed in the Step 03_0 acceptance)"
RELAXED = []
_original = ctl._load_into_memory


def _relaxed(path, kind):
    dataset = _original(path, kind)
    commit = str(dataset.attrs.get("casspian_git_commit", ""))
    if commit.endswith("-dirty"):
        dataset.attrs["casspian_git_commit"] = commit[: -len("-dirty")] + MARK
        RELAXED.append(str(path))
    return dataset


ctl._load_into_memory = _relaxed

if __name__ == "__main__":
    suite = sys.argv[1]
    sys.argv = [suite]
    print(f"[runner] -dirty input refusal relaxed in memory for {suite}")
    runpy.run_path(suite, run_name="__main__")
