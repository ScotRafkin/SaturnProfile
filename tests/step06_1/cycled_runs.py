"""SPEC_06 Step 1, reported not checked: the two runs that cycled in SPEC_05 Step 4.

Run 6 at 5e4 and run 7 at 2.5e4 stopped at the outer loop's tolerance of 1e-8 with the residual
held at a constant 1.5e-4 and 1.4e-4 (REPORT_05_step4 section 7). Their namelists are now at 1e-3,
so each is rerun here in a copy under `reports/step06_1/cycled/` at 1e-8, with the new kernel and
the per-pass residual printed, to see whether the mesh blend at the kink set them cycling. The
`-dirty` refusal is relaxed in this process only. Run from the repository root.
"""

import shutil
import time
import traceback
import warnings
from pathlib import Path

from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.tools.run import run_inputs

HERE = Path("reports/step06_1/cycled")
CASES = (("shear_r6_decay40", "5.0e4"), ("shear_r7_decay_linp", "2.5e4"))


def copy_at_1e8(run, spacing):
    source = Path("forward") / run
    copy = HERE / f"{run}_{spacing}"
    if copy.exists():
        shutil.rmtree(copy)
    copy.mkdir(parents=True)
    for suffix in ("_build.toml", ".toml"):
        text = (source / f"{run}{suffix}").read_text(encoding="utf-8").replace('"../../', '"../../../../')
        if suffix == ".toml":
            for old, new in (("relative_tolerance_ln_p = 1.0e-3", "relative_tolerance_ln_p = 1.0e-8"),
                             ("geopotential_spacing_m2s2 = 5.0e4",
                              f"geopotential_spacing_m2s2 = {spacing}"),
                             ("figures = true", "figures = false")):
                if text.count(old) != 1:
                    raise RuntimeError(f"{run}{suffix}: {old!r} matched {text.count(old)} times")
                text = text.replace(old, new)
        (copy / f"{run}{suffix}").write_bytes(text.encode("utf-8"))
    run_inputs.build(copy / f"{run}_build.toml")
    return copy / f"{run}.toml"


original = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
    attrs.get("casspian_git_commit", ""))
print("RELAXATION: casspian.lib.control._refuse_dirty_commit replaced for this process only\n")
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for run, spacing in CASES:
            started = time.time()
            try:
                product = forward_transfer.run(copy_at_1e8(run, spacing)).product
                print(f"{run} at {spacing}, tolerance 1e-8: completed in {time.time() - started:.0f} s, "
                      f"{product}\n", flush=True)
            except Exception as exc:  # a named failure is a result
                print(f"{run} at {spacing}, tolerance 1e-8: {type(exc).__name__}: {exc}\n"
                      + traceback.format_exc().strip().splitlines()[-1] + "\n", flush=True)
finally:
    ctl._refuse_dirty_commit = original
