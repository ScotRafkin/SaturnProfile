# The candidate rebuild: the chain, closure and transfer inputs rebuilt -dirty in the working tree, for a step that changes an input file.
"""Rebuild the reduction chain, the closure run and the transfer run's inputs.

SPEC_04 Step 0 deliverable 1, the cascade of SPEC_04 section 10 finding 2. Kind W changes, kind
N records the wind file's hash in `input_hashes` and `reduction_record`, and the reader refuses
a mismatch, so the whole chain follows: the six reduction inputs, the kind N product, the
closure run's four inputs and the closure product, and then the transfer run's four inputs.

**The relaxation, named here and in the output** (SPEC_03 section 0, the author's ruling of 14
September 2026, carried into SPEC_04 section 0). This step changes an input file, so it is
accepted on candidates built on the working tree, which carries the step's own uncommitted
source. Every product this script writes therefore carries `casspian_git_commit` ending in
`-dirty`, and `refrac` and `forward` refuse a `-dirty` input. The single refusal point,
`casspian.lib.control._refuse_dirty_commit`, is relaxed for the duration of this script and
nowhere else; nothing else is patched, no commit string is altered on disk, and every product
keeps the honest `-dirty` mark it earns. The registered products are rebuilt from the clean
tree at the sweep, and their hashes recorded in the report then.

Run from the repository root.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

from casspian.lib import control as ctl  # noqa: E402
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


def relax_dirty_refusal():
    """Relax the one `-dirty` refusal point, and say so. Returns the original function."""
    original = ctl._refuse_dirty_commit

    def permissive(path, attrs, consumer, what=None):
        return str(attrs.get("casspian_git_commit", ""))

    ctl._refuse_dirty_commit = permissive
    print("RELAXATION: casspian.lib.control._refuse_dirty_commit is replaced for this script "
          "only, so that refrac and forward accept the -dirty products this working tree "
          "builds (SPEC_03 section 0, in-memory candidate rule). Nothing else is patched.")
    return original


def main() -> int:
    original = relax_dirty_refusal()
    try:
        written = []
        print("\n--- the reduction chain ---")
        build = LINDAL / "lindal_build.toml"
        written.append(build_raw.build(LINDAL / "raw", LINDAL / "raw" / "lindal_raw.nc"))
        written.append(build_gravity.build(build))
        written.append(build_rotation.build(build))
        written.append(build_wind.build(build))
        written.append(build_composition.build(build))
        written.extend(build_inputs.build(build))
        written.append(refrac_product.build_product(LINDAL / "lindal_reduction.toml"))

        print("\n--- the closure run ---")
        written.extend(run_inputs.build(CLOSURE / "lindal_closure_build.toml"))
        written.append(forward_production.run(CLOSURE / "lindal_closure.toml").product)

        print("\n--- the transfer run's inputs ---")
        written.extend(run_inputs.build(TRANSFER / "lindal_transfer_build.toml"))
    finally:
        ctl._refuse_dirty_commit = original

    print("\n--- what was written ---")
    for path in written:
        path = Path(path)
        print(f"| `{path.relative_to(ROOT).as_posix()}` | {cio.sha256(path)} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
