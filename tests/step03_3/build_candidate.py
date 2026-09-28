# The in-memory candidate kind N (SPEC_03 section 0 rule): used when a step's rebuilt inputs are -dirty before its acceptance commit.
"""The candidate kind N for the SPEC_03 Step 3 acceptance, built in memory (SPEC_03 v0.6 section 0).

Copies the manifest and the six rebuilt inputs to `reports/step03_3/candidate/` and runs
`refrac.product.build_product` there with the `-dirty` refusal relaxed by `relaxed.py`, then
renders the standard figures from the file written. Nothing under `occul_data/` is written.
"""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relaxed  # noqa: E402

from casspian.refrac.product import build_product  # noqa: E402
from casspian.tools.plots import describe, render  # noqa: E402

D = Path("occul_data/lindal")
CANDIDATE = Path("reports/step03_3/candidate")

if CANDIDATE.exists():
    shutil.rmtree(CANDIDATE)
CANDIDATE.mkdir(parents=True)
for name in ("lindal_reduction.toml", "lindal_thermo.nc", "lindal_composition.nc",
             "lindal_geodesy.nc", "lindal_gravity.nc", "lindal_rotation.nc", "lindal_wind.nc"):
    shutil.copy2(D / name, CANDIDATE / name)
product = build_product(CANDIDATE / "lindal_reduction.toml")
print(f"candidate product {product}")
print(f"inputs relaxed in memory: {[Path(p).name for p in relaxed.RELAXED]}")
for line in describe(render(product, CANDIDATE / "figures", "png", 150)):
    print(f"  {line}")
