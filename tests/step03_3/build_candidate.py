# The in-memory candidate kind N (SPEC_03 section 0 rule): used when a step's rebuilt inputs are -dirty before its acceptance commit.
"""The candidate kind N for the SPEC_03 Step 3 acceptance, built in memory (SPEC_03 v0.6 section 0).

Copies the manifest and the six rebuilt inputs to `reports/step03_3/candidate/` and runs
`refrac.product.build_product` there, then renders the standard figures from the file written.
Nothing under `occul_data/` is written. SPEC_R1 removed the `-dirty` refusal, and with it the
relaxation this script used to apply.
"""

import shutil
from pathlib import Path

from casspian.refrac.product import build_product
from casspian.tools.plots import describe, render

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
for line in describe(render(product, CANDIDATE / "figures", "png", 150)):
    print(f"  {line}")
