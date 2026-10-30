"""The three new targets' planetocentric latitudes, converted as SPEC_12 converted 31.2 S.

SPEC_15 section 1 deliverable 5. Lindal labels each profile in planetographic latitude (SPEC_15
section 2 ruling 1). Each label becomes planetocentric by the reduction's own fixed point
(`lib.latitude.planetocentric_fixed_point`) on the wind-included geoid of the reduction's inputs, as
`tests/step11_2/target_latitude.py` converted the egress label; that script's two conversions are
repeated first, as the check.

Prints the latitudes and writes `reports/step15_1/target_latitudes.txt`. Run from the repository root.
"""

from pathlib import Path

import numpy as np

from casspian.lib import control as ctl
from casspian.lib import latitude as lat
from casspian.refrac.anchor import geoid_setup

HERE = Path("reports/step15_1")
HERE.mkdir(parents=True, exist_ok=True)
LABELS_DEG = {"36.3 N (the anchor, the check)": 36.3, "31.2 S (SPEC_12, the check)": -31.2,
              "74.2 S": -74.2, "3.0 S": -3.0, "11.9 S": -11.9}

manifest = ctl.read_reduction_manifest("occul_data/lindal/lindal_reduction.toml")
inputs = ctl.load_reduction_inputs(manifest)
setup = geoid_setup(inputs, manifest)


def convert(label_deg):
    phi_c, _, _, count = lat.planetocentric_fixed_point(
        np.radians(np.array([label_deg])), lambda phi: setup.march(phi).radius, setup.u_of_phi,
        *setup.constants, tol_rad=float(np.radians(manifest.fixed_point_tolerance_deg)),
        max_iter=manifest.max_iterations, flattening=setup.flattening)
    return float(np.degrees(phi_c[0])), int(count[0])


lines = []
for name, label in LABELS_DEG.items():
    value, count = convert(label)
    lines.append(f"{name}: {label:g} deg planetographic -> {value!r} deg planetocentric ({count} iterations)")
text = "\n".join(lines) + f"\ntolerance {manifest.fixed_point_tolerance_deg:g} deg"
print(text)
(HERE / "target_latitudes.txt").write_text(text + "\n", encoding="utf-8", newline="\n")
