"""The run's target latitude from Lindal's egress label, converted as the anchor's label is.

SPEC_12 section 2 item 1: Lindal labels each profile with the latitude where the ray's lowest point
touched the 100 mbar surface, 36.3 N for ingress (the anchor) and 31.2 S for egress, planetographic.
The anchor's label became its planetocentric `phi_c` in the reduction, by `refrac.anchor`'s fixed
point (`lib.latitude.planetocentric_fixed_point`) on the wind-included geoid of the reduction's own
inputs, anchored by the manifest's rule, to the manifest's tolerance. This applies the same fixed
point to the egress label, and first reproduces the anchor's own `phi_c` from 36.3 N as the check.

Prints both and writes `reports/step11_2/target_latitude.txt`. Run from the repository root.
"""

from pathlib import Path

import numpy as np

from casspian.lib import control as ctl
from casspian.lib import latitude as lat
from casspian.refrac.anchor import geoid_setup

HERE = Path("reports/step11_2")
HERE.mkdir(parents=True, exist_ok=True)
EGRESS_LABEL_DEG = -31.2

manifest = ctl.read_reduction_manifest("occul_data/lindal/lindal_reduction.toml")
inputs = ctl.load_reduction_inputs(manifest)
setup = geoid_setup(inputs, manifest)


def convert(label_deg):
    phi_c, _, _, count = lat.planetocentric_fixed_point(
        np.radians(np.array([label_deg])), lambda phi: setup.march(phi).radius, setup.u_of_phi,
        *setup.constants, tol_rad=float(np.radians(manifest.fixed_point_tolerance_deg)),
        max_iter=manifest.max_iterations, flattening=setup.flattening)
    return float(np.degrees(phi_c[0])), int(count[0])


ingress_label = float(inputs.thermo.attrs["latitude_planetographic_deg"])
ingress, n_in = convert(ingress_label)
egress, n_out = convert(EGRESS_LABEL_DEG)
text = (f"ingress label {ingress_label:g} deg planetographic -> {ingress!r} deg planetocentric "
        f"({n_in} iterations); the anchor's frozen phi_c, for the check, is in the kind N product\n"
        f"egress label {EGRESS_LABEL_DEG:g} deg planetographic -> {egress!r} deg planetocentric "
        f"({n_out} iterations); tolerance {manifest.fixed_point_tolerance_deg:g} deg")
print(text)
(HERE / "target_latitude.txt").write_text(text + "\n", encoding="utf-8")
