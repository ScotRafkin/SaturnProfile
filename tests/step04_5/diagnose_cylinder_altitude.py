"""Why the cylinder-extended run's delivered altitudes miss their stated values by 54 m.

SPEC_04 section 7's third run states altitudes 416,300, 99,290 and -15,456 m to 10 m (v0.10). Check
4 of the Step 5 acceptance measures 416,354, 99,309 and -15,450 at the pinned 0.05 degrees by
50,000 m2/s2, a departure of 53.79 m, while the same run's T, p and N pass at decision Q's floor and
the closure-wind runs 1, 2 and 4 reproduce their altitudes to 2.88 m.

Two hypotheses were killed before this script was written. The altitude does not depend on the node
set the target column is integrated on: subdividing the 67 arrival levels by up to 16 moves the
closure-wind altitudes by less than 0.1 m. And the cylinder field is not the difference: its wind on
the construction is 9.488, 3.716, 2.167 and 1.925 m/s at the four named levels, which is what the
Step 4 acceptance measured and reported for the same construction.

What this script measures is the geopotential spacing, the one thing decision R changed. Under the
closure wind `u` depends on latitude alone, so the isobar map never enters the target column: the
column asks the wind at `exp(ln_p_at(Phi))` and the answer does not depend on the pressure. Under
the cylinder wind `u` varies along a column, so the map's resolution in geopotential enters the
altitude through `|g_eff|`. The run is repeated at 50,000, 25,000 and 5,000 m2/s2 on one cylinder
field, and the arrival geopotential, the reference radius and the altitudes are reported at each, so
that what moves can be told from what does not.

The cylinder field is built once, at the namelist's latitude spacing, and written to
`reports/step04_5/step04_5_cylinder_wind.nc` so that a later run of this script can be read from it
rather than rebuilt; the construction takes about half an hour.
"""

import math
import sys
import time
from pathlib import Path

import numpy as np

import fields
from casspian.forward import estimate as fe
from casspian.forward import propagate as fprop
from casspian.forward import transfer as tr
from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib import windfield as wf

HERE = Path("reports/step04_5")
CYLINDER = HERE / "step04_5_cylinder_wind.nc"
STATED = {"top": 416300.0, "gauge": 99290.0, "bottom": -15456.0}

namelist = ctl.read_run_namelist("forward/lindal_transfer/lindal_transfer.toml")
inputs = ctl.load_run_inputs(namelist)
anchors = tuple(fprop.propagate(a, namelist.solar_longitude_deg) for a in inputs.anchors)
lindal = anchors[0]
closure_field = wf.WindField(inputs.wind)
p_b = tr.boundary_pressure(anchors)
target10 = math.radians(float(namelist.target_latitude_deg))
target60 = math.radians(60.0)
LAT = math.radians(float(namelist.grid["latitude_spacing_deg"]))
GEO = float(namelist.grid["geopotential_spacing_m2s2"])
LOOP = namelist.numerics["outer_loop"]
CONSTANTS = (float(inputs.rotation["angular_rate_rad_s"]), float(inputs.gravity["GM_m3s2"]),
             np.asarray(inputs.gravity["J"].values), np.asarray(inputs.gravity["degree"].values),
             float(inputs.gravity["normalization_radius_m"]))
p_tab = np.asarray(lindal.label_pressure_Pa, dtype="float64")

started = time.time()
if CYLINDER.exists() and "--rebuild" not in sys.argv:
    handle = cio.read(CYLINDER, "wind")
    try:
        cylinder_field = wf.WindField(handle.load())
    finally:
        handle.close()
    print(f"cylinder field read from {CYLINDER.name}", flush=True)
else:
    built = fields.cylinder_extended(inputs, lindal, closure_field, p_b, CONSTANTS,
                                     namelist.gauge_isobar_Pa, LAT,
                                     extra_latitudes=(target10, target60))
    cylinder_field = built["field"]
    written = inputs.wind.copy(deep=True)
    written["u_total_ms"] = (("latitude_planetocentric", "pressure"),
                             np.asarray(cylinder_field.u_total_ms),
                             dict(inputs.wind["u_total_ms"].attrs))
    written["u_shear_ms"] = (("latitude_planetocentric", "pressure"),
                             np.asarray(cylinder_field.u_total_ms)
                             - np.asarray(inputs.wind["u_reference_ms"].values)[:, None],
                             dict(inputs.wind["u_shear_ms"].attrs))
    written.attrs["title"] = "Cylinder-extended wind for the SPEC_04 Step 5 diagnosis"
    cio.write(CYLINDER, written, "wind", created_by="diagnose_cylinder_altitude")
    print(f"cylinder field built in {built['passes']} passes over {built['latitudes']} latitudes, "
          f"{time.time() - started:.0f} s, written to {CYLINDER.name}\n"
          f"on the construction u is "
          + ", ".join(f"{float(built['on_construction'][k]):.3f}"
                      for k in (0, lindal.gauge_level_index,
                                int(np.argmin(np.abs(p_tab - 1.0e5))), p_tab.size - 1))
          + " m/s (Step 4 measured 9.488, 3.716, 2.167, 1.925)", flush=True)

print(f"\nstated at 10 N under the cylinder wind: top {STATED['top']:.0f}, "
      f"gauge {STATED['gauge']:.0f}, bottom {STATED['bottom']:.0f} m, bound 10 m", flush=True)
print("the closure-wind run at the same latitude reproduces 411132, 98186, -15290 to 2.88 m\n",
      flush=True)


def run_at(geopotential_spacing, field, label):
    """One run, and what the altitude at the three named levels comes out as."""
    at = time.time()
    carried = tr.chain(
        inputs, anchors, gauge_latitude_rad=fe.gauge_latitude(anchors),
        target_latitude_rad=target10, gauge_isobar_Pa=namelist.gauge_isobar_Pa, p_b=p_b,
        latitude_spacing_rad=LAT, geopotential_spacing_m2s2=geopotential_spacing,
        relative_tolerance_ln_p=LOOP["relative_tolerance_ln_p"],
        max_iterations=LOOP["max_iterations"],
        kernel_uncertainty_per_rad=float(namelist.estimation["kernel_uncertainty_per_rad"]),
        datum_isobar_Pa=float(namelist.datum_isobar_Pa), field=field)
    produced = carried.target
    altitude = np.asarray(produced["altitude_m"], dtype="float64")
    gauge = int(produced["gauge_level_index"])
    Phi = np.asarray(produced["geopotential_m2s2"], dtype="float64")
    index = carried.state.mesh.latitude_index(target10)
    got = {"top": float(altitude[0]), "gauge": float(altitude[gauge]),
           "bottom": float(altitude[-1])}
    worst = max(abs(got[k] - STATED[k]) for k in STATED)
    print(f"{label}: {geopotential_spacing:.4g} m2/s2, mesh {carried.state.mesh.shape[0]} by "
          f"{carried.state.mesh.shape[1]}, {carried.state.record['passes']} passes, "
          f"{time.time() - at:.0f} s", flush=True)
    print(f"    altitude top {got['top']:.1f} ({got['top'] - STATED['top']:+.1f} from stated), "
          f"gauge {got['gauge']:.1f} ({got['gauge'] - STATED['gauge']:+.1f}), "
          f"bottom {got['bottom']:.1f} ({got['bottom'] - STATED['bottom']:+.1f}), "
          f"largest departure {worst:.2f} m", flush=True)
    print(f"    arrival Phi top {Phi[0]:.6e}, gauge {Phi[gauge]:.6e}, bottom {Phi[-1]:.6e} m2/s2; "
          f"r0 at the target {float(carried.state.reference_radius_m[index]):.3f} m; "
          f"datum Phi {float(produced['datum_geopotential_m2s2']):.6e}", flush=True)
    print(f"    largest |p/p_label - 1| "
          f"{float(np.abs(np.asarray(produced['pressure_identity_residual'])).max()):.3e}",
          flush=True)
    return got, Phi, float(carried.state.reference_radius_m[index])


cylinder = {}
for spacing, label in ((GEO, "cylinder wind at the namelist's spacing"),
                       (0.5 * GEO, "cylinder wind, halved"),
                       (0.1 * GEO, "cylinder wind at 5,000, the spacing Steps 0 to 4 ran at")):
    cylinder[spacing] = run_at(spacing, cylinder_field, label)

print("\nfor comparison, the closure wind, whose u does not vary along a column:", flush=True)
closure = {}
for spacing, label in ((GEO, "closure wind at the namelist's spacing"),
                       (0.1 * GEO, "closure wind at 5,000")):
    closure[spacing] = run_at(spacing, closure_field, label)

print(f"\nwhat moves with the spacing, cylinder wind, 50,000 to 5,000 m2/s2:", flush=True)
coarse, fine = cylinder[GEO], cylinder[0.1 * GEO]
for key in ("top", "gauge", "bottom"):
    print(f"    altitude at {key}: {coarse[0][key]:.1f} to {fine[0][key]:.1f} m, "
          f"moved {fine[0][key] - coarse[0][key]:+.1f}", flush=True)
print(f"    arrival Phi at the top level: {coarse[1][0]:.6e} to {fine[1][0]:.6e}, "
      f"moved {fine[1][0] - coarse[1][0]:+.3e} m2/s2", flush=True)
print(f"    r0 at the target: {coarse[2]:.3f} to {fine[2]:.3f} m, moved {fine[2] - coarse[2]:+.3f}",
      flush=True)
c_coarse, c_fine = closure[GEO], closure[0.1 * GEO]
print(f"the same for the closure wind:", flush=True)
for key in ("top", "gauge", "bottom"):
    print(f"    altitude at {key}: {c_coarse[0][key]:.1f} to {c_fine[0][key]:.1f} m, "
          f"moved {c_fine[0][key] - c_coarse[0][key]:+.1f}", flush=True)
print(f"\n{time.time() - started:.0f} s in all", flush=True)
