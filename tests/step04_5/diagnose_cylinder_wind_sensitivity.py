"""How large a difference in `u` at the target column the cylinder run's 42 m of altitude is.

`diagnose_cylinder_altitude.py` establishes that the departure is not the mesh: from 50,000 to
5,000 m2/s2 the top altitude moves 12.1 m and the remaining 41.7 m does not converge, while under
the closure wind the same path is insensitive to the spacing to 0.0 m and reproduces its own stated
altitudes to 2.88 m. What is left is `|g_eff|` along the target column, and the only thing in it that
the cylinder wind changes is `u`.

This measures `d(altitude) / d(u)` at the target directly: the same column, integrated with `u`
scaled, using the cylinder field that `diagnose_cylinder_altitude.py` wrote. It says what change in
`u` a 42 m change in the delivered altitude corresponds to, so that the size of the departure can be
compared with the size of the differences between cylinder constructions.
"""

import math

import numpy as np

from casspian.forward import estimate as fe
from casspian.forward import propagate as fprop
from casspian.forward import transfer as tr
from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf

CYLINDER = "reports/step04_5/step04_5_cylinder_wind.nc"
STATED_TOP = 416300.0

namelist = ctl.read_run_namelist("forward/lindal_transfer/lindal_transfer.toml")
inputs = ctl.load_run_inputs(namelist)
anchors = tuple(fprop.propagate(a, namelist.solar_longitude_deg) for a in inputs.anchors)
closure_field = wf.WindField(inputs.wind)
handle = cio.read(CYLINDER, "wind")
try:
    cylinder_field = wf.WindField(handle.load())
finally:
    handle.close()

target10 = math.radians(float(namelist.target_latitude_deg))
LAT = math.radians(float(namelist.grid["latitude_spacing_deg"]))
GEO = float(namelist.grid["geopotential_spacing_m2s2"])
LOOP = namelist.numerics["outer_loop"]
CONSTANTS = (float(inputs.rotation["angular_rate_rad_s"]), float(inputs.gravity["GM_m3s2"]),
             np.asarray(inputs.gravity["J"].values), np.asarray(inputs.gravity["degree"].values),
             float(inputs.gravity["normalization_radius_m"]))

carried = tr.chain(
    inputs, anchors, gauge_latitude_rad=fe.gauge_latitude(anchors),
    target_latitude_rad=target10, gauge_isobar_Pa=namelist.gauge_isobar_Pa,
    p_b=tr.boundary_pressure(anchors), latitude_spacing_rad=LAT, geopotential_spacing_m2s2=GEO,
    relative_tolerance_ln_p=LOOP["relative_tolerance_ln_p"],
    max_iterations=LOOP["max_iterations"],
    kernel_uncertainty_per_rad=float(namelist.estimation["kernel_uncertainty_per_rad"]),
    datum_isobar_Pa=float(namelist.datum_isobar_Pa), field=cylinder_field)
produced = carried.target
state = carried.state
index = state.mesh.latitude_index(target10)
Phi = np.asarray(produced["geopotential_m2s2"], dtype="float64")
datum_Phi = float(produced["datum_geopotential_m2s2"])
labels = np.asarray(produced["label_pressure_Pa"], dtype="float64")

u_cylinder = np.asarray(cylinder_field.wind_at(np.full(labels.shape, float(target10)), labels),
                        dtype="float64")
u_closure = np.asarray(closure_field.wind_at(np.full(labels.shape, float(target10)), labels),
                       dtype="float64")
print(f"at the target, {math.degrees(target10):g} degrees, on the delivered labels:")
print(f"    cylinder u from {u_cylinder.min():.3f} to {u_cylinder.max():.3f} m/s, "
      f"top level {u_cylinder[0]:.3f}, bottom {u_cylinder[-1]:.3f}")
print(f"    closure  u from {u_closure.min():.3f} to {u_closure.max():.3f} m/s, "
      f"top level {u_closure[0]:.3f}, bottom {u_closure[-1]:.3f}")
print(f"    the two differ by up to {np.abs(u_cylinder - u_closure).max():.3f} m/s, and the "
      f"delivered top altitude differs by 416354 - 411135 = 5219 m")


def top_altitude(scale=1.0, offset=0.0):
    """The delivered top altitude with `u` at the target scaled and offset."""
    nodes = np.unique(np.concatenate([Phi, [0.0, datum_Phi]]))
    column = lm.column(
        float(target10), nodes, float(state.reference_radius_m[index]),
        (lambda Phi_at: scale * np.asarray(cylinder_field.wind_at(
            float(target10), np.exp(state.isobar_map.ln_p_at(index, Phi_at))),
            dtype="float64") + offset),
        *CONSTANTS)
    at = np.searchsorted(nodes, Phi)
    z = column.z_local_vertical_m[at]
    z_datum = float(np.interp(datum_Phi, column.geopotential_m2s2, column.z_local_vertical_m))
    return float((z - z_datum)[0])


base = top_altitude()
print(f"\nthe top altitude on this state is {base:.1f} m, stated {STATED_TOP:.0f}, "
      f"departure {base - STATED_TOP:+.1f} m")
print("with u at the target offset by a constant:")
for offset in (-1.0, -0.5, -0.1, 0.1, 0.5, 1.0):
    moved = top_altitude(offset=offset)
    print(f"    u + {offset:+.2f} m/s: top {moved:.1f} m, moved {moved - base:+.2f}, "
          f"so {abs((moved - base) / offset):.1f} m per m/s")
per = (top_altitude(offset=0.1) - top_altitude(offset=-0.1)) / 0.2
print(f"\nd(top altitude)/du = {per:.1f} m per m/s at the target column, so the {base - STATED_TOP:.1f} m "
      f"departure is {(base - STATED_TOP) / per:+.3f} m/s in u, against a cylinder u there of "
      f"{u_cylinder[0]:.3f} m/s at the top level and {np.abs(u_cylinder).max():.3f} m/s at most")
