"""Is the M = 2 residual decision G, or my interpolation of the column onto the traced Phi?

The anchor's `radius_m` and height were read off the state's column by interpolating from the
mesh's geopotential nodes onto the traced `Phi`. That interpolation is mine and it is O(dPhi^2).
This integrates the 60 N column directly on the traced `Phi` instead, which removes it, and
measures what is left.
"""
import math
import shutil
from pathlib import Path

import netCDF4
import numpy as np

from casspian.forward import estimate as fe
from casspian.forward import production as fp
from casspian.forward import propagate as fprop
from casspian.forward import transfer as tr
from casspian.lib import control as ctl
from casspian.lib import geoid
from casspian.lib import hydrostatic as hs
from casspian.lib import kernel as lk
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf

HERE = Path("reports/step04_4")
LAT = math.radians(float(__import__("sys").argv[1]) if len(__import__("sys").argv) > 1 else 0.5)
PHI = float(__import__("sys").argv[2]) if len(__import__("sys").argv) > 2 else 5.0e4

namelist = ctl.read_run_namelist("forward/lindal_transfer/lindal_transfer.toml")
inputs = ctl.load_run_inputs(namelist)
lindal = fprop.propagate(inputs.anchors[0], namelist.solar_longitude_deg)
field = wf.WindField(inputs.wind)
CONSTANTS = (float(inputs.rotation["angular_rate_rad_s"]), float(inputs.gravity["GM_m3s2"]),
             np.asarray(inputs.gravity["J"].values), np.asarray(inputs.gravity["degree"].values),
             float(inputs.gravity["normalization_radius_m"]))
phi_a = math.radians(lindal.latitude_planetocentric_deg)
target10, target60 = math.radians(10.0), math.radians(60.0)
p_tab = np.asarray(lindal.label_pressure_Pa, dtype="float64")
p_b = tr.boundary_pressure([lindal])
phi_r = 0.5 * (phi_a + target60)

Phi_lindal = tr.place(lindal, inputs, field, p_b).geopotential_m2s2
mesh = lm.build_mesh([target10, phi_a, phi_r, target60], Phi_lindal, LAT, PHI)
state = tr.outer_loop(inputs, [lindal], gauge_latitude_rad=phi_r,
                      target_latitude_rad=target60, gauge_isobar_Pa=namelist.gauge_isobar_Pa,
                      p_b=p_b, latitude_spacing_rad=LAT, geopotential_spacing_m2s2=PHI,
                      relative_tolerance_ln_p=1e-8, max_iterations=50, field=field, mesh=mesh)
curves = state.curves[0]
latitude_c, slope = lk.composition_term(inputs.composition,
                                        state.placements[0].label_pressure_Pa)
ln_N = tr.transfer(mesh, state.s_over_g, curves, lindal.ln_N, latitude_c, slope)
index60 = mesh.latitude_index(target60)
Phi_traced = curves.geopotential_m2s2[:, index60]
N60 = np.exp(ln_N[:, index60])
column = state.columns[index60]
R_bar60, _ = fp.mean_properties(inputs.composition, 60.0)
relative = (np.asarray(lindal.tree.to_dataset(inherit=False)["refractivity_uncertainty"].values,
                       dtype="float64") / np.exp(lindal.ln_N))

# A: the column interpolated from the mesh's nodes onto the traced Phi, which is what the
# acceptance script does now.
rising = np.argsort(Phi_traced)
back = np.argsort(rising)
radius_A = np.interp(Phi_traced[rising], column.geopotential_m2s2, column.radius_m)[back]
z_A = np.interp(Phi_traced[rising], column.geopotential_m2s2, column.z_local_vertical_m)[back]

# B: the column integrated on the traced Phi itself, which removes that interpolation. The nodes
# must contain Phi = 0, where the column starts, and be strictly increasing.
nodes = np.unique(np.concatenate([Phi_traced, [0.0]]))
direct = lm.column(
    target60, nodes, float(state.reference_radius_m[index60]),
    (lambda Phi: field.wind_at(target60, np.exp(state.isobar_map.ln_p_at(index60, Phi)))),
    *CONSTANTS)
at_traced = np.searchsorted(nodes, Phi_traced)
radius_B, z_B = direct.radius_m[at_traced], direct.z_local_vertical_m[at_traced]

print("spacings %.3g degrees by %.3g m2/s2, mesh %s"
      % (math.degrees(LAT), PHI, (mesh.shape,)))
print("radius A against B: largest %.4e m" % float(np.max(np.abs(radius_A - radius_B))))
print("z_lv   A against B: largest %.4e m" % float(np.max(np.abs(z_A - z_B))))

_original = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda p, a, c, what=None: str(a.get("casspian_git_commit", ""))
own_at_gauge = fe.Arrival(
    slug="lindal", latitude_rad=phi_a, weight=1.0,
    geopotential_m2s2=curves.geopotential_m2s2[:, mesh.latitude_index(phi_r)],
    C_i=ln_N[:, mesh.latitude_index(phi_r)], sigma_ln_N=fe.sigma_ln_N(lindal),
    label_pressure_Pa=state.placements[0].label_pressure_Pa)

for name, radius, z_lv in (("A, interpolated from the mesh's nodes", radius_A, z_A),
                           ("B, integrated on the traced Phi", radius_B, z_B)):
    path = HERE / f"probe{'A' if name.startswith('A') else 'B'}_refractivity.nc"
    shutil.copy(lindal.path, path)
    with netCDF4.Dataset(path, "a") as handle:
        handle["radius_m"][:] = radius
        handle["height_above_anchor_isobar_m"][:] = z_lv
        handle["refractivity"][:] = N60
        handle["refractivity_uncertainty"][:] = relative * N60
        handle["latitude_planetocentric_deg"][...] = 60.0
        handle["psi_deg"][...] = float(np.degrees(column.psi_rad[mesh.gauge_geopotential_index]))
        handle["anchor_isobar_radius_m"][...] = float(
            column.radius_m[mesh.gauge_geopotential_index])
        handle["inputs/thermo"]["temperature_K"][:] = hs.temperature(p_tab, N60, R_bar60)
        handle.setncattr("profile_or_run", path.name.split("_refractivity")[0])
    handle = ctl.cio.read(path, "refractivity")
    try:
        tree = handle.load()
    finally:
        handle.close()
    loaded = ctl._load_anchor_object(
        ctl.RunAnchor(slug=path.name.split("_refractivity")[0], path=path, weight=1.0,
                      measurement_uncertainty_scale=1.0), namelist, tree)
    placed = tr.place(loaded, inputs, field, p_b)
    gap = float(np.max(np.abs(placed.geopotential_m2s2 - Phi_traced)))
    back_curves = tr.trace(mesh, state.shear_integral, placed.geopotential_m2s2, target60, phi_r)
    at_gauge = int(np.flatnonzero(back_curves.latitude_rad == phi_r)[0])
    lc, sl = lk.composition_term(inputs.composition, placed.label_pressure_Pa)
    back_lnN = tr.transfer(mesh, state.s_over_g, back_curves, np.log(N60), lc, sl)
    E = fe.estimate([own_at_gauge, fe.Arrival(
        slug="probe", latitude_rad=target60, weight=1.0,
        geopotential_m2s2=back_curves.geopotential_m2s2[:, at_gauge],
        C_i=back_lnN[:, at_gauge], sigma_ln_N=fe.sigma_ln_N(lindal),
        label_pressure_Pa=placed.label_pressure_Pa)], 0.02, phi_r)
    D = np.asarray(list(E.differences.values())[0], dtype="float64")
    print("%s: place()'s Phi against the traced %.4e m2/s2, |D_12| %.4e"
          % (name, gap, float(np.max(np.abs(D[np.isfinite(D)])))))
ctl._refuse_dirty_commit = _original
