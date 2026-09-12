"""The two constants frozen per profile: the anchor latitude `phi_c` and the anchor radius `r0`.

SPEC_02 Step 2; handoff section 9A; SPEC_01 Steps 5 to 7. Built from `lib` only, on the
source's own G, R, W, D and T, as `lib.control.load_reduction_inputs` delivered them.

This module is the `refrac` boundary for angles. The manifest states its fixed point tolerance
in degrees and kind T states the label latitude in degrees; both are converted to radians here,
once, and everything passed into `lib` and kept in the record is in radians.

**The anchor radius is marched, never interpolated** (SPEC_01 v0.7). Every surface evaluation
inside the fixed point is a full anchored march of Eq. B3 with the iterate's latitude inserted
as a node, and `r0` is the radius at the node `phi_c` of one final march. That final march is
also the surface whose polar radii, asymmetry and anchoring residual the record carries, so
every geoid number in the record describes the one surface `r0` was read from.

The no wind values are computed alongside for the record and are not used: the same fixed
point on `lib.geoid.reference_geoid`, solved at each iterate, which are the handoff section
9A.8 numbers. The difference between the two pairs is the effect of the wind on the anchor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from casspian.lib import geoid as gd
from casspian.lib import latitude as lat
from casspian.lib.control import ReductionInputs, ReductionManifest
from casspian.lib.gravity import g_eff_vector

__all__ = ["AnchorConvergenceError", "FrozenAnchor", "freeze_anchor", "wind_of_latitude"]


class AnchorConvergenceError(RuntimeError):
    """The anchor latitude fixed point did not reach the manifest tolerance."""


@dataclass(frozen=True)
class FrozenAnchor:
    """The frozen pair and everything needed to say how it was reached. Angles in radians."""

    phi_g_rad: float
    phi_c_rad: float
    psi_rad: float
    r0_m: float
    u_at_phi_c_ms: float
    iterates_rad: np.ndarray
    iteration_count: int
    final_step_rad: float
    tolerance_rad: float
    anchor_rule: str
    anchor_surface_Pa: float
    anchor_quantity: str
    r_anchor_m: float
    north_start_m: float
    polar_north_m: float
    polar_south_m: float
    polar_asymmetry_m: float
    anchor_residual_m: float
    nowind_phi_c_rad: float
    nowind_psi_rad: float
    nowind_r0_m: float
    nowind_iterates_rad: np.ndarray


def wind_of_latitude(wind):
    """`u(phi_c)` on the reference level of a kind W dataset, as a callable in radians.

    Linear in planetocentric latitude on the file's own grid. Kind W carries both poles as nodes
    with exactly zero wind (SPEC_00 section 6.6, checked by `read`), so the callable is zero at
    the poles by the file's own values and needs no rule of its own. Returns an array shaped like
    its argument, as `lib.geoid.wind_geoid` requires.
    """
    latitude = np.radians(np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64"))
    pressure = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    reference = float(wind["reference_level_pressure_Pa"])
    column = int(np.flatnonzero(pressure == reference)[0])
    u = np.asarray(wind["u_total_ms"].values, dtype="float64")[:, column]
    order = np.argsort(latitude)
    latitude, u = latitude[order], u[order]

    def u_of_phi(phi):
        phi = np.asarray(phi, dtype="float64")
        return np.reshape(np.interp(phi, latitude, u), phi.shape)

    return u_of_phi


def _final_step(iterates):
    """The size of the last fixed point update, in radians."""
    if len(iterates) < 2:
        return float("inf")
    return float(abs(iterates[-1][0] - iterates[-2][0]))


def freeze_anchor(inputs: ReductionInputs, manifest: ReductionManifest) -> FrozenAnchor:
    """Solve for `phi_c` and `r0` on the wind included geoid. SPEC_02 Step 2.

    1. `u_of_phi` from the kind W reference level.
    2. The wind geoid by `lib.geoid.wind_geoid`, anchored by the manifest's `anchor_rule` on the
       geodesy file's `anchor_quantity` at `anchor_surface_Pa`, converged to `convergence_m`.
    3. `phi_c` by `lib.latitude.planetocentric_fixed_point` from the kind T label, the surface
       marched at each iterate, seeded on the anchor surface's oblateness, to the manifest
       tolerance in radians within `max_iterations`.
    4. `r0` from one final march with `phi_c` a node, and `psi` evaluated there.
    5. The no wind pair by the same fixed point on the reference geoid, for the record.

    Refuses, with `AnchorConvergenceError`, a fixed point that exhausts `max_iterations`.
    """
    gravity, rotation, geodesy, thermo = (inputs.gravity, inputs.rotation, inputs.geodesy,
                                          inputs.thermo)
    degrees = np.asarray(gravity["degree"].values)
    J = np.asarray(gravity["J"].values, dtype="float64")
    GM = float(gravity["GM_m3s2"])
    R_norm = float(gravity["normalization_radius_m"])
    Omega = float(rotation["angular_rate_rad_s"])
    constants = (Omega, GM, J, degrees, R_norm)

    surfaces = np.asarray(geodesy["surface_pressure_Pa"].values, dtype="float64")
    row = int(np.flatnonzero(surfaces == manifest.anchor_surface_Pa)[0])
    r_anchor = float(geodesy[manifest.anchor_quantity].values[row])
    flattening = float(geodesy["oblateness"].values[row])

    # The refrac boundary: the only degree to radian conversions of the reduction.
    phi_g = np.radians(np.array([float(thermo.attrs["latitude_planetographic_deg"])]))
    tolerance = float(np.radians(manifest.fixed_point_tolerance_deg))
    max_iter = manifest.max_iterations
    rule = manifest.anchor_rule
    tol_m = manifest.convergence_m

    u_of_phi = wind_of_latitude(inputs.wind)

    def march(phi):
        return gd.wind_geoid(np.atleast_1d(np.asarray(phi, dtype="float64")), r_anchor, rule,
                             u_of_phi, *constants, tol_m=tol_m)

    def wind_surface(phi):
        # A full anchored march with these latitudes as nodes; wind_geoid returns the radius
        # at a node exactly, so nothing here interpolates the anchor radius.
        return march(phi).radius

    def zero_wind(phi):
        return np.zeros_like(np.asarray(phi, dtype="float64"))

    def reference_surface(phi):
        radius, _, _ = gd.reference_geoid(np.atleast_1d(np.asarray(phi, dtype="float64")),
                                          r_anchor, *constants, tol_m=tol_m)
        return radius

    solved = {}
    for label, surface, wind in (("wind", wind_surface, u_of_phi),
                                 ("no wind", reference_surface, zero_wind)):
        phi_c, _, iterates, count = lat.planetocentric_fixed_point(
            phi_g, surface, wind, *constants, tol_rad=tolerance, max_iter=max_iter,
            flattening=flattening)
        step = _final_step(iterates)
        if not step < tolerance:
            raise AnchorConvergenceError(
                f"the {label} anchor latitude fixed point took {int(count[0])} of {max_iter} "
                f"iterations and its last step was {np.degrees(step):.3e} deg, above the "
                f"manifest tolerance {manifest.fixed_point_tolerance_deg:g} deg."
            )
        solved[label] = (phi_c, iterates, int(count[0]), step)

    phi_c, iterates, count, step = solved["wind"]
    final = march(phi_c)
    r0 = final.radius
    psi = g_eff_vector(u_of_phi(phi_c), r0, phi_c, *constants)[3]

    nw_phi_c, nw_iterates, _, _ = solved["no wind"]
    nw_r0 = reference_surface(nw_phi_c)
    nw_psi = g_eff_vector(0.0, nw_r0, nw_phi_c, *constants)[3]

    return FrozenAnchor(
        phi_g_rad=float(phi_g[0]),
        phi_c_rad=float(phi_c[0]),
        psi_rad=float(psi[0]),
        r0_m=float(r0[0]),
        u_at_phi_c_ms=float(u_of_phi(phi_c)[0]),
        iterates_rad=np.asarray(iterates[:, 0], dtype="float64"),
        iteration_count=count,
        final_step_rad=step,
        tolerance_rad=tolerance,
        anchor_rule=rule,
        anchor_surface_Pa=manifest.anchor_surface_Pa,
        anchor_quantity=manifest.anchor_quantity,
        r_anchor_m=r_anchor,
        north_start_m=final.north_start_m,
        polar_north_m=final.polar_north_m,
        polar_south_m=final.polar_south_m,
        polar_asymmetry_m=final.polar_asymmetry_m,
        anchor_residual_m=final.anchor_residual_m,
        nowind_phi_c_rad=float(nw_phi_c[0]),
        nowind_psi_rad=float(nw_psi[0]),
        nowind_r0_m=float(nw_r0[0]),
        nowind_iterates_rad=np.asarray(nw_iterates[:, 0], dtype="float64"),
    )
