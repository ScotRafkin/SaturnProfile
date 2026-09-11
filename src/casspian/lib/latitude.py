"""Planetographic to planetocentric conversion. Handoff section 9A.2, Lindal Eq. 5.

SPEC_01 Step 6. Pure functions, NumPy in and out, no file access, no module level Saturn
quantity.

The two latitudes differ because the local vertical on an oblate planet is not radial. The
planetographic latitude is the angle the local vertical makes with the equatorial plane; the
planetocentric latitude is the angle the radius vector makes with it. They differ by `psi`, the
tilt of the local vertical from the radial direction, which `lib.gravity.g_eff_vector` returns
as `arctan2(-G_phi, g)`, positive in the north:

    phi_g = phi_c + psi(phi_c)

`psi` depends on where you stand, so recovering `phi_c` from `phi_g` is a fixed point:
`phi_c <- phi_g - psi(phi_c)`, iterated from an ellipsoid seed. It converges quickly because
`psi` varies slowly with latitude.

**Angles are radians**, in and out, with no exceptions, consistent with `lib.gravity` and
`lib.geoid`. That includes the convergence tolerance: SPEC_01 v0.7 named it `tol_deg`, the only
degree valued quantity anywhere in `lib`, and v0.8 corrected it to `tol_rad` at the Step 6
review. A caller who thinks in degrees converts at the call site.

**`surface` is a callable, and it must solve, not interpolate.** It takes planetocentric
latitudes in radians and returns the radius of the constructed surface there. For the no wind
geoid that is a direct Newton solve at those latitudes; for the wind geoid it is a march with
those latitudes inserted as nodes. SPEC_01 v0.7 forbids passing an interpolated radius here,
because this is where the frozen anchor radius is decided.
"""

from __future__ import annotations

import numpy as np

from casspian.lib.gravity import g_eff_vector

__all__ = [
    "planetocentric_from_ellipsoid",
    "planetocentric_fixed_point",
    "planetographic_from_planetocentric",
]


def planetocentric_from_ellipsoid(phi_g, flattening):
    """The ellipsoid relation `tan(phi_c) = (1 - f)^2 tan(phi_g)`, used as the seed.

    Exact for an ellipsoid of flattening `f` and a good starting point for the fixed point on
    a real geoid. `arctan2` is used rather than `arctan` so that the poles, where `tan(phi_g)`
    is infinite, return exactly `+-pi/2` instead of a NaN.
    """
    phi_g = np.asarray(phi_g, dtype="float64")
    scale = (1.0 - np.asarray(flattening, dtype="float64")) ** 2
    return np.arctan2(scale * np.sin(phi_g), np.cos(phi_g))


def planetographic_from_planetocentric(
    phi_c, surface, u_of_phi, Omega, GM, J, degrees, R_norm
):
    """`phi_g = phi_c + psi(phi_c)`. A direct evaluation, no iteration.

    This is the inverse of `planetocentric_fixed_point` and the reason that function needs to
    iterate at all: `psi` is known as a function of `phi_c`, not of `phi_g`.
    """
    phi_c = np.asarray(phi_c, dtype="float64")
    r = np.asarray(surface(phi_c), dtype="float64")
    _, _, _, psi = g_eff_vector(
        u_of_phi(phi_c), r, phi_c, Omega, GM, J, degrees, R_norm
    )
    return phi_c + psi


def planetocentric_fixed_point(
    phi_g,
    surface,
    u_of_phi,
    Omega,
    GM,
    J,
    degrees,
    R_norm,
    tol_rad=1.0e-8,
    max_iter=20,
    flattening=None,
):
    """Solve `phi_c = phi_g - psi(phi_c)` by fixed point iteration from an ellipsoid seed.

    `flattening` seeds the iteration and is optional: when it is not given, the surface itself
    is asked for its equatorial and polar radii and the flattening is taken from those, so the
    seed comes from the same object the iteration converges on. The converged result does not
    depend on the seed; only the iteration count does.

    Vectorized over `phi_g`, as the wind tool needs. Each element stops on its own tolerance.

    Returns `(phi_c, psi, iterates, count)`: the converged planetocentric latitude in radians,
    the tilt, the full history of iterates as an array of shape `(n_steps + 1,) + phi_g.shape`
    with the seed first, and the number of steps each element took.

    **The returned `psi` is the tilt at the previous iterate, not at the returned `phi_c`.**
    That is deliberate and is the more useful of the two: `phi_c + psi` equals `phi_g`
    **exactly**, because the last update was `phi_c = phi_g - psi`. The tilt evaluated afresh
    at the returned `phi_c` differs from it by less than the final step times `dpsi/dphi`,
    which is under a tenth of the tolerance. The round trip check of the Step 6 acceptance
    agrees to 2.3e-11 degrees. Do not "fix" this into an extra evaluation of `psi`: it would
    cost a surface solve and would break the exact identity.
    """
    phi_g = np.asarray(phi_g, dtype="float64")
    tol = float(tol_rad)

    if flattening is None:
        equator = float(np.asarray(surface(np.array([0.0])))[0])
        pole = float(np.asarray(surface(np.array([np.pi / 2])))[0])
        flattening = (equator - pole) / equator

    phi_c = np.array(
        np.broadcast_to(planetocentric_from_ellipsoid(phi_g, flattening), phi_g.shape),
        dtype="float64",
    )
    iterates = [phi_c.copy()]
    count = np.zeros(phi_g.shape, dtype="int64")
    active = np.ones(phi_g.shape, dtype=bool)
    psi = np.zeros(phi_g.shape, dtype="float64")

    for _ in range(int(max_iter)):
        if not active.any():
            break
        here = phi_c[active]
        r = np.asarray(surface(here), dtype="float64")
        _, _, _, psi_here = g_eff_vector(
            u_of_phi(here), r, here, Omega, GM, J, degrees, R_norm
        )
        psi[active] = psi_here
        updated = phi_g[active] - psi_here
        step = updated - here
        phi_c[active] = updated
        count[active] += 1
        iterates.append(phi_c.copy())
        idx = np.flatnonzero(active)
        active[idx[np.abs(step) < tol]] = False

    return phi_c, psi, np.array(iterates), count
