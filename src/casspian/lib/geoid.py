"""The reference surface: Lindal Appendix Eqs. 11 to 14, and manuscript Eq. B3.

SPEC_01 Step 5. Pure functions, NumPy in and out, no file access, no module level Saturn
quantity. Every physical input is an argument, including the harmonic set, the rotation rate,
and the wind.

Two surfaces are built here and they are not the same object.

`reference_geoid` is the **no wind** surface: the equipotential of a uniformly rotating fluid
that passes through the stated polar radius. It is found by Newton iteration on `U`, which is
Lindal Eqs. 12 to 14, and it is an equipotential by construction.

`wind_geoid` is the surface a **latitude dependent** zonal wind produces. With `u` varying in
latitude the effective field is not exactly conservative, which is the whole point of
manuscript section A6, so there is no potential whose level set it is. It is therefore obtained
by integrating the slope equation of Eq. B3,

    g dr0/dphi = r0 G_phi,

from the pole inward with a fourth order Runge Kutta scheme, and the residual
`U(r0(phi), phi) - U_ref` is returned alongside as a diagnostic of how far the wind has taken
the surface from an equipotential. Both are returned because the size of that residual is
itself a result.

Sign conventions follow `lib.gravity`. `g` is positive inward and `G_phi` is the component
along increasing planetocentric latitude, negative in the northern hemisphere. With those,
`dU/dr = g` and `(1/r) dU/dphi = -G_phi`, so Eq. B3 needs no hemisphere dependent sign: north
of the equator `G_phi < 0` gives `dr0/dphi < 0` and the surface falls away from the equatorial
bulge, and south of it both signs reverse together.
"""

from __future__ import annotations

import numpy as np

from casspian.lib.gravity import G_phi_eff, g_eff_radial, potential_V

__all__ = ["U_rigid", "reference_geoid", "wind_geoid", "radius_at", "ellipsoid_seed"]


def U_rigid(r, phi_c, Omega, GM, J, degrees, R_norm):
    """Potential of a uniformly rotating fluid, `U = V - (1/2) Omega^2 r^2 cos^2(phi_c)`.

    Lindal Appendix Eq. 11, **term for term and in the same sign**. His Eqs. 9 and 11 give
    `U = -GM/r + (GM/r) sum_i J_2i (R/r)^2i P_2i(sin phi) - (1/2) omega^2 r^2 cos^2(phi)` with
    `g = -grad U`, which is exactly `V - (1/2) Omega^2 r^2 cos^2(phi_c)` with the `V` of
    `lib.gravity`. SPEC_01 v0.6 claimed his form carried the opposite overall sign; it does
    not, and v0.7 withdraws the claim.

    In this convention `g_eff = -grad U`, so `dU/dr = g` with `g` positive inward, and `U`
    tends to zero on the rotation axis at infinity, where both terms vanish.
    """
    r = np.asarray(r, dtype="float64")
    phi_c = np.asarray(phi_c, dtype="float64")
    V = potential_V(r, phi_c, GM, J, degrees, R_norm)
    centrifugal = 0.5 * np.asarray(Omega, dtype="float64") ** 2 * r**2 * np.cos(phi_c) ** 2
    return V - centrifugal


def ellipsoid_seed(phi_c, r_polar, Omega, GM, J, degrees):
    """A first order hydrostatic ellipsoid, used only to start the Newton iteration.

    The flattening is estimated as `f = (3/2) J2 + q/2` with `q = Omega^2 r_polar^3 / GM`,
    which for Saturn gives 0.082 against the true 0.098. That is close enough to start from and
    is never used as a result.
    """
    phi_c = np.asarray(phi_c, dtype="float64")
    degrees = np.asarray(degrees).ravel()
    J = np.asarray(J, dtype="float64").ravel()
    J2 = float(J[degrees == 2][0]) if np.any(degrees == 2) else 0.0
    q = float(np.asarray(Omega) ** 2 * r_polar**3 / np.asarray(GM))
    flattening = 1.5 * J2 + 0.5 * q
    a = r_polar / (1.0 - flattening)
    return a * r_polar / np.sqrt(
        (a * np.sin(phi_c)) ** 2 + (r_polar * np.cos(phi_c)) ** 2
    )


def reference_geoid(
    phi_c_grid, r_polar, Omega, GM, J, degrees, R_norm, tol_m=1.0, max_iter=50
):
    """The no wind reference geoid through `r_polar`. Lindal Appendix Eqs. 12 to 14.

    `U_ref = U_rigid(r_polar, pi/2)`, and at each latitude the radius solving
    `U(r, phi) = U_ref` is found by Newton iteration from an ellipsoid seed.

    SPEC_01 Step 5 writes the update as `r <- r + (U - U_ref) / g_r`, where `g_r` is Lindal's
    **outward** radial gravity. Handoff section 9A.5 fixes `g_r = -g`, and `dU/dr = g` in the
    convention of `U_rigid`, so the Newton step in the convention of this package is
    `r <- r - (U - U_ref) / g`. The two are the same update written twice.

    Returns `(r_ref, iterations, residual_U)`: the radius at each latitude, the number of
    Newton steps that latitude took, and `|U(r_ref, phi) - U_ref|` there.
    """
    phi = np.asarray(phi_c_grid, dtype="float64")
    U_ref = float(U_rigid(r_polar, np.pi / 2, Omega, GM, J, degrees, R_norm))

    r = np.array(ellipsoid_seed(phi, r_polar, Omega, GM, J, degrees), dtype="float64")
    r = np.broadcast_to(r, phi.shape).astype("float64").copy()
    iterations = np.zeros(phi.shape, dtype="int64")
    active = np.ones(phi.shape, dtype=bool)

    for _ in range(int(max_iter)):
        if not active.any():
            break
        U = U_rigid(r[active], phi[active], Omega, GM, J, degrees, R_norm)
        g = g_eff_radial(0.0, r[active], phi[active], Omega, GM, J, degrees, R_norm)
        step = -(U - U_ref) / g
        r[active] = r[active] + step
        iterations[active] += 1
        # A latitude stops when its own step falls below the tolerance, so a slowly
        # converging latitude does not make the others iterate needlessly.
        still = np.abs(step) >= float(tol_m)
        idx = np.flatnonzero(active)
        active[idx[~still]] = False

    residual = np.abs(U_rigid(r, phi, Omega, GM, J, degrees, R_norm) - U_ref)
    return r, iterations, residual


def _rk4_slope(phi, r, u_of_phi, Omega, GM, J, degrees, R_norm):
    """`dr0/dphi = r0 G_phi / g`, manuscript Eq. B3, with the Step 4 components."""
    u = u_of_phi(phi)
    g = g_eff_radial(u, r, phi, Omega, GM, J, degrees, R_norm)
    G_phi = G_phi_eff(u, r, phi, Omega, GM, J, degrees, R_norm)
    return r * G_phi / g


def _integrate_hemisphere(nodes, r_start, u_of_phi, Omega, GM, J, degrees, R_norm):
    """RK4 from `nodes[0]` (a pole) along `nodes`, returning the radius at each node."""
    radii = np.empty(nodes.shape, dtype="float64")
    radii[0] = r_start
    r = float(r_start)
    # A loop over grid nodes, which an ODE march requires: each step needs the previous
    # result. The SPEC_00 section 3.1 array rule is about element wise arithmetic that could
    # be whole array, and an initial value problem is not that.
    for i in range(nodes.size - 1):
        phi0, phi1 = float(nodes[i]), float(nodes[i + 1])
        h = phi1 - phi0
        args = (u_of_phi, Omega, GM, J, degrees, R_norm)
        k1 = float(_rk4_slope(phi0, r, *args))
        k2 = float(_rk4_slope(phi0 + 0.5 * h, r + 0.5 * h * k1, *args))
        k3 = float(_rk4_slope(phi0 + 0.5 * h, r + 0.5 * h * k2, *args))
        k4 = float(_rk4_slope(phi1, r + h * k3, *args))
        r = r + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        radii[i + 1] = r
    return radii


def wind_geoid(phi_c_grid, r_polar, u_of_phi, Omega, GM, J, degrees, R_norm):
    """The surface a latitude dependent zonal wind produces. Manuscript Eq. B3.

    `u_of_phi` is a callable taking planetocentric latitude in radians and returning the zonal
    wind in m/s. It must return zero at the pole: a nonzero zonal wind there is not a wind
    field but a defect, and `lib.gravity.omega_abs` will not rescue it (see its docstring).

    Integrates `g dr0/dphi = r0 G_phi` from each pole inward with RK4, anchored at `r_polar`.
    Each hemisphere is integrated from its own pole, so a wind that is not symmetric about the
    equator is handled without a reflection assumption. The pole is inserted into the march if
    the supplied grid does not contain it; the returned arrays are on the supplied grid.

    Returns `(r0, closure)`: the radius at each supplied latitude, and
    `U(r0(phi), phi) - U_ref` there against the **no wind** potential through the same polar
    radius. That is the dynamical height of the wind surface above the no wind geoid in
    potential units; divide by `g` for meters.

    No pseudo-potential built from `Omega_abs` is evaluated. With `u` varying in latitude the
    gradient of such a quantity is not the effective gravity, so its variation along the
    surface restates the wind kinetic term and measures nothing about conservativeness. The
    non conservativeness of the field is the shear kernel of Eq. A15 and belongs to the forward
    model. Ruled at the Step 5 review; SPEC_01 v0.7 states it.

    There is no `tol_m` or `max_iter`: an initial value problem does not iterate.
    """
    phi = np.asarray(phi_c_grid, dtype="float64")
    pole = np.pi / 2
    result = np.empty(phi.shape, dtype="float64")

    for sign in (+1, -1):
        selected = phi >= 0 if sign > 0 else phi < 0
        if not selected.any():
            continue
        here = phi[selected]
        order = np.argsort(-sign * here)  # from the pole inward
        march = here[order]
        if abs(abs(march[0]) - pole) > 1e-12:
            march = np.concatenate(([sign * pole], march))
            radii = _integrate_hemisphere(
                march, r_polar, u_of_phi, Omega, GM, J, degrees, R_norm
            )[1:]
        else:
            radii = _integrate_hemisphere(
                march, r_polar, u_of_phi, Omega, GM, J, degrees, R_norm
            )
        restored = np.empty_like(radii)
        restored[order] = radii
        result[selected] = restored

    U_ref = float(U_rigid(r_polar, pole, Omega, GM, J, degrees, R_norm))
    closure = U_rigid(result, phi, Omega, GM, J, degrees, R_norm) - U_ref
    return result, closure


def radius_at(phi_c, phi_c_grid, radii):
    """Interpolate a constructed geoid to a latitude.

    **`1/r^2` linear in `sin^2(phi_c)`**, between the target's two neighboring grid nodes.
    SPEC_01 v0.7 Step 5. For an ellipse `1/r^2 = sin^2(phi)/b^2 + cos^2(phi)/a^2` holds
    exactly, so the scheme interpolates only the harmonic and wind residual and is about eight
    times more accurate at equal spacing than linear in `sin(phi_c)`, which v0.6 specified.

    **Not for the frozen anchor radius.** `reference_geoid` solves every latitude
    independently and `wind_geoid` marches any node set, so a latitude at which the radius is
    actually wanted is put in the grid and solved or marched. This function exists for the
    vectorized uses of Step 7 and later, where a whole binned latitude axis is wanted at once.

    The bracketing is done in `phi_c`, not in `sin^2(phi_c)`, because the latter is not
    monotonic across the equator: a global search on it would fold the hemispheres onto each
    other and silently return a northern radius for a southern latitude on any surface that is
    not symmetric. An interval that straddles the equator is refused for the same reason, so a
    grid is to contain the equator as a node. Refuses to extrapolate.
    """
    phi_c = np.asarray(phi_c, dtype="float64")
    grid = np.asarray(phi_c_grid, dtype="float64")
    radii = np.asarray(radii, dtype="float64")
    if grid.ndim != 1 or grid.size < 2:
        raise ValueError("the latitude grid must be one dimensional with at least two nodes")
    order = np.argsort(grid)
    grid, radii = grid[order], radii[order]

    tol = 1e-12
    if np.any(phi_c < grid[0] - tol) or np.any(phi_c > grid[-1] + tol):
        raise ValueError(
            f"latitude outside the constructed grid [{np.degrees(grid[0]):.4f}, "
            f"{np.degrees(grid[-1]):.4f}] degrees; this function does not extrapolate."
        )

    index = np.clip(np.searchsorted(grid, phi_c, side="right") - 1, 0, grid.size - 2)
    lo, hi = grid[index], grid[index + 1]

    straddles = (lo < -tol) & (hi > tol)
    if np.any(straddles):
        bad = np.atleast_1d(np.degrees(phi_c))[np.atleast_1d(straddles)]
        raise ValueError(
            f"the interval bracketing {bad.tolist()} degrees straddles the equator, where "
            "sin^2(phi_c) is not monotonic. Put the equator in the grid as a node."
        )

    sin2 = np.sin(phi_c) ** 2
    sin2_lo, sin2_hi = np.sin(lo) ** 2, np.sin(hi) ** 2
    inverse_lo, inverse_hi = 1.0 / radii[index] ** 2, 1.0 / radii[index + 1] ** 2
    span = sin2_hi - sin2_lo
    # `span` is nonzero once the straddling case is refused: two distinct nodes on the same
    # side of the equator cannot share a value of sin^2. The guard keeps the divide quiet if
    # a caller supplies duplicate nodes.
    weight = np.where(span != 0.0, (sin2 - sin2_lo) / np.where(span != 0.0, span, 1.0), 0.0)
    return (inverse_lo + weight * (inverse_hi - inverse_lo)) ** -0.5
