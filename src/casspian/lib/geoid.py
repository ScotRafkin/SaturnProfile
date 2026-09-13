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

as **one** march from the north pole through the equator to the south pole with a fourth order
Runge Kutta scheme. Eq. B3 is first order and carries one constant, so there is one surface and
one constant to fix; marching each hemisphere from its own pole with the same polar radius
gives two surfaces that do not meet where the wind is not symmetric. The constant is set by a
declared `anchor_rule`, and the residual `U(r0(phi), phi) - U_ref` is returned alongside as the
dynamical height above the no wind geoid.

Sign conventions follow `lib.gravity`. `g` is positive inward and `G_phi` is the component
along increasing planetocentric latitude, negative in the northern hemisphere. With those,
`dU/dr = g` and `(1/r) dU/dphi = -G_phi`, so Eq. B3 needs no hemisphere dependent sign: north
of the equator `G_phi < 0` gives `dr0/dphi < 0` and the surface falls away from the equatorial
bulge, and south of it both signs reverse together.
"""

from __future__ import annotations

from dataclasses import dataclass

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
    phi_c_grid, r_polar, Omega, GM, J, degrees, R_norm, tol_m=1.0, max_iter=50,
    anchor_latitude=np.pi / 2,
):
    """The no wind reference geoid through `r_polar`. Lindal Appendix Eqs. 12 to 14.

    `U_ref = U_rigid(r_polar, pi/2)`, and at each latitude the radius solving
    `U(r, phi) = U_ref` is found by Newton iteration from an ellipsoid seed.

    SPEC_01 Step 5 writes the update as `r <- r + (U - U_ref) / g_r`, where `g_r` is Lindal's
    **outward** radial gravity. Handoff section 9A.5 fixes `g_r = -g`, and `dU/dr = g` in the
    convention of `U_rigid`, so the Newton step in the convention of this package is
    `r <- r - (U - U_ref) / g`. The two are the same update written twice.

    `anchor_latitude` (radians, default the pole) is where the surface passes through the
    given radius, so `U_ref = U_rigid(r_polar, anchor_latitude)`. The default keeps the polar
    anchor of Lindal Eq. 12; `0.0` anchors the no wind surface on an equatorial radius, the
    no wind counterpart of the `equatorial_radius` rule of `wind_geoid` (SPEC_02 v0.8 Step 6).
    The argument keeps its name `r_polar` because the default is the pole.

    Returns `(r_ref, iterations, residual_U)`: the radius at each latitude, the number of
    Newton steps that latitude took, and `|U(r_ref, phi) - U_ref|` there.
    """
    phi = np.asarray(phi_c_grid, dtype="float64")
    U_ref = float(U_rigid(r_polar, anchor_latitude, Omega, GM, J, degrees, R_norm))

    # The seed ellipsoid is parameterized by its polar radius. Away from the pole that radius
    # is estimated by scaling, so the seed passes near the anchor; Newton does the rest.
    seed_polar = float(r_polar)
    if float(anchor_latitude) != np.pi / 2:
        at_anchor = float(ellipsoid_seed(np.array(float(anchor_latitude)), seed_polar, Omega,
                                         GM, J, degrees))
        seed_polar = seed_polar * seed_polar / at_anchor
    r = np.array(ellipsoid_seed(phi, seed_polar, Omega, GM, J, degrees), dtype="float64")
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


def _march_nodes(nodes, r_start, u_of_phi, Omega, GM, J, degrees, R_norm):
    """RK4 from `nodes[0]` along `nodes`, returning the radius at each node.

    The nodes run north pole to south pole, so this is the whole surface in one pass. It was
    named for a hemisphere when the march was per hemisphere; that was the error SPEC_01 v0.16
    corrects.
    """
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


@dataclass
class WindGeoidResult:
    """The marched wind surface and the numbers that describe how it was anchored."""

    radius: np.ndarray
    closure: np.ndarray
    polar_north_m: float
    polar_south_m: float
    polar_asymmetry_m: float
    anchor_residual_m: float
    north_start_m: float
    anchor_rule: str
    #: The marched radius at the exact equator node, whatever the rule.
    equator_radius_m: float = float("nan")
    #: The node the anchor was read at: 0.0 for `equatorial_radius`, the poles for the pole
    #: rules, `anchor_latitude` for `latitude`, None for `mean_polar_radius` (two nodes).
    anchor_node_latitude_rad: float | None = None

    def __iter__(self):
        """Unpack as `(radius, closure)`, which is what most callers want."""
        return iter((self.radius, self.closure))


def wind_geoid(phi_c_grid, r_anchor, anchor_rule, u_of_phi, Omega, GM, J, degrees, R_norm,
               anchor_latitude=None, tol_m=1.0e-6, max_iter=40, march_step_deg=0.05):
    """The surface a latitude dependent zonal wind produces. Manuscript Eq. B3.

    `u_of_phi` is a callable taking planetocentric latitude in radians and returning the zonal
    wind in m/s. It must return zero at both poles: a nonzero zonal wind there is not a wind
    field but a defect, SPEC_00 section 6.6 refuses a kind W file that carries one, and
    `lib.gravity.omega_abs` will not rescue it (see its docstring). It must also return an
    array shaped like its argument; under NumPy 2 a one element array no longer converts to a
    scalar, and a callable that always returns an array breaks the march.

    **One march, one constant.** Eq. B3, `g dr0/dphi = r0 G_phi`, is first order and carries a
    single constant of integration, so the surface is integrated once from the north pole
    through the equator to the south pole. Marching from each pole separately with the same
    polar radius produces two surfaces which, for a wind that is not symmetric about the
    equator, do not meet: the Step 7 v0.15 figure showed a 38 km step at the equator from
    exactly that mistake. The constant is fixed by `anchor_rule`:

    * `"mean_polar_radius"` (the default): the north polar start is found so that the two
      polar radii of the marched surface average to `r_anchor`. Lindal states a **mean** polar
      radius, which is what this reproduces.
    * `"north_pole"` or `"south_pole"`: `r_anchor` is the radius at that pole.
    * `"latitude"`: `r_anchor` is the radius at `anchor_latitude` (radians), for a surface
      anchored on an observed radius.
    * `"equatorial_radius"`: `r_anchor` is the radius at planetocentric latitude zero, the
      `latitude` rule at the equator under its own name, so a manifest needs no latitude key
      (SPEC_01 v0.20 Step 5 amendment).

    The root find is a secant on the north polar start. The map from start to outcome is very
    nearly affine, because Eq. B3 is linear in `r0` to the accuracy that matters here, so it
    converges in a few steps; the residual it reached is returned rather than assumed.

    Returns a `WindGeoidResult`, which unpacks as `(radius, closure)` for callers that want
    only those. `closure` is `U(r0(phi), phi) - U_ref` against the **no wind** potential
    through `r_anchor`, the dynamical height of the wind surface above the no wind geoid in
    potential units; divide by `g` for meters. No pseudo-potential built from `Omega_abs` is
    evaluated: with `u` varying in latitude its gradient is not the effective gravity, so its
    variation along the surface restates the wind kinetic term and measures nothing.
    """
    phi = np.asarray(phi_c_grid, dtype="float64")
    pole = np.pi / 2

    # The march runs on its own dense grid spanning both poles, unioned with whatever the
    # caller asked for, so every requested latitude is a node and the accuracy of the march
    # never depends on the caller's grid. A caller who passes a single hemisphere would
    # otherwise get one enormous step across the other one, which is silent and wrong.
    dense = np.radians(np.arange(-90.0, 90.0 + 0.5 * march_step_deg, march_step_deg))
    # The equator is always an exact node, 0.0 and not the floating point neighbor the dense
    # grid lands on, so a radius read there is marched and never interpolated (SPEC_01 v0.20).
    dense = dense[np.abs(dense) > 1.0e-9]
    wanted = [phi, dense, [-pole, 0.0, pole]]
    if anchor_rule == "latitude":
        if anchor_latitude is None:
            raise ValueError("anchor_rule 'latitude' needs anchor_latitude in radians")
        # The anchor latitude is a node, so the radius the anchor is set from is marched and
        # never interpolated. Without this it would be the one latitude in the function that
        # broke that rule.
        wanted.append([float(anchor_latitude)])
    nodes = np.unique(np.concatenate(wanted))
    order = np.argsort(-nodes)          # north pole first, south pole last
    march_nodes = nodes[order]
    equator_index = int(np.flatnonzero(march_nodes == 0.0)[0])

    def march(start):
        return _march_nodes(march_nodes, start, u_of_phi, Omega, GM, J, degrees,
                                     R_norm)

    def outcome(radii):
        if anchor_rule == "mean_polar_radius":
            return 0.5 * (float(radii[0]) + float(radii[-1]))
        if anchor_rule == "north_pole":
            return float(radii[0])
        if anchor_rule == "south_pole":
            return float(radii[-1])
        if anchor_rule == "latitude":
            if anchor_latitude is None:
                raise ValueError("anchor_rule 'latitude' needs anchor_latitude in radians")
            return float(np.interp(float(anchor_latitude), march_nodes[::-1], radii[::-1]))
        if anchor_rule == "equatorial_radius":
            return float(radii[equator_index])
        raise ValueError(
            f"unknown anchor_rule {anchor_rule!r}; SPEC_01 v0.20 Step 5 defines "
            "'mean_polar_radius', 'north_pole', 'south_pole', 'latitude' and 'equatorial_radius'"
        )

    start = float(r_anchor)
    radii = march(start)
    residual = outcome(radii) - float(r_anchor)
    if anchor_rule != "north_pole":
        previous_start, previous_residual = start, residual
        # The first correction scales the start rather than shifting it. Eq. B3 is nearly
        # homogeneous of degree one in r0, so every radius of the march scales with the start;
        # shifting by the residual assumes the anchored radius moves one for one with the north
        # polar start, which holds at a pole and is about 10 percent wrong at the equator
        # (SPEC_02 v0.8 Step 6). The secant then finishes from a start already close.
        start = start * float(r_anchor) / outcome(radii)
        for _ in range(int(max_iter)):
            radii = march(start)
            residual = outcome(radii) - float(r_anchor)
            if abs(residual) <= float(tol_m):
                break
            slope = (residual - previous_residual) / (start - previous_start)
            if slope == 0.0:
                break
            previous_start, previous_residual = start, residual
            start = start - residual / slope

    U_ref = float(U_rigid(r_anchor, pole, Omega, GM, J, degrees, R_norm))
    on_grid = np.interp(phi, march_nodes[::-1], radii[::-1])
    closure = U_rigid(on_grid, phi, Omega, GM, J, degrees, R_norm) - U_ref
    north, south = float(radii[0]), float(radii[-1])
    return WindGeoidResult(
        radius=on_grid,
        closure=closure,
        polar_north_m=north,
        polar_south_m=south,
        polar_asymmetry_m=south - north,
        anchor_residual_m=float(residual),
        north_start_m=float(start),
        anchor_rule=anchor_rule,
        equator_radius_m=float(radii[equator_index]),
        anchor_node_latitude_rad={"equatorial_radius": 0.0, "north_pole": pole,
                                  "south_pole": -pole,
                                  "latitude": anchor_latitude}.get(anchor_rule),
    )


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
