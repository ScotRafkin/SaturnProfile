"""Newtonian potential and effective gravity for an oblate, rotating, zonally windy planet.

SPEC_01 Step 4. Implements manuscript Eqs. A2, A3, A4, A5 and handoff section 9A.5, in the
convention `lib.schema.HARMONIC_CONVENTION` names as `CASSPIAN-J1`: even zonal harmonics
`J_l` of degree `l`, unscaled, in

    V = -(GM/r) [1 - sum_l J_l (R/r)^l P_l(sin phi_c)]

with `R` the normalization radius of the set, `phi_c` planetocentric, and the Newtonian
gravity positive inward.

Pure functions: every physical input is an argument, nothing is read from a file, and there is
no module level Saturn quantity (SPEC_00 principle 7). Scalars and broadcastable arrays are
both accepted; `r` and `phi_c` broadcast against each other, and the harmonic sum is taken over
a leading axis that is removed before the result is returned.

**The coefficient on the harmonic term of `g` is `l + 1` with `l` the degree: 3 on J2, 5 on J4,
7 on J6.** Reading the handoff's `(2n+1) J_2n` with `n` taken as the degree rather than as the
half degree gives 5, 9, 13, and produces an equatorial gravity of 9.29 rather than 8.95, a four
percent error that looks entirely plausible. The trap is sprung once already in this project
and `data_static/harmonics/null1981.toml` records it. The Step 4 acceptance runs the wrong
coefficient deliberately to show that the check discriminates.

Sign conventions, stated once because two of them are easy to invert:

* `g` and `g_N` are **positive inward**. Lindal's Appendix writes an outward `g_r`; handoff
  section 9A.5 fixes `g_r = -g`.
* `G_phi` is the **latitudinal component of the acceleration in the direction of increasing
  planetocentric latitude**, `-(1/r) dV/dphi`, as SPEC_01 Step 4 writes it. On an oblate planet
  the bulge pulls a mid latitude point toward the equator, so in the northern hemisphere this
  component is **negative**. See `psi` below and the Step 4 report, section 3.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "legendre_even",
    "legendre_even_derivative",
    "potential_V",
    "g_newton",
    "G_phi_newton",
    "omega_abs",
    "g_eff_radial",
    "G_phi_eff",
    "g_eff_vector",
    "MAX_DEGREE",
]

#: SPEC_01 Step 4: even degrees to 12, which covers Null (1981) and Iess et al. (2019).
MAX_DEGREE = 12


def _legendre_all(max_degree: int, x):
    """Return `P_n(x)` and `dP_n/dx` for n = 0 to `max_degree`, as arrays over a leading axis.

    Uses the standard recurrences

        (n + 1) P_{n+1} = (2n + 1) x P_n - n P_{n-1}
        P'_{n+1} = x P'_n + (n + 1) P_n

    The derivative recurrence is used rather than the closed form
    `P'_n = n (x P_n - P_{n-1}) / (x^2 - 1)` because that form is singular at `x = +-1`, which
    is exactly the pole, where the acceptance evaluates.
    """
    x = np.asarray(x, dtype="float64")
    values = np.empty((max_degree + 1,) + x.shape, dtype="float64")
    derivatives = np.empty_like(values)
    values[0] = 1.0
    derivatives[0] = 0.0
    if max_degree >= 1:
        values[1] = x
        derivatives[1] = 1.0
    # A loop over degrees, not over array elements: the SPEC_00 section 3.1 array rule is about
    # element wise Python loops, and this runs at most twelve times whatever the grid size.
    for n in range(1, max_degree):
        values[n + 1] = ((2 * n + 1) * x * values[n] - n * values[n - 1]) / (n + 1)
        derivatives[n + 1] = x * derivatives[n] + (n + 1) * values[n]
    return values, derivatives


def _as_degrees(degrees) -> np.ndarray:
    degrees = np.asarray(degrees, dtype="int64").ravel()
    if degrees.size == 0:
        raise ValueError("at least one harmonic degree is required")
    if np.any(degrees < 2) or np.any(degrees % 2 != 0):
        raise ValueError(
            f"degrees must be even and at least 2; got {degrees.tolist()}. The formulation is "
            "axisymmetric and equatorially symmetric (handoff section 9A.5), so odd degrees "
            "are not carried."
        )
    if np.any(degrees > MAX_DEGREE):
        raise ValueError(f"degrees above {MAX_DEGREE} are not supported; got {degrees.tolist()}")
    return degrees


def legendre_even(degrees, x) -> np.ndarray:
    """`P_l(x)` for the even degrees supplied, stacked on a leading axis."""
    degrees = _as_degrees(degrees)
    values, _ = _legendre_all(int(degrees.max()), x)
    return values[degrees]


def legendre_even_derivative(degrees, x) -> np.ndarray:
    """`dP_l/dx` for the even degrees supplied, stacked on a leading axis."""
    degrees = _as_degrees(degrees)
    _, derivatives = _legendre_all(int(degrees.max()), x)
    return derivatives[degrees]


def _harmonic_terms(r, phi_c, J, degrees, R_norm):
    """Common factors: the ratio `(R/r)^l`, `P_l(sin phi_c)`, and its derivative, broadcast."""
    degrees = _as_degrees(degrees)
    J = np.asarray(J, dtype="float64").ravel()
    if J.shape != degrees.shape:
        raise ValueError(
            f"J has {J.size} values but {degrees.size} degrees were supplied; they must "
            "correspond one to one."
        )
    r = np.asarray(r, dtype="float64")
    phi_c = np.asarray(phi_c, dtype="float64")
    r, phi_c = np.broadcast_arrays(r, phi_c)
    sin_phi = np.sin(phi_c)

    shape = (degrees.size,) + (1,) * r.ndim
    degrees_b = degrees.reshape(shape)
    J_b = J.reshape(shape)
    ratio = (np.asarray(R_norm, dtype="float64") / r) ** degrees_b
    P = legendre_even(degrees, sin_phi)
    dP = legendre_even_derivative(degrees, sin_phi)
    return r, phi_c, sin_phi, degrees_b, J_b, ratio, P, dP


def potential_V(r, phi_c, GM, J, degrees, R_norm):
    """`V = -(GM/r) [1 - sum_l J_l (R/r)^l P_l(sin phi_c)]`. Manuscript Eq. A1, handoff 9A.5."""
    r, _, _, _, J_b, ratio, P, _ = _harmonic_terms(r, phi_c, J, degrees, R_norm)
    total = np.sum(J_b * ratio * P, axis=0)
    return -(np.asarray(GM, dtype="float64") / r) * (1.0 - total)


def g_newton(r, phi_c, GM, J, degrees, R_norm):
    """`g_N = (GM/r^2) [1 - sum_l (l+1) J_l (R/r)^l P_l(sin phi_c)]`, positive inward.

    The coefficient is `l + 1` with `l` the degree: 3 on J2, 5 on J4, 7 on J6.
    """
    r, _, _, degrees_b, J_b, ratio, P, _ = _harmonic_terms(r, phi_c, J, degrees, R_norm)
    total = np.sum((degrees_b + 1) * J_b * ratio * P, axis=0)
    return (np.asarray(GM, dtype="float64") / r**2) * (1.0 - total)


def G_phi_newton(r, phi_c, GM, J, degrees, R_norm):
    """`-(1/r) dV/dphi = -(GM/r^2) sum_l J_l (R/r)^l cos(phi_c) P'_l(sin phi_c)`.

    The component along increasing planetocentric latitude. Negative in the northern
    hemisphere for an oblate planet, the bulge pulling a mid latitude point equatorward.
    """
    r, phi_c, _, _, J_b, ratio, _, dP = _harmonic_terms(r, phi_c, J, degrees, R_norm)
    total = np.sum(J_b * ratio * np.cos(phi_c) * dP, axis=0)
    return -(np.asarray(GM, dtype="float64") / r**2) * total


def omega_abs(u, r, phi_c, Omega):
    """`Omega_abs = Omega + u / (r cos phi_c)`. Manuscript Eq. A2.

    The absolute angular rate of a parcel carrying zonal wind `u` eastward in a frame rotating
    at `Omega`.

    The guard on `cos phi_c` catches an **exact** zero only. `np.cos(np.pi/2)` is 6.1e-17, not
    zero, so a caller passing `np.radians(90.0)` takes the division path. With `u = 0` at the
    pole, which is what a zonal wind field must supply there, that is `0 / 6e-17 = 0` and the
    result is correct. With a nonzero `u` at a numerically near polar latitude it would give an
    enormous rate rather than an error. Do not rely on this guard to make a nonzero polar wind
    safe; a wind field that is nonzero at the pole is the defect, and the tool that writes it
    is where that is caught. Noted at the Step 4 review.
    """
    u = np.asarray(u, dtype="float64")
    r = np.asarray(r, dtype="float64")
    phi_c = np.asarray(phi_c, dtype="float64")
    cos_phi = np.cos(phi_c)
    with np.errstate(divide="ignore", invalid="ignore"):
        contribution = np.where(cos_phi != 0.0, u / (r * cos_phi), 0.0)
    return np.asarray(Omega, dtype="float64") + contribution


def g_eff_radial(u, r, phi_c, Omega, GM, J, degrees, R_norm):
    """Effective gravity, radial component, positive inward. Manuscript Eq. A3, equally A4.

    `g = g_N - Omega_abs^2 r cos^2(phi_c)`. The centrifugal term is written by Lindal as
    `(2/3) w^2 r (1 - P_2(sin phi_c))`, which is the same quantity, since
    `1 - P_2(sin phi) = (3/2) cos^2(phi)`.
    """
    r_b = np.asarray(r, dtype="float64")
    phi_b = np.asarray(phi_c, dtype="float64")
    w = omega_abs(u, r_b, phi_b, Omega)
    centrifugal = w**2 * r_b * np.cos(phi_b) ** 2
    return g_newton(r_b, phi_b, GM, J, degrees, R_norm) - centrifugal


def G_phi_eff(u, r, phi_c, Omega, GM, J, degrees, R_norm):
    """Effective gravity, latitudinal component. Manuscript Eq. A5.

    `G_phi = G_phi_N - Omega_abs^2 r cos(phi_c) sin(phi_c)`, along increasing latitude. The
    centrifugal contribution is the latitudinal projection of the outward, axis perpendicular
    centrifugal acceleration, and like the Newtonian term it points toward the equator.
    """
    r_b = np.asarray(r, dtype="float64")
    phi_b = np.asarray(phi_c, dtype="float64")
    w = omega_abs(u, r_b, phi_b, Omega)
    centrifugal = w**2 * r_b * np.cos(phi_b) * np.sin(phi_b)
    return G_phi_newton(r_b, phi_b, GM, J, degrees, R_norm) - centrifugal


def g_eff_vector(u, r, phi_c, Omega, GM, J, degrees, R_norm):
    """Both components of the effective gravity, its magnitude, and the angle `psi`.

    Returns `(g, G_phi, magnitude, psi)` with `g` positive inward, `G_phi` along increasing
    latitude, `magnitude = hypot(g, G_phi)`, and

        psi = arctan(-G_phi / g)

    `psi` is the tilt of the local vertical from the radial direction toward the pole, positive
    in the northern hemisphere, so that the planetographic to planetocentric relation of Step 6
    reads `phi_c = phi_g - psi`. The minus sign on `G_phi` is part of that definition, not a
    correction: with `G_phi` along increasing latitude it is negative in the north, and `-G_phi`
    is the equatorward magnitude the tilt is measured from. This is Lindal Eq. 5 read with
    handoff section 9A.5, where Lindal's outward `g_r` is `-g`.

    `arctan2` rather than `arctan` so that the southern hemisphere needs no case split: there
    `G_phi` is positive and `psi` comes out negative, which is what `phi_c = phi_g - psi`
    requires. Convention ruled at the Step 4 review; SPEC_01 v0.6 Step 4 states it.
    """
    g = g_eff_radial(u, r, phi_c, Omega, GM, J, degrees, R_norm)
    G_phi = G_phi_eff(u, r, phi_c, Omega, GM, J, degrees, R_norm)
    magnitude = np.hypot(g, G_phi)
    psi = np.arctan2(-G_phi, g)
    return g, G_phi, magnitude, psi
