"""Pressure and temperature from refractivity by hydrostatic integration. Manuscript Eqs. B4 to B6.

SPEC_03 Step 2. Pure functions, NumPy in and out, no file access, no module level Saturn
quantity. The constants of nature come from `lib.constants`.

**The layer integral is exact for a log-linear density** (SPEC_00 section 7.3 rule i). Density
is close to exponential in geopotential across a layer, so the mass per unit area between levels
`k` and `k + 1` is the exact integral of the interpolant that is linear in `ln rho`:

    I_{k+1/2} = integral from Phi_{k+1} to Phi_k of rho dPhi
              = (rho_{k+1} - rho_k) (Phi_k - Phi_{k+1}) / ln(rho_{k+1} / rho_k)

never the trapezoid of the endpoint densities, which at Table I's spacing is in error by about
1e-2 and would hide the local vertical correction of SPEC_03 Step 0.

**Direction.** Levels are top down, as kind N stores them: `Phi` decreases with the level index,
and the pressure is built from the top by `p_0 = p_b`, `p_{k+1} = p_k + I_{k+1/2}` (B7.2), the only
direction this module integrates.
"""

from __future__ import annotations

import numpy as np

from casspian.lib.constants import AVOGADRO_CONSTANT, BOLTZMANN_CONSTANT

__all__ = [
    "LIMIT_LOG_RATIO",
    "density",
    "layer_mass",
    "pressure_from_top",
    "temperature",
]

#: SPEC_03 Step 2: below this magnitude of `ln(rho_{k+1} / rho_k)` the layer integral is taken in
#: its limit `rho_k (Phi_k - Phi_{k+1})`, which differs from the exact form by half the log ratio
#: relative, under 5e-11.
LIMIT_LOG_RATIO = 1.0e-10


def density(N, R_bar, m_bar_kg_mol):
    """`rho = (N / R_bar) m_bar`, Eq. B4 with the mass, in kg m-3.

    `m_bar_kg_mol` is the mean molar mass in kg/mol, as kind N and kind C carry it, converted to
    kg per molecule with the Avogadro constant.
    """
    N = np.asarray(N, dtype="float64")
    R_bar = np.asarray(R_bar, dtype="float64")
    m = np.asarray(m_bar_kg_mol, dtype="float64") / AVOGADRO_CONSTANT
    return N / R_bar * m


def layer_mass(rho, Phi):
    """`I_{k+1/2}`, the exact integral of the log-linear density across each layer, in Pa.

    `rho` and `Phi` are one dimensional on the same top down levels. Refuses a non-positive or
    non-finite density, a non-finite geopotential, and a layer whose geopotential does not
    decrease from level `k` to level `k + 1` (zero thickness, or a column that is not top down),
    naming the level or the layer.

    Written as `rho_k (Phi_k - Phi_{k+1}) expm1(x) / x` with `x = ln(rho_{k+1} / rho_k)`, which is
    the specification's `(rho_{k+1} - rho_k)(Phi_k - Phi_{k+1}) / ln(rho_{k+1} / rho_k)` exactly,
    since `rho_{k+1} - rho_k = rho_k expm1(x)`, and does not lose digits to the difference of two
    nearly equal densities when `x` is small. Where `|x| < LIMIT_LOG_RATIO` the limit
    `rho_k (Phi_k - Phi_{k+1})` is returned.
    """
    rho = np.asarray(rho, dtype="float64")
    Phi = np.asarray(Phi, dtype="float64")
    if rho.ndim != 1 or Phi.ndim != 1 or rho.shape != Phi.shape:
        raise ValueError(
            f"layer_mass needs rho and Phi one dimensional on the same levels; got shapes "
            f"{rho.shape} and {Phi.shape}."
        )
    if rho.size < 2:
        raise ValueError(f"layer_mass needs at least two levels; got {rho.size}.")
    bad = np.flatnonzero(~np.isfinite(rho) | (rho <= 0.0))
    if bad.size:
        raise ValueError(
            f"layer_mass: the density is not positive and finite at levels {bad.tolist()} "
            f"(values {rho[bad].tolist()}); a log-linear layer needs a positive density."
        )
    bad = np.flatnonzero(~np.isfinite(Phi))
    if bad.size:
        raise ValueError(f"layer_mass: non-finite geopotential at levels {bad.tolist()}.")
    thickness = Phi[:-1] - Phi[1:]
    bad = np.flatnonzero(thickness <= 0.0)
    if bad.size:
        k = int(bad[0])
        more = f", and {bad.size - 1} more layers" if bad.size > 1 else ""
        raise ValueError(
            f"layer_mass: the geopotential does not decrease across layer {k} (levels {k} and "
            f"{k + 1}, Phi {float(Phi[k])!r} and {float(Phi[k + 1])!r}){more}; a layer of zero "
            "thickness has no mass, and the column must run top down."
        )
    x = np.log(rho[1:] / rho[:-1])
    small = np.abs(x) < LIMIT_LOG_RATIO
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(small, 1.0, np.expm1(x) / np.where(small, 1.0, x))
    return rho[:-1] * thickness * factor


def pressure_from_top(p_b, I):
    """`p_0 = p_b`, `p_{k+1} = p_k + I_{k+1/2}`: the pressure on the levels from the boundary at
    the top (B7.2). Refuses a boundary pressure that is not positive and finite, and a layer mass
    that is not positive and finite, naming the layer."""
    p_b = float(p_b)
    I = np.asarray(I, dtype="float64")
    if I.ndim != 1:
        raise ValueError(f"pressure_from_top needs one dimensional layer masses; got {I.shape}.")
    if not np.isfinite(p_b) or p_b <= 0.0:
        raise ValueError(f"pressure_from_top: the boundary pressure {p_b!r} Pa is not positive.")
    bad = np.flatnonzero(~np.isfinite(I) | (I <= 0.0))
    if bad.size:
        raise ValueError(
            f"pressure_from_top: the layer mass is not positive and finite in layers "
            f"{bad.tolist()}; integrating from the top, every layer adds pressure."
        )
    return p_b + np.concatenate(([0.0], np.cumsum(I)))


def temperature(p, N, R_bar):
    """`T = p R_bar / (k_B N)`, Eq. B6."""
    p = np.asarray(p, dtype="float64")
    return p * np.asarray(R_bar, dtype="float64") / (
        BOLTZMANN_CONSTANT * np.asarray(N, dtype="float64"))
