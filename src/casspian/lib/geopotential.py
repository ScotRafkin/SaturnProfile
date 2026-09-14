"""Geopotential along a profile's own vertical. Manuscript Eq. B2 restated along the local vertical.

SPEC_03 Step 1. Pure functions, NumPy in and out, no file access, no module level Saturn
quantity. Every angle is in radians.

**Which vertical, and which gravity** (SPEC_03 decision 1). A source profile's altitude `h` is
measured along the local vertical, the field line of the effective gravity, so the geopotential
between two levels is the line integral of the **magnitude** of the effective gravity along
that line, `dPhi = |g_eff| dh`, not the radial component times a radial increment. At the
Lindal latitude the two differ by `1 / cos psi - 1 = 4.6e-3`.

**The staggering** (SPEC_03 Step 1 table, SPEC_00 section 7.3 rule ii). Levels are the
tabulated levels, index `k = 0` at the top, radius decreasing with `k`, as kind N stores them.
`|g_eff|_k`, `psi_k` and `Phi_k` live on levels; the increment `dPhi_{k+1/2}` lives on the layer
between levels `k` and `k + 1`:

    dPhi_{k+1/2} = (|g_eff|_k + |g_eff|_{k+1}) / 2 * (h_{k+1} - h_k)

the trapezoid, which is negative for top-down levels. `|g_eff|` is smooth and varies by under two
percent over a profile, and the trapezoid error of a layer is about `(dh)^2 / (2 r^2)` of it.

**The gauge.** `Phi = 0` at the gauge level `a`. Every other level is the sum of the layer
increments between `a` and itself, taken outward from `a`, so that `Phi_a` is exactly zero and
no level near the gauge carries the round-off of a long cumulative sum. `Phi` increases upward.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from casspian.lib.gravity import g_eff_vector

__all__ = [
    "GeopotentialProfile",
    "effective_gravity_magnitude",
    "geopotential",
    "geopotential_along_profile",
    "layer_increments",
]


@dataclass(frozen=True)
class GeopotentialProfile:
    """Everything `geopotential_along_profile` forms, on levels and layers. Angles in radians."""

    geopotential_m2s2: np.ndarray
    layer_increments_m2s2: np.ndarray
    g_magnitude_ms2: np.ndarray
    g_radial_ms2: np.ndarray
    G_phi_ms2: np.ndarray
    psi_rad: np.ndarray
    gauge_index: int


def effective_gravity_magnitude(u, r, phi_c, Omega, GM, J, degrees, R_norm):
    """`(|g_eff|, g, G_phi, psi)` at every radius `r` of a column at latitude `phi_c`.

    `u` is the zonal wind at `phi_c`, a scalar or an array broadcastable to `r`. From
    `lib.gravity.g_eff_vector`; `|g_eff| = hypot(g, G_phi) = g / cos psi`.
    """
    r = np.asarray(r, dtype="float64")
    shape = r.shape
    u = np.broadcast_to(np.asarray(u, dtype="float64"), shape)
    phi = np.broadcast_to(np.asarray(phi_c, dtype="float64"), shape)
    g, G_phi, magnitude, psi = g_eff_vector(u, r, phi, Omega, GM, J, degrees, R_norm)
    return magnitude, g, G_phi, psi


def layer_increments(g_mag, h):
    """`dPhi_{k+1/2} = (|g_eff|_k + |g_eff|_{k+1}) / 2 * (h_{k+1} - h_k)`, the trapezoid.

    `g_mag` and `h` are one dimensional on the same levels. Refuses fewer than two levels, a shape
    mismatch, a non-finite value, and heights that are not strictly monotonic: a layer of zero
    or reversed thickness has no place on a profile's vertical.
    """
    g_mag = np.asarray(g_mag, dtype="float64")
    h = np.asarray(h, dtype="float64")
    if g_mag.ndim != 1 or h.ndim != 1 or g_mag.shape != h.shape:
        raise ValueError(
            f"layer_increments needs g_mag and h one dimensional on the same levels; got shapes "
            f"{g_mag.shape} and {h.shape}."
        )
    if h.size < 2:
        raise ValueError(f"layer_increments needs at least two levels; got {h.size}.")
    bad = np.flatnonzero(~(np.isfinite(g_mag) & np.isfinite(h)))
    if bad.size:
        raise ValueError(f"layer_increments: non-finite g_mag or h at levels {bad.tolist()}.")
    dh = np.diff(h)
    if not (np.all(dh > 0.0) or np.all(dh < 0.0)):
        layers = np.flatnonzero(np.sign(dh) != np.sign(dh[0])).tolist() if dh[0] != 0.0 else [0]
        raise ValueError(
            f"layer_increments: h is not strictly monotonic; the layers below levels {layers} "
            "have zero thickness or reverse direction."
        )
    return 0.5 * (g_mag[:-1] + g_mag[1:]) * dh


def geopotential(dPhi, gauge_index):
    """`Phi_k` on the levels, zero at `gauge_index`, from the layer increments `dPhi`.

    For `k > a`, `Phi_k = sum_{j=a}^{k-1} dPhi_j`; for `k < a`, `Phi_k = -sum_{j=k}^{a-1} dPhi_j`.
    Both sums run outward from the gauge. Refuses a gauge index that is not an integer, is
    negative, or is not a level of the `len(dPhi) + 1` levels.
    """
    dPhi = np.asarray(dPhi, dtype="float64")
    if dPhi.ndim != 1:
        raise ValueError(f"geopotential needs one dimensional layer increments; got {dPhi.shape}.")
    levels = dPhi.size + 1
    if isinstance(gauge_index, bool) or not isinstance(gauge_index, Integral):
        raise ValueError(
            f"geopotential: gauge_index {gauge_index!r} is not an integer level index."
        )
    a = int(gauge_index)
    if not 0 <= a < levels:
        raise ValueError(
            f"geopotential: gauge_index {a} is off the grid of {levels} levels (0 to "
            f"{levels - 1}); the gauge must be a level."
        )
    Phi = np.empty(levels, dtype="float64")
    Phi[a] = 0.0
    Phi[a + 1:] = np.cumsum(dPhi[a:])
    Phi[:a] = -np.cumsum(dPhi[:a][::-1])[::-1]
    return Phi


def geopotential_along_profile(u, r, h, phi_c, gauge_index, Omega, GM, J, degrees, R_norm):
    """The three steps in sequence on one column: gravity on the levels, the layer
    increments, the geopotential with its gauge. Returns a `GeopotentialProfile`."""
    g_mag, g, G_phi, psi = effective_gravity_magnitude(u, r, phi_c, Omega, GM, J, degrees,
                                                       R_norm)
    dPhi = layer_increments(g_mag, h)
    Phi = geopotential(dPhi, gauge_index)
    return GeopotentialProfile(
        geopotential_m2s2=Phi,
        layer_increments_m2s2=dPhi,
        g_magnitude_ms2=g_mag,
        g_radial_ms2=g,
        G_phi_ms2=G_phi,
        psi_rad=psi,
        gauge_index=int(gauge_index),
    )
