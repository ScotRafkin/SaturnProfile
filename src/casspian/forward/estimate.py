"""The anchor constant from M anchors at the gauge latitude. SPEC_04 Step 4 deliverable 2.

Eq. A29 gives each anchor's estimate `C_i(Phi)` of the anchor constant: its own `ln N`
transferred along its isobars to the gauge latitude. Eq. A30 combines them by inverse variance,
A33 is the difference of any two, and A34 places the gauge latitude itself.

**Where the uncertainty comes from (decision J).** Each anchor arrives carrying its own two
columns, the measurement term and the season term, and this module reads them from the anchor
and from no file. The transfer term is the placeholder `P_i = (sigma_K (phi_i - phi_r))^2` of a
declared per-run `sigma_K`, which the combination specification replaces with the covariance of
Eq. A35.

**What M = 1 is.** The identity: `phi_r` is the one anchor's latitude, the union of arrival
levels is its own, `C` is its `ln N` with no interpolation, its weight is the only one, and no
difference `D` exists. Nothing in the code takes a different path to get there.

**The gauge latitude is fixed before the mesh (decision K).** A34's weights carry `P_i`, which
depends on `phi_r` itself and on `Phi`; this module takes the centroid weighted by
`1 / mean(sigma_i^2)` over each anchor's levels instead, which is one number per anchor and
keeps `phi_r` out of the fixed point and out of the vertical. The departure is recorded here and
in the product.

Pure functions and NumPy arrays, no file access. Every angle is in radians.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

__all__ = ["Arrival", "Estimate", "GAUGE_LATITUDE_RULE", "sigma_ln_N", "gauge_latitude",
           "estimate"]

#: The one rule implemented, and what the namelist's `gauge_latitude_rule` must name.
GAUGE_LATITUDE_RULE = "weighted_centroid"


def sigma_ln_N(anchor) -> np.ndarray:
    """`sigma_i(Phi)` on an anchor's own levels: the measurement and season terms in quadrature.

    Decision J: the two columns arrive on the anchor object and are read from there. The season
    term is zero until the propagator specification fills it, and the anchor's `season_term`
    says so.
    """
    measurement = np.asarray(anchor.sigma_ln_N_measurement, dtype="float64")
    season = np.asarray(anchor.sigma_ln_N_season, dtype="float64")
    return np.sqrt(measurement ** 2 + season ** 2)


def gauge_latitude(anchors) -> float:
    """`phi_r` in radians, the weighted centroid of the construction anchors (decision K).

    The weight of an anchor is `1 / mean(sigma_i^2)` over its own levels. At M = 1 the result is
    the one anchor's latitude exactly, with no arithmetic that could move it.
    """
    construction = [a for a in anchors if a.is_construction]
    if not construction:
        raise ValueError(
            "no anchor has weight 1: the gauge latitude is the centroid of the construction "
            "anchors and there is nothing to place it with (SPEC_04 Step 4 deliverable 2)"
        )
    latitude = np.radians(np.asarray([a.latitude_planetocentric_deg for a in construction],
                                     dtype="float64"))
    if len(construction) == 1:
        return float(latitude[0])
    weight = np.asarray([1.0 / float(np.mean(sigma_ln_N(a) ** 2)) for a in construction],
                        dtype="float64")
    return float(np.sum(weight * latitude) / np.sum(weight))


@dataclass(frozen=True)
class Arrival:
    """One anchor as it reaches the gauge latitude, on its own arrival levels.

    `geopotential_m2s2` is where its isobars arrive, `C_i` its transferred `ln N` there (A29),
    `sigma_ln_N` its own uncertainty on those levels, and `label_pressure_Pa` the labels the
    isobars carry. A validation anchor (`weight = 0`) arrives the same way and is reported the
    same way; it enters neither `phi_r` nor `C`.
    """

    slug: str
    latitude_rad: float
    weight: float
    geopotential_m2s2: np.ndarray
    C_i: np.ndarray
    sigma_ln_N: np.ndarray
    label_pressure_Pa: np.ndarray

    @property
    def is_construction(self) -> bool:
        return self.weight == 1.0


@dataclass(frozen=True)
class Estimate:
    """The anchor constant on the union of the anchors' arrival levels.

    Every array on the union is `(level,)` and every per-anchor array `(anchor, level)`, in the
    order the arrivals were given. `present` says which levels each anchor covers; outside its
    own span an anchor is absent, not extrapolated, and its entries are NaN.
    """

    gauge_latitude_rad: float
    slug: tuple
    weight: np.ndarray
    geopotential_m2s2: np.ndarray
    present: np.ndarray
    C_i: np.ndarray
    sigma_ln_N: np.ndarray
    transfer_variance: np.ndarray
    estimate_weight: np.ndarray
    C: np.ndarray
    variance: np.ndarray
    chi_square_reduced: np.ndarray
    variance_inflated: np.ndarray
    label_pressure_Pa: np.ndarray
    label_spread_ln_p: np.ndarray
    differences: MappingProxyType
    record: MappingProxyType


def _onto(union, own, values):
    """One anchor's column on the union, linear in `Phi`, absent outside its own span."""
    rising = np.argsort(own)
    x, y = own[rising], values[rising]
    out = np.interp(union, x, y)
    out[(union < x[0]) | (union > x[-1])] = np.nan
    return out


def estimate(arrivals, kernel_uncertainty_per_rad, gauge_latitude_rad) -> Estimate:
    """`C(Phi)` and its companions at the gauge latitude. SPEC_04 Step 4 deliverable 2.

    The union of the arrival levels; each anchor's `C_i`, `sigma_i` and label interpolated onto
    it within its own span (decision C: the one interpolation of `N` in the chain, and none of
    it at M = 1, where the union is the anchor's own levels); weights
    `w_i = 1 / (sigma_i^2 + P_i)` with `P_i = (sigma_K (phi_i - phi_r))^2`; `C` by A30 over the
    construction anchors present at each level, its variance `1 / sum w_i` inflated by the
    reduced chi-square where that exceeds one; `D_ij = C_i - C_j` (A33) for every pair,
    validation anchors included, on the levels both cover.

    The isobar label of a union level is the weighted mean of the labels the anchors that reach
    it carry, in `ln p` and with the same weights, and the spread between them is recorded per
    level; at M = 1 it is the anchor's own label exactly and the spread is zero.

    The union is the construction anchors' arrival levels. A validation anchor is placed, traced
    and reported like any other, but it enters neither the levels the estimate is made on nor the
    mean, which is what leaves `C` at the values it would have without it.
    """
    arrivals = list(arrivals)
    if not arrivals:
        raise ValueError("the estimate needs at least one anchor")
    sigma_K = float(kernel_uncertainty_per_rad)
    building_anchors = [a for a in arrivals if a.is_construction]
    if not building_anchors:
        raise ValueError(
            "no anchor has weight 1: there is nothing to make the estimate from "
            "(SPEC_04 Step 4 deliverable 2)"
        )
    union = np.unique(np.concatenate([np.asarray(a.geopotential_m2s2, dtype="float64")
                                      for a in building_anchors]))

    C_i = np.stack([_onto(union, np.asarray(a.geopotential_m2s2, dtype="float64"),
                          np.asarray(a.C_i, dtype="float64")) for a in arrivals])
    sigma = np.stack([_onto(union, np.asarray(a.geopotential_m2s2, dtype="float64"),
                            np.asarray(a.sigma_ln_N, dtype="float64")) for a in arrivals])
    ln_label = np.stack([_onto(union, np.asarray(a.geopotential_m2s2, dtype="float64"),
                               np.log(np.asarray(a.label_pressure_Pa, dtype="float64")))
                         for a in arrivals])
    present = np.isfinite(C_i)
    weight = np.asarray([float(a.weight) for a in arrivals], dtype="float64")
    latitude = np.asarray([float(a.latitude_rad) for a in arrivals], dtype="float64")
    P = (sigma_K * (latitude - float(gauge_latitude_rad))) ** 2

    # An absent level carries NaN, and `np.where` evaluates both branches, so the arithmetic on
    # the absent entries is done and thrown away rather than warned about.
    with np.errstate(invalid="ignore"):
        w = np.where(present, 1.0 / (sigma ** 2 + P[:, None]), 0.0)
    building = present & (weight[:, None] == 1.0)
    w_build = np.where(building, w, 0.0)
    total = np.sum(w_build, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        C = np.sum(np.where(building, w * C_i, 0.0), axis=0) / total
        variance = 1.0 / total
        count = np.sum(building, axis=0)
        chi = np.sum(np.where(building, w * (C_i - C) ** 2, 0.0), axis=0) / np.maximum(count - 1,
                                                                                       1)
        chi = np.where(count >= 2, chi, np.nan)
        label = np.exp(np.sum(np.where(building, w * ln_label, 0.0), axis=0) / total)
    C[total == 0.0] = np.nan
    variance[total == 0.0] = np.nan
    label[total == 0.0] = np.nan
    inflated = np.where(np.isfinite(chi) & (chi > 1.0), variance * chi, variance)
    spread = np.where(np.sum(building, axis=0) >= 2,
                      np.nanmax(np.where(building, ln_label, np.nan), axis=0)
                      - np.nanmin(np.where(building, ln_label, np.nan), axis=0), 0.0)

    differences = {}
    for i in range(len(arrivals)):
        for j in range(i + 1, len(arrivals)):
            both = present[i] & present[j]
            D = np.where(both, C_i[i] - C_i[j], np.nan)
            differences[(arrivals[i].slug, arrivals[j].slug)] = D

    return Estimate(
        gauge_latitude_rad=float(gauge_latitude_rad),
        slug=tuple(a.slug for a in arrivals),
        weight=weight,
        geopotential_m2s2=union,
        present=present,
        C_i=C_i,
        sigma_ln_N=sigma,
        transfer_variance=P,
        estimate_weight=w,
        C=C,
        variance=variance,
        chi_square_reduced=chi,
        variance_inflated=inflated,
        label_pressure_Pa=label,
        label_spread_ln_p=spread,
        differences=MappingProxyType(differences),
        record=MappingProxyType({
            "gauge_latitude_rule": GAUGE_LATITUDE_RULE,
            "gauge_latitude_deg": float(np.degrees(gauge_latitude_rad)),
            "kernel_uncertainty_per_rad": sigma_K,
            "anchors": len(arrivals),
            "construction_anchors": int(np.sum(weight == 1.0)),
            "union_levels": int(union.size),
            "interpolated": bool(len(arrivals) > 1),
            "gauge_latitude_departure_from_A34": (
                "the centroid is weighted by 1 / mean(sigma_i^2) over each anchor's levels, "
                "not by A34's weights, which carry P_i and depend on Phi and on phi_r itself "
                "(SPEC_04 decision K)"),
        }),
    )
