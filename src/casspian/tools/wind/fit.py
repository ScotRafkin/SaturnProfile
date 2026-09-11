"""Constrained penalized spline fit of point data against latitude, with quality checks.

SPEC_01 v0.11 Step 7, item 2. A cubic B-spline with a second derivative penalty, fitted by
least squares to scattered points subject to equality constraints at the endpoints, with the
smoothing parameter set by matching the residual scatter to the scatter the data itself shows.

Manuscript equations implemented: none. This is a fitting method.

The quality rules of SPEC_01 v0.11 are written to be general, to be applied to any latitude fit
of point data in this project. They live here rather than in `lib` because Step 7 is the only
caller today and SPEC_00 section 3.1 does not list a fitting module; when a second latitude fit
appears this module is what moves to `lib`.

Why a penalty and not an interpolant. The 323 digitized points carry a reading error of about
2.5 m/s and, in the merged clusters, up to 2 degrees of latitude error. An interpolant would
chase that noise. A penalized fit with the smoothing chosen so that the residual scatter equals
the within bin scatter follows the jets and leaves the noise in the residuals, which is then
what the uncertainty is built from.

Why equality constraints and not weighted pseudo-points at the poles. A zonal wind is zero at a
pole as a matter of geometry, not of evidence, and SPEC_00 section 6.6 makes a nonzero polar
wind a defect in any kind W file. A heavy weight would make it nearly zero; a constraint makes
it exactly zero, which is what the load time check demands.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import BSpline, PchipInterpolator
from scipy.optimize import brentq

__all__ = ["FitResult", "fit_constrained_spline", "pooled_within_bin_scatter"]

_DEGREE = 3
_GAUSS_POINTS = 4


@dataclass
class FitResult:
    """A fitted wind curve and everything the report has to quote about it."""

    evaluate: object
    method: str
    smoothing: float
    smoothing_scatter_match: float
    smoothing_steps: int
    forced_by: str
    rms_residual: float
    pooled_scatter: float
    hard_checks: dict = field(default_factory=dict)
    soft_checks: dict = field(default_factory=dict)
    warnings: tuple = ()


def pooled_within_bin_scatter(values, bin_index, counts) -> float:
    """Pooled within bin sample standard deviation of the points about their bin means.

    This is the scatter the data shows at a scale the fit is not expected to resolve, and it is
    the target the smoothing parameter is matched to. Bins holding one point contribute no
    degrees of freedom and are skipped.
    """
    total, dof = 0.0, 0
    for i in np.flatnonzero(counts >= 2):
        members = values[bin_index == i]
        total += float(np.sum((members - members.mean()) ** 2))
        dof += members.size - 1
    if dof == 0:
        raise ValueError("no bin holds two or more points, so no scatter can be pooled")
    return float(np.sqrt(total / dof))


def _knot_vector(lo: float, hi: float, spacing: float) -> np.ndarray:
    interior = np.arange(lo, hi + 0.5 * spacing, spacing)
    return np.concatenate([np.full(_DEGREE, lo), interior, np.full(_DEGREE, hi)])


def _basis(knots: np.ndarray, x: np.ndarray, derivative: int = 0) -> np.ndarray:
    """Design matrix of the B-spline basis, or of its `derivative`-th derivative."""
    n_basis = knots.size - _DEGREE - 1
    out = np.empty((np.size(x), n_basis), dtype="float64")
    identity = np.eye(n_basis)
    for i in range(n_basis):
        spline = BSpline(knots, identity[i], _DEGREE, extrapolate=False)
        values = spline(np.atleast_1d(x), nu=derivative)
        out[:, i] = np.nan_to_num(values, nan=0.0)
    return out


def _penalty(knots: np.ndarray) -> np.ndarray:
    """Exact `integral of (d2f/dx2)^2`, assembled by Gauss-Legendre on each knot span.

    The second difference penalty of a P-spline is the usual shortcut; the specification says
    second derivative, so the integral is assembled rather than approximated.
    """
    nodes, weights = np.polynomial.legendre.leggauss(_GAUSS_POINTS)
    spans = np.unique(knots)
    n_basis = knots.size - _DEGREE - 1
    total = np.zeros((n_basis, n_basis), dtype="float64")
    for lo, hi in zip(spans[:-1], spans[1:]):
        half = 0.5 * (hi - lo)
        x = 0.5 * (lo + hi) + half * nodes
        second = _basis(knots, x, derivative=2)
        total += (second * (weights * half)[:, None]).T @ second
    return total


def _solve(design, penalty, values, constraint, smoothing):
    """Penalized least squares subject to `constraint @ c = 0`, by the KKT system."""
    n = design.shape[1]
    normal = design.T @ design + smoothing * penalty
    rhs = design.T @ values
    m = constraint.shape[0]
    kkt = np.block([[normal, constraint.T], [constraint, np.zeros((m, m))]])
    full = np.concatenate([rhs, np.zeros(m)])
    return np.linalg.solve(kkt, full)[:n]


def _extrema(x, y):
    """Indices of interior local extrema, with +1 for a maximum and -1 for a minimum."""
    rising = np.diff(y)
    sign = np.sign(rising)
    nonzero = sign != 0
    idx, kinds = [], []
    last = None
    for i in range(sign.size):
        if not nonzero[i]:
            continue
        if last is not None and sign[i] != last[1]:
            idx.append(last[0] + 1)
            kinds.append(1 if last[1] > 0 else -1)
        last = (i, sign[i])
    return np.array(idx, dtype=int), np.array(kinds, dtype=int)


def _hard_checks(evaluate, points_x, points_y, centers, means, counts, min_count,
                 width, scatter, dense):
    """SPEC_01 v0.11 checks (a) to (d). Returns a dict of name to (passed, measured)."""
    fitted = evaluate(dense)
    populated = counts > 0
    first, last = float(centers[populated][0]), float(centers[populated][-1])
    inside = (dense >= first) & (dense <= last)

    # (a) Range
    hi, lo = float(points_y.max()), float(points_y.min())
    over = float(np.max(fitted[inside]) - hi)
    under = float(lo - np.min(fitted[inside]))
    local_excess = 0.0
    for x in dense[inside]:
        near = np.abs(points_x - x) <= width
        if not near.any():
            continue
        value = float(evaluate(np.array([x]))[0])
        local_excess = max(
            local_excess,
            value - (float(points_y[near].max()) + scatter),
            (float(points_y[near].min()) - scatter) - value,
        )
    range_ok = over <= scatter and under <= scatter and local_excess <= 0.0
    checks = {
        "a_range": (
            range_ok,
            f"overshoot above the point maximum {max(over, 0.0):.3f} m/s and below the "
            f"minimum {max(under, 0.0):.3f} m/s against a tolerance of s = {scatter:.3f}; "
            f"worst local excess beyond the points within +-{width} deg widened by s: "
            f"{max(local_excess, 0.0):.3f} m/s",
        )
    }

    # (b) Gaps: no extremum in an interior interval holding no points within w
    gap_ok, gap_detail = True, "no interior gap wider than the bin width"
    empty = populated.copy()
    empty[:] = False
    for i in np.flatnonzero(~(counts > 0)):
        if centers[i] < first or centers[i] > last:
            continue
        empty[i] = True
    if empty.any():
        offenders = []
        for i in np.flatnonzero(empty):
            window = (dense >= centers[i] - 0.5 * width) & (dense <= centers[i] + 0.5 * width)
            if window.sum() < 3:
                continue
            idx, _ = _extrema(dense[window], fitted[window])
            if idx.size:
                offenders.append(float(centers[i]))
        gap_ok = not offenders
        gap_detail = (f"{int(empty.sum())} interior empty bins; extrema found in "
                      f"{len(offenders)}: {offenders}")
    checks["b_gaps"] = (gap_ok, gap_detail)

    # (c) Polar caps: |u| decreases monotonically to zero, no sign change
    cap_ok, cap_parts = True, []
    for side, mask in (("north", dense >= last), ("south", dense <= first)):
        x = dense[mask]
        if x.size < 3:
            continue
        y = evaluate(x)
        if side == "south":
            x, y = x[::-1], y[::-1]
        magnitude = np.abs(y)
        rising = float(np.max(np.diff(magnitude)))
        signs = np.sign(y[np.abs(y) > 1e-12])
        changed = bool(signs.size and not np.all(signs == signs[0]))
        ok = rising <= 1e-9 and not changed
        cap_ok = cap_ok and ok
        cap_parts.append(
            f"{side}: largest rise in |u| {max(rising, 0.0):.3e} m/s, sign change {changed}, "
            f"endpoint {float(y[-1]):.3e} m/s"
        )
    checks["c_polar_caps"] = (cap_ok, "; ".join(cap_parts))

    # (d) Extrema correspondence
    span = inside
    idx, kinds = _extrema(dense[span], fitted[span])
    good = (counts >= min_count) & np.isfinite(means)
    bin_idx, bin_kinds = _extrema(centers[good], means[good])
    unmatched = []
    for position, kind in zip(dense[span][idx], kinds):
        near = np.abs(centers[good][bin_idx] - position) <= 2.0 * width
        if not np.any(near & (bin_kinds == kind)):
            unmatched.append((round(float(position), 2), int(kind)))
    checks["d_extrema"] = (
        not unmatched,
        f"{idx.size} extrema in the fit over the observed span, {bin_idx.size} in the bin "
        f"means; unmatched within +-{2 * width} deg: {unmatched if unmatched else 'none'}",
    )
    return checks


def _soft_checks(evaluate, points_x, points_y, bin_index, centers, counts, scatter,
                 refit, min_count):
    """SPEC_01 v0.11 checks (e) and (f). Check (g) is the report figure."""
    residual = points_y - evaluate(points_x)
    populated = np.flatnonzero(counts > 0)
    within, failures = 0, []
    for i in populated:
        members = residual[bin_index == i]
        n = members.size
        bias = float(members.mean())
        if abs(bias) <= 2.0 * scatter / np.sqrt(n):
            within += 1
        else:
            failures.append((round(float(centers[i]), 1), round(bias, 2)))
    fraction = within / populated.size
    checks = {
        "e_bias": (
            fraction >= 0.95,
            f"{within} of {populated.size} populated bins ({100 * fraction:.1f} percent) have "
            f"a mean residual within 2 s / sqrt(n); failures "
            f"{failures if failures else 'none'}",
        )
    }

    keep = counts >= min_count
    predicted = []
    for i in np.flatnonzero(keep):
        hold = bin_index == i
        if not hold.any():
            continue
        curve = refit(~hold)
        if curve is None:
            continue
        predicted.append(points_y[hold] - curve(points_x[hold]))
    if predicted:
        joined = np.concatenate(predicted)
        rmse = float(np.sqrt(np.mean(joined**2)))
    else:
        rmse = float("nan")
    checks["f_generalization"] = (
        rmse <= 1.5 * scatter,
        f"leave one bin out prediction RMSE {rmse:.3f} m/s against 1.5 s = "
        f"{1.5 * scatter:.3f} m/s",
    )
    return checks


def fit_constrained_spline(
    points_x,
    points_y,
    bin_index,
    centers,
    means,
    counts,
    *,
    lo,
    hi,
    knot_spacing,
    bin_width,
    min_count,
    escalation_factor,
    max_steps,
    dense_step=0.05,
):
    """Fit the wind curve and run the quality checks, escalating the smoothing if needed."""
    points_x = np.asarray(points_x, dtype="float64")
    points_y = np.asarray(points_y, dtype="float64")
    knots = _knot_vector(lo, hi, knot_spacing)
    design = _basis(knots, points_x)
    penalty = _penalty(knots)
    constraint = _basis(knots, np.array([lo, hi]))
    scatter = pooled_within_bin_scatter(points_y, bin_index, counts)
    dense = np.arange(lo, hi + 0.5 * dense_step, dense_step)

    def curve_for(smoothing, mask=None):
        d = design if mask is None else design[mask]
        y = points_y if mask is None else points_y[mask]
        coefficients = _solve(d, penalty, y, constraint, smoothing)
        spline = BSpline(knots, coefficients, _DEGREE, extrapolate=False)

        def evaluate(x):
            x = np.clip(np.asarray(x, dtype="float64"), lo, hi)
            return np.nan_to_num(spline(x), nan=0.0)

        return evaluate

    def rms_for(smoothing):
        evaluate = curve_for(smoothing)
        return float(np.sqrt(np.mean((points_y - evaluate(points_x)) ** 2)))

    # The residual RMS rises with the smoothing, so the scatter match is a bracketed root.
    low, high = 1e-8, 1e12
    if rms_for(low) > scatter:
        matched = low
    elif rms_for(high) < scatter:
        matched = high
    else:
        matched = float(brentq(lambda t: rms_for(np.exp(t)) - scatter,
                               np.log(low), np.log(high), xtol=1e-6))
        matched = float(np.exp(matched))

    smoothing = matched
    for step in range(int(max_steps) + 1):
        evaluate = curve_for(smoothing)
        hard = _hard_checks(evaluate, points_x, points_y, centers, means, counts, min_count,
                            bin_width, scatter, dense)
        failed = [name for name, (ok, _) in hard.items() if not ok]
        if not failed:
            soft = _soft_checks(
                evaluate, points_x, points_y, bin_index, centers, counts, scatter,
                lambda mask: curve_for(smoothing, mask), min_count,
            )
            warnings = tuple(f"{name}: {detail}" for name, (ok, detail) in soft.items()
                             if not ok)
            return FitResult(
                evaluate=evaluate,
                method="constrained_penalized_spline",
                smoothing=smoothing,
                smoothing_scatter_match=matched,
                smoothing_steps=step,
                forced_by="scatter match" if step == 0 else f"raised to clear {failed_last}",
                rms_residual=rms_for(smoothing),
                pooled_scatter=scatter,
                hard_checks={k: v for k, v in hard.items()},
                soft_checks=soft,
                warnings=warnings,
            )
        failed_last = ", ".join(failed)
        smoothing *= float(escalation_factor)

    # Declared fallback: a shape preserving interpolant through the populated bin means with
    # the polar zeros. PCHIP cannot overshoot its data, so checks (a) to (c) hold by
    # construction.
    keep = counts >= min_count
    x = np.concatenate([[lo], centers[keep], [hi]])
    y = np.concatenate([[0.0], means[keep], [0.0]])
    order = np.argsort(x)
    interpolant = PchipInterpolator(x[order], y[order], extrapolate=False)

    def evaluate(xx):
        xx = np.clip(np.asarray(xx, dtype="float64"), lo, hi)
        return np.nan_to_num(interpolant(xx), nan=0.0)

    return FitResult(
        evaluate=evaluate,
        method="pchip_fallback",
        smoothing=float("nan"),
        smoothing_scatter_match=matched,
        smoothing_steps=int(max_steps),
        forced_by=f"no smoothing up to {max_steps} steps cleared {failed_last}",
        rms_residual=float(np.sqrt(np.mean((points_y - evaluate(points_x)) ** 2))),
        pooled_scatter=scatter,
        hard_checks={},
        soft_checks={},
        warnings=("fell back to the shape preserving interpolant",),
    )
