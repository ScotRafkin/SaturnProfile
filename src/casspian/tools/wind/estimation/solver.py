"""The solver, written once. SPEC_19 section 0: it does not change as observation types are added.

It sees the data `y`, the model's equivalent `h` at the current state, the operator `H`, the error
covariance `R` and a group label per row. The linearized `J` of the design note section 2 is solved in
observation space, `x = x_b + B H^T S^-1 (d + H (x_k - x_b))`, `S = H B H^T + R`, `d = y - h(x_k)`, which
gives the posterior, the consistency check and leave-one-out from the same factors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq


@dataclass(frozen=True)
class Blocks:
    """All observations at one linearization."""

    y: np.ndarray
    h: np.ndarray
    H: np.ndarray
    R: np.ndarray
    groups: np.ndarray


def concatenate(blocks) -> Blocks:
    """The observation types' terms as one system. Types are independent of each other."""
    sizes = [b.y.size for b in blocks]
    R = np.zeros((sum(sizes), sum(sizes)))
    start = 0
    for b, n in zip(blocks, sizes):
        R[start:start + n, start:start + n] = b.R
        start += n
    return Blocks(np.concatenate([b.y for b in blocks]), np.concatenate([b.h for b in blocks]),
                  np.vstack([b.H for b in blocks]), R, np.concatenate([b.groups for b in blocks]))


@dataclass(frozen=True)
class Solution:
    state: np.ndarray
    BHt: np.ndarray          #: B H^T, `(state, m)`
    factor: tuple
    innovation: np.ndarray   #: d + H (x_k - x_b)


def solve(prior, H, R, innovation) -> Solution:
    BHt = prior.apply(H).T
    S = H @ BHt + R
    factor = cho_factor(S, lower=True)
    return Solution(prior.mean + BHt @ cho_solve(factor, innovation), BHt, factor, innovation)


def linearized(prior, state, blocks: Blocks) -> Solution:
    return solve(prior, blocks.H, blocks.R, blocks.y - blocks.h + blocks.H @ (state - prior.mean))


def posterior_variance(prior, solution, matrix) -> np.ndarray:
    """The posterior variance of `matrix x`, `(m, state)`, row by row."""
    prior_part = np.einsum("mi,mi->m", prior.apply(matrix), matrix)
    cross = matrix @ solution.BHt
    return prior_part - np.einsum("mi,im->m", cross, cho_solve(solution.factor, cross.T))


def consistency(blocks: Blocks, solution: Solution, x_b=None, exclude=()) -> dict:
    """Each group's weighted misfit against its expected value.

    With `r = innovation - H (x - x_b)` the posterior residual, `E[r r^T] = R S^-1 R` if the stated errors
    are right, so a group's `r^T R_gg^-1 r` has the expected value `tr(R_gg^-1 (R S^-1 R)_gg)`. Groups in
    `exclude`, constraints rather than observations, are not reported, and `all` is the observations' alone,
    `d^T S_dd^-1 d / m` over their rows."""
    r = blocks.R @ cho_solve(solution.factor, solution.innovation)
    expected_cov = blocks.R @ cho_solve(solution.factor, blocks.R)
    out = {}
    for group in dict.fromkeys(blocks.groups.tolist()):
        if group in exclude:
            continue
        g = blocks.groups == group
        factor = cho_factor(blocks.R[np.ix_(g, g)], lower=True)
        misfit = float(r[g] @ cho_solve(factor, r[g]))
        expected = float(np.trace(cho_solve(factor, expected_cov[np.ix_(g, g)])))
        out[group] = {"count": int(g.sum()), "misfit": misfit, "expected": expected,
                      "ratio": misfit / expected if expected > 0 else float("nan")}
    data = ~np.isin(blocks.groups, list(exclude))
    innovation = solution.innovation[data]
    if data.all():
        total = float(innovation @ cho_solve(solution.factor, innovation))
    else:
        L = np.tril(solution.factor[0])
        S = (L @ L.T)[np.ix_(data, data)]
        total = float(innovation @ np.linalg.solve(S, innovation))
    out["all"] = {"count": int(innovation.size), "misfit": total, "expected": float(innovation.size),
                  "ratio": total / innovation.size}
    return out


def leave_out(prior, blocks: Blocks, innovation, group) -> dict:
    """`group` withheld: the rest's solution and its prediction of the group, on this linearization.

    `error` is the data minus the prediction, `variance` the stated error's plus the prediction's."""
    g = blocks.groups == group
    rest = solve(prior, blocks.H[~g], blocks.R[np.ix_(~g, ~g)], innovation[~g])
    delta = rest.state - prior.mean
    # The withheld rows' errors are correlated with the rest's through the shared anchor term.
    cross = blocks.R[np.ix_(g, ~g)]
    from_errors = cross @ cho_solve(rest.factor, innovation[~g])
    error = innovation[g] - blocks.H[g] @ delta - from_errors
    prediction = blocks.H[g]
    BHt_g = prior.apply(prediction).T
    a = prediction @ BHt_g
    c = blocks.H[~g] @ BHt_g + cross.T
    covariance = a + blocks.R[np.ix_(g, g)] - c.T @ cho_solve(rest.factor, c)
    return {"state": rest.state, "error": error, "variance": np.diag(covariance),
            "weighted": float(error @ np.linalg.solve(covariance, error)), "count": int(g.sum())}


def systematic_term(prior, blocks: Blocks, innovation, group, shape) -> float:
    """The size of a term `sigma^2 shape` added to `group`'s errors that brings its weighted misfit to
    its expected value, re-solving as it grows. Zero if the group is already at or below it."""
    g = np.flatnonzero(blocks.groups == group)

    def excess(sigma):
        R = blocks.R.copy()
        R[np.ix_(g, g)] += sigma ** 2 * shape
        solved = solve(prior, blocks.H, R, innovation)
        check = consistency(Blocks(blocks.y, blocks.h, blocks.H, R, blocks.groups), solved)[group]
        return check["misfit"] - check["expected"]

    if excess(0.0) <= 0.0:
        return 0.0
    high = max(float(np.sqrt(np.mean(np.diag(blocks.R)[g]))), 1e-12)
    for _ in range(40):
        if excess(high) < 0.0:
            return float(brentq(excess, 0.0, high, xtol=high * 1e-3))
        high *= 2.0
    return float("nan")


def robust(prior, state, blocks: Blocks, threshold, iterations=30, tolerance=1e-4, exempt=()):
    """Huber's rule on the whitened residuals within each group (SPEC_21 section 1 item 8, REVIEW_21_step0
    ruling 5). With `R_gg = L L^T`, the posterior residual `r = R S^-1 innovation` is whitened, `z = L^-1 r`; a row
    with `|z|` above `threshold` gets the weight `threshold / |z|`, and the group's block becomes
    `L diag(1 / weight) L^T`, which only adds a positive part, so `R` stays positive definite. Repeated until the
    weights settle. Returns the solution, the blocks with the weighted `R`, and the weights (1 where untouched).
    Groups in `exempt`, constraints rather than observations, keep their weight."""
    groups = [g for g in dict.fromkeys(blocks.groups.tolist()) if g not in exempt]
    factors = {g: np.linalg.cholesky(blocks.R[np.ix_(blocks.groups == g, blocks.groups == g)]) for g in groups}
    weights = np.ones(blocks.y.size)
    for _ in range(int(iterations)):
        R = blocks.R.copy()
        for g in groups:
            rows = blocks.groups == g
            L = factors[g]
            R[np.ix_(rows, rows)] = L @ np.diag(1.0 / weights[rows]) @ L.T
        weighted = Blocks(blocks.y, blocks.h, blocks.H, R, blocks.groups)
        solution = linearized(prior, state, weighted)
        r = R @ cho_solve(solution.factor, solution.innovation)
        new = np.ones_like(weights)
        for g in groups:
            rows = blocks.groups == g
            z = np.linalg.solve(factors[g], r[rows])
            new[rows] = np.minimum(1.0, float(threshold) / np.maximum(np.abs(z), 1e-300))
        if np.max(np.abs(new - weights)) < tolerance:
            weights = new
            break
        weights = new
    return solution, weighted, weights


def cost_parts(prior, state, blocks: Blocks, solution: Solution, shared=None) -> dict:
    """The cost `J` of one solve in its parts, before (at `state`, the linearization's) and after (at the solve's
    linear prediction). `J = (x - x_b)^T B^-1 (x - x_b) + r^T R^-1 r`. With `v = R^-1 r`, `r^T R^-1 r = v^T R v`
    splits over any sum `R = sum_k R_k`: the shared anchor term `shared` on the first rows, then each group's own
    block less its share of it. After the solve `v` is `S^-1 innovation` and the prior's part `v^T H B H^T v`.
    `blocks` is the one the solve used (Huber's weighting in `R`)."""
    m = 0 if shared is None else shared.shape[0]
    S_full = np.zeros_like(blocks.R)
    if m:
        S_full[:m, :m] = shared

    def split(v):
        parts = {}
        if m:
            parts["shared anchor"] = float(v @ S_full @ v)
        for group in dict.fromkeys(blocks.groups.tolist()):
            g = blocks.groups == group
            parts[group] = float(v[g] @ (blocks.R[np.ix_(g, g)] - S_full[np.ix_(g, g)]) @ v[g])
        return parts

    dx = state - prior.mean
    B = prior.covariance()
    before = {"prior": float(dx @ np.linalg.lstsq(B, dx, rcond=1e-12)[0])}
    before.update(split(cho_solve(cho_factor(blocks.R, lower=True), blocks.y - blocks.h)))
    v = cho_solve(solution.factor, solution.innovation)
    after = {"prior": float(v @ (blocks.H @ (solution.BHt @ v)))}
    after.update(split(v))
    return {"before": before, "after": after}
