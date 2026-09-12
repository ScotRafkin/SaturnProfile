"""Number density, absolute radius and refractivity. Manuscript Appendix B3.1, B1 and B3.3.

SPEC_02 Step 3. Pure functions, NumPy in and out, no file access, no module level Saturn
quantity. The one constant of nature used, the Boltzmann constant, comes from `lib.constants`.

**Uncertainty is first order propagation of declared uncertainties, and nothing else**
(SPEC_02 v0.4). The functions here supply the partial derivatives of the algebra and one rule
for combining terms, `quadrature`: terms are added in quadrature with no covariance, and a term
that is NaN is **left out**, not counted as zero. Which terms entered and which were left out is
returned beside the total, so a product can say so.

**Mole fractions obey a closure, and the derivatives respect it.** For a composition declared
as a share `s` of a remainder, with some species assigned first (Lindal: ammonia first, then the
remainder split `s` to H2 and `1 - s` to He),

    x_share   = s (1 - x_first)
    x_partner = (1 - s)(1 - x_first)

so a mean property `q_bar = sum_i x_i q_i` has

    d q_bar / d s       = (1 - x_first)(q_share - q_partner)
    d q_bar / d x_k     = q_k - s q_share - (1 - s) q_partner      for a species k assigned first

and never the derivative obtained by moving one mole fraction with the others held fixed, which
would break the sum rule.
"""

from __future__ import annotations

import numpy as np

from casspian.lib.constants import BOLTZMANN_CONSTANT

__all__ = [
    "absolute_radius",
    "mean_over_species",
    "number_density",
    "quadrature",
    "refractivity",
    "share_of_remainder_partials",
    "temperature_from_refractivity",
    "to_one_sigma",
]

#: SPEC_00 section 5 v0.13: how each declared uncertainty kind becomes one standard deviation
#: before quadrature. `1sigma` needs no conversion. `stated` is taken as `1sigma`, and that is
#: an assumption, so it is listed. `range` is the half width of a uniform distribution.
ONE_SIGMA_RULES = {
    "1sigma": (1.0, None),
    "stated": (1.0, "stated taken as 1sigma"),
    "range": (1.0 / np.sqrt(3.0), "range taken as a uniform half width, divided by sqrt(3)"),
}


def to_one_sigma(value, kind):
    """Convert a declared uncertainty to one standard deviation under SPEC_00 section 5 v0.13.

    Returns `(value, conversion)`, where `conversion` is the text of the rule applied, or None
    when the kind is already `1sigma`. An unknown kind is refused: guessing how to read a
    declared uncertainty is exactly what the rule exists to prevent.
    """
    if kind not in ONE_SIGMA_RULES:
        raise ValueError(
            f"uncertainty_kind {kind!r} has no conversion to one standard deviation; SPEC_00 "
            f"section 5 v0.13 defines {sorted(ONE_SIGMA_RULES)}."
        )
    factor, conversion = ONE_SIGMA_RULES[kind]
    return np.asarray(value, dtype="float64") * factor, conversion


def number_density(p_Pa, T_K):
    """`n = p / (k_B T)`. Appendix B3.1."""
    p = np.asarray(p_Pa, dtype="float64")
    T = np.asarray(T_K, dtype="float64")
    return p / (BOLTZMANN_CONSTANT * T)


def absolute_radius(h_m, h_ref_m, r0_m):
    """`r = r0 + (h - h_ref)`. Appendix B1.

    Written with the difference taken first, so that at the anchor level, where `h` is `h_ref`
    exactly, the result is `r0` exactly.
    """
    h = np.asarray(h_m, dtype="float64")
    return float(r0_m) + (h - float(h_ref_m))


def mean_over_species(x, per_molecule):
    """`q_bar = sum_i x_i q_i` at every point. `x` has the species axis first."""
    x = np.asarray(x, dtype="float64")
    q = np.asarray(per_molecule, dtype="float64")
    return np.tensordot(q, x, axes=(0, 0))


def refractivity(n_m3, mean_refractivity_m3):
    """`N = n R_bar`, unscaled. Appendix B3.3."""
    return np.asarray(n_m3, dtype="float64") * np.asarray(mean_refractivity_m3, dtype="float64")


def temperature_from_refractivity(p_Pa, mean_refractivity_m3, N):
    """`T = p R_bar / (k_B N)`, the inverse of B3.1 and B3.3 together."""
    p = np.asarray(p_Pa, dtype="float64")
    return p * np.asarray(mean_refractivity_m3, dtype="float64") / (
        BOLTZMANN_CONSTANT * np.asarray(N, dtype="float64"))


def share_of_remainder_partials(x, per_molecule, share, partner):
    """Partial derivatives of a mean property under a share of remainder closure.

    `x` has the species axis first; `share` and `partner` are the indices of the two species
    that split the remainder. Every other species is taken as assigned first.

    Returns `(s, d_ds, d_dx_first)`: the share at every point, `d q_bar / d s`, and a dict from
    the index of each species assigned first to `d q_bar / d x_k` with the closure applied.
    """
    x = np.asarray(x, dtype="float64")
    q = np.asarray(per_molecule, dtype="float64")
    remainder = x[share] + x[partner]
    s = x[share] / remainder
    d_ds = remainder * (q[share] - q[partner])
    d_dx_first = {
        k: q[k] - s * q[share] - (1.0 - s) * q[partner]
        for k in range(x.shape[0]) if k not in (share, partner)
    }
    return s, d_ds, d_dx_first


def quadrature(terms):
    """Add the stated terms in quadrature and say which were stated.

    `terms` maps a label to an array of absolute contributions (already multiplied by their
    partial derivative), broadcastable to a common shape. At each point the total is the root
    sum of squares of the finite terms there, and NaN where no term is finite.

    Returns `(total, included, unstated)`. `included` lists, in the order given, every label
    finite at any point; `unstated` every label NaN at any point. A label finite at some points
    and NaN at others appears in both, which is the honest description of a partly stated term.
    """
    labels = list(terms)
    arrays = np.broadcast_arrays(*[np.asarray(terms[k], dtype="float64") for k in labels])
    stack = np.stack(arrays) if arrays else np.empty((0,))
    finite = np.isfinite(stack)
    squares = np.where(finite, stack, 0.0) ** 2
    total = np.sqrt(squares.sum(axis=0))
    total = np.where(finite.any(axis=0), total, np.nan)
    included = tuple(k for k, f in zip(labels, finite) if f.any())
    unstated = tuple(k for k, f in zip(labels, finite) if not f.all())
    return total, included, unstated
