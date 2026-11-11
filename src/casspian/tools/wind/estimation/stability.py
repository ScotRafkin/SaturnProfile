"""A constraint, not an observation: no delivered profile less stable than its source. SPEC_21 v0.5 section 1
item 1b.

**Layers at the forward run's resolution.** A layer is one interval of the wind grid's pressure nodes. The
forward run holds the wind linear in `ln p` between nodes (`lib.windfield`, decision L), so the shear is constant
across an interval and two levels inside one move together: a slope between close levels is one the fit cannot
change (REVIEW_21_step1_preliminary_v04). Each target's delivered `ln T` is taken at the nodes, linear in `ln p`
between its levels, and a layer's slope is the difference over the interval. Its row is formed the same way from
`H_forward`, the delivered rows with the shear as the forward run sees it (`Basis.x_matrix_forward`).

**Relative to the source.** The excess over the dry adiabat, `e = d ln T / d ln p - kappa` with `kappa = R / c_p`
for normal hydrogen, `c_p` from the run's composition at the layer's mean temperature (`lib.heat_capacity`), is
held to no more than the source anchor's own excess in the same layer: `z = e - max(e_anchor, 0)`. Excess already
in the anchor is part of the observation and no wind can remove it; it is reported, not penalized.

**A smooth one-sided penalty.** Every layer enters, always, as a row `y = 0`, `h = s(z)`, `H = s'(z) dz/dstate`,
with `s(z) = width * log(1 + exp(z / width))`, a softened hinge, and the control file's `sigma`: far below the
allowed excess a layer adds nothing, above it the penalty is `z / sigma`, and in between it turns on smoothly, so no
layer switches the iteration on or off. `kappa`'s own dependence on temperature is left out of the row.

Not a measurement, so Huber's rule and the consistency check leave it out (`GROUP`). A layer is reported active
where `z > 0`.
"""

from __future__ import annotations

import numpy as np
from scipy.special import expit

from casspian.lib import heat_capacity as hc
from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.tools.wind.estimation.solver import Blocks

GROUP = "stability"


def _at_nodes(p, values, nodes):
    """Linear in `ln p` from levels `p` (increasing) to `nodes` inside them: `(values at nodes, weights)`."""
    x, xn = np.log(p), np.log(nodes)
    k = np.clip(np.searchsorted(x, xn, side="right") - 1, 0, x.size - 2)
    a = (xn - x[k]) / (x[k + 1] - x[k])
    W = np.zeros((nodes.size, p.size))
    W[np.arange(nodes.size), k] = 1.0 - a
    W[np.arange(nodes.size), k + 1] = a
    return W @ values, W


class Stability:
    def __init__(self, sigma, width, composition, latitudes_deg: dict, anchor, anchor_latitude_deg, nodes_Pa):
        self.sigma = float(sigma)
        self.width = float(width)
        self.composition = composition
        self.latitudes_deg = dict(latitudes_deg)
        self.anchor_p, self.anchor_lnT = np.asarray(anchor[0], dtype="float64"), np.log(np.asarray(anchor[1], dtype="float64"))
        self.anchor_latitude_deg = float(anchor_latitude_deg)
        self.nodes = np.sort(np.asarray(nodes_Pa, dtype="float64"))
        self.active = {}
        self.source = {}

    def kappa(self, p_lo, p_hi, lnT_lo, lnT_hi, latitude_deg):
        mid_p, mid_T = np.sqrt(p_lo * p_hi), np.exp(0.5 * (lnT_lo + lnT_hi))
        return MOLAR_GAS_CONSTANT / hc.mixture_cp(self.composition, mid_p, latitude_deg, mid_T, "normal")

    def layers(self, p, lnT):
        """The nodes inside both the delivered profile and the anchor, and each layer's `(slope, excess)` for the
        delivered profile and for the anchor."""
        nodes = self.nodes[(self.nodes >= max(p[0], self.anchor_p[0])) & (self.nodes <= min(p[-1], self.anchor_p[-1]))]
        dx = np.diff(np.log(nodes))
        return nodes, dx

    def anchor_excess(self, nodes, dx):
        lnT_a, _ = _at_nodes(self.anchor_p, self.anchor_lnT, nodes)
        slope = np.diff(lnT_a) / dx
        return slope - self.kappa(nodes[:-1], nodes[1:], lnT_a[:-1], lnT_a[1:], self.anchor_latitude_deg)

    def block(self, delivered: dict, n_state: int) -> Blocks:
        """Every layer's row at this linearization; `delivered` is `RadioTargets.delivered`'s with the wind nodes."""
        ys, hs, rows = [], [], []
        self.active, self.source = {}, {}
        for name, (p, lnT, _, _, H) in delivered.items():
            if H is None:
                raise ValueError("the stability constraint needs the delivered rows on the wind grid (H_forward)")
            good = np.isfinite(lnT)
            p, lnT, H = p[good], lnT[good], H[good]
            nodes, dx = self.layers(p, lnT)
            lnT_n, W = _at_nodes(p, lnT, nodes)
            H_n = W @ H
            slope = np.diff(lnT_n) / dx
            excess = slope - self.kappa(nodes[:-1], nodes[1:], lnT_n[:-1], lnT_n[1:], self.latitudes_deg[name])
            e_a = self.anchor_excess(nodes, dx)
            z = excess - np.maximum(e_a, 0.0)
            u = z / self.width
            s = self.width * np.logaddexp(0.0, u)
            ds = expit(u)
            on = z > 0
            self.active[name] = {"p_lo": nodes[:-1][on], "p_hi": nodes[1:][on], "z": z[on], "penalty": s[on]}
            src = e_a > 0
            self.source[name] = {"p_lo": nodes[:-1][src], "p_hi": nodes[1:][src], "anchor": e_a[src], "delivered": excess[src]}
            ys.append(np.zeros(z.size))
            hs.append(s)
            rows.append(ds[:, None] * (H_n[1:] - H_n[:-1]) / dx[:, None])
        y = np.concatenate(ys)
        return Blocks(y, np.concatenate(hs), np.vstack(rows), np.eye(y.size) * self.sigma ** 2, np.array([GROUP] * y.size))

    def report(self) -> str:
        parts = []
        for name, a in self.active.items():
            if a["z"].size:
                parts.append(f"{name}: " + ", ".join(f"{lo / 100:.0f} to {hi / 100:.0f} mbar +{z:.4f}"
                                                     for lo, hi, z in zip(a["p_lo"], a["p_hi"], a["z"])))
        return "; ".join(parts) if parts else "none"

    def report_source(self) -> str:
        """Excess already in the anchor, reported not penalized: the anchor's and the delivered profile's."""
        parts = []
        for name, a in self.source.items():
            if a["anchor"].size:
                parts.append(f"{name}: " + ", ".join(f"{lo / 100:.0f} to {hi / 100:.0f} mbar anchor +{ea:.4f} delivered {ed:+.4f}"
                                                     for lo, hi, ea, ed in zip(a["p_lo"], a["p_hi"], a["anchor"], a["delivered"])))
        return "; ".join(parts) if parts else "none"
