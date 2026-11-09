"""Observation type: radio targets of one anchor. SPEC_19 section 2 item 2, section 3 deliverable 2.

**Data.** A profile `pressure_mbar, temperature_K` (`tools/lindal/radio_profiles.read_csv`), read at the
target's delivered pressures.
**Operator.** `ln T_delivered` from a transfer run to the target; linearized as `int w Delta s dphi` along
the path from the anchor, `w` from `lindal_wind.balance_weight` (SPEC_17).
**Error.** Each target its own error, every target sharing the anchor's: `R` is the targets' blocks with
the anchor's term across them. Correlation in `ln p` is Matern 3/2, `(1 + sqrt(3) d/L) exp(-sqrt(3) d/L)`
with `d = |d ln p|` (SPEC_19 v0.3): the exponential expected rough noise that the printed profiles do not
have, so the group ratio could not test sigma_T (REPORT_19_step1 section 5), and this form stays well
conditioned on 66 closely spaced levels, as a Gaussian would not. Each target's diagonal also carries an
uncorrelated floor (SPEC_19 v0.5 section 2 item 2), the control file's `floor_K`, with its basis there.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from casspian.lib.constants import BOLTZMANN_CONSTANT
from casspian.tools.lindal import radio_profiles as rp
from casspian.tools.wind.estimation.solver import Blocks


@dataclass(frozen=True)
class Target:
    name: str                #: the group label
    run: str                 #: the base run, whose inputs and namelist the transfer copies
    latitude_deg: float      #: planetocentric
    observed: str            #: the data file
    extra_K: tuple = ()      #: added errors: (kelvins, top Pa, bottom Pa)


def correlation(a, b, length):
    """Matern 3/2 between the pressures `a` and `b`, `(a, b)`."""
    d = np.sqrt(3.0) * np.abs(np.log(a)[:, None] - np.log(b)[None, :]) / length
    return (1.0 + d) * np.exp(-d)


def interp_ln_p(p, values, at):
    p, values = np.asarray(p, dtype="float64"), np.asarray(values, dtype="float64")
    order = np.argsort(p)
    return np.interp(np.log(at), np.log(p[order]), values[order], left=np.nan, right=np.nan)


def sigma_curve(path):
    """Lindal's sigma_T against pressure from a profile that carries it at every row (Fig. 4)."""
    table = rp.read_csv(path)
    good = np.isfinite(table["temperature_uncertainty_K"])
    p, s = table["pressure_mbar"][good] * 100.0, table["temperature_uncertainty_K"][good]
    order = np.argsort(p)
    return lambda at: np.interp(np.log(at), np.log(p[order]), s[order])


def temperature_on_wind_grid(root, isobars, lat_wind):
    """T on each isobar of a run, `(level, wind latitude)`, NaN off the path: the isobars' `ln N` with the
    target's mean refractivity and the isobar's label (SPEC_17)."""
    lat_m = np.asarray(isobars["latitude_planetocentric_deg"].values)
    labels = np.asarray(root["pressure_label_Pa"].values)
    r_bar = np.asarray(root["mean_refractivity_m3"].values)
    T = labels[:, None] * r_bar[:, None] / (BOLTZMANN_CONSTANT * np.exp(np.asarray(isobars["ln_refractivity_lindal"].values)))
    step = float(np.max(np.diff(lat_wind)))
    off = (lat_wind < lat_m.min() - step) | (lat_wind > lat_m.max() + step)
    out = np.stack([np.interp(lat_wind, lat_m, row) for row in T])
    out[:, off] = np.nan
    return labels, out


def path_integral(lat_wind, weight, basis, start, end):
    """`int_start^end w(phi) B_i(phi) dphi` (radians), `(level, n_lat)`, oriented from `start` to `end`.
    `weight` is `(level, wind latitude)`; the path's ends are nodes of the integral."""
    lo, hi = min(start, end), max(start, end)
    nodes = np.concatenate([[lo], lat_wind[(lat_wind > lo) & (lat_wind < hi)], [hi]])
    w_nodes = np.stack([np.interp(nodes, lat_wind, row) for row in weight])
    b = basis.lat_matrix(nodes)
    out = np.trapezoid(w_nodes[:, :, None] * b[None, :, :], np.radians(nodes), axis=1)
    return out if end > start else -out


class RadioTargets:
    """The targets of one anchor, and their errors."""

    def __init__(self, targets, anchor_latitude_deg, sigma_path, correlation_x, floor_K):
        self.targets = tuple(targets)
        self.anchor_latitude_deg = float(anchor_latitude_deg)
        self.sigma = sigma_curve(sigma_path)
        self.correlation_x = float(correlation_x)
        self.floor_K = float(floor_K)
        self.observed = {t.name: rp.read_csv(t.observed) for t in self.targets}

    def block(self, basis, products, unit_weight, lat_wind) -> Blocks:
        """The term at one linearization. `products` maps a target's name to its run's product;
        `unit_weight` is `balance_weight` with T = 1 on the run's labels."""
        ys, hs, rows, groups, sets = [], [], [], [], []
        for target in self.targets:
            product = products[target.name]
            root = product.to_dataset(inherit=False)
            p = np.asarray(root["pressure_Pa"].values, dtype="float64")
            T = np.asarray(root["temperature_K"].values, dtype="float64")
            table = self.observed[target.name]
            T_obs = interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], p)
            labels, T_wind = temperature_on_wind_grid(root, product["isobars"].to_dataset(inherit=False), lat_wind)
            weight = unit_weight / T_wind
            lat_part = path_integral(lat_wind, weight, basis, self.anchor_latitude_deg, target.latitude_deg)
            x_part = basis.x_matrix(labels)
            H = basis.row(lat_part, x_part)
            keep = np.isfinite(T_obs)
            own = self.sigma(p) ** 2
            for size, top, bottom in target.extra_K:
                own = own + np.where((p >= top) & (p <= bottom), float(size) ** 2, 0.0)
            T_anchor = np.asarray(product["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)["temperature_K"].values)
            anchor_p = np.asarray(product["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)["pressure_Pa"].values)
            T_a = interp_ln_p(anchor_p, T_anchor, labels)
            sets.append({"p": p[keep], "own": np.sqrt(own[keep]) / T_obs[keep], "anchor": self.sigma(labels[keep]) / T_a[keep],
                         "floor": self.floor_K / T_obs[keep]})
            ys.append(np.log(T_obs[keep]))
            hs.append(np.log(T[keep]))
            rows.append(H[keep])
            groups += [target.name] * int(keep.sum())
        size = sum(s["p"].size for s in sets)
        R = np.zeros((size, size))
        rho = lambda a, b: correlation(a, b, self.correlation_x)
        start_i = 0
        for i, a in enumerate(sets):
            start_j = 0
            for j, b in enumerate(sets):
                block = np.outer(a["anchor"], b["anchor"]) * rho(a["p"], b["p"])
                if i == j:
                    block = block + np.outer(a["own"], a["own"]) * rho(a["p"], a["p"]) + np.diag(a["floor"] ** 2)
                R[start_i:start_i + a["p"].size, start_j:start_j + b["p"].size] = block
                start_j += b["p"].size
            start_i += a["p"].size
        return Blocks(np.concatenate(ys), np.concatenate(hs), np.vstack(rows), R, np.array(groups))

    def shape(self, blocks: Blocks, group, products):
        """The correlation shape a systematic term in this group would take."""
        target = next(t for t in self.targets if t.name == group)
        root = products[target.name].to_dataset(inherit=False)
        p = np.asarray(root["pressure_Pa"].values, dtype="float64")
        table = self.observed[target.name]
        keep = np.isfinite(interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], p))
        return correlation(p[keep], p[keep], self.correlation_x)
