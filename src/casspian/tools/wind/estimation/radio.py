"""Observation type: radio targets of one anchor. SPEC_19 section 2 item 2, section 3 deliverable 2.

**Data.** Each profile's kind T, the reduction's own input, read at the target's delivered pressures; with
its temperature uncertainty from kind T, or from a file named beside it where its registered kind T carries none
(Lindal's, REVIEW_21_step0). The fit holds none of a source's errors (SPEC_21 section 1 item 5).
**Operator.** `ln T_delivered` from a transfer run to the target; linearized as `int w Delta s dphi` along
the path from the anchor, `w` from `lindal_wind.balance_weight` (SPEC_17).
**Error.** Each target its own error, every target sharing the anchor's: `R` is the targets' blocks with
the anchor's term across them. Each profile's error is its stated uncertainty and, since every profile is an
observation made at one place, the representativeness `sigma_rep(p) = sigma_0 (p_0 / p)^e`, capped, with its own
correlation (SPEC_21 section 1 item 4); IRIS, a zonal mean, carries none. All of these are the control file's. Correlation in `ln p` is Matern 3/2, `(1 + sqrt(3) d/L) exp(-sqrt(3) d/L)`
with `d = |d ln p|` (SPEC_19 v0.3): the exponential expected rough noise that the printed profiles do not
have, so the group ratio could not test sigma_T (REPORT_19_step1 section 5), and this form stays well
conditioned on 66 closely spaced levels, as a Gaussian would not. Each target's diagonal also carries an
floor (SPEC_19 v0.5 section 2 item 2), the control file's `floor_K`, with its basis there. The floor is correlated
in `ln p` over the control file's `floor_correlation_ln_p` in its `floor_correlation_form`, so that a profile's weight
does not depend on how finely it is sampled (SPEC_21 v0.7, after the check 5 diagnosis); with no length given it is
uncorrelated, row by row, as SPEC_19's.
**Layers.** With `average_layers`, every profile is compared at the forward model's resolution: the observed
profile on its own levels, and the delivered profile with its rows and errors (the reference's shared term
included), are averaged in `ln p` into the layers between the forward wind grid's nodes, the grid the stability
constraint uses; one row per layer, the floor once per layer (SPEC_21 v0.7, after the check 5 diagnosis).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

from casspian.lib.constants import BOLTZMANN_CONSTANT
from casspian.tools.lindal import radio_profiles as rp
from casspian.tools.wind.estimation.prior import correlation as prior_correlation
from casspian.tools.wind.estimation.solver import Blocks


@dataclass(frozen=True)
class Profile:
    """A radio profile: its name, its latitude and its kind T, and where its uncertainty is if not in kind T."""
    name: str                #: the group label
    latitude_deg: float      #: planetocentric
    thermo: str              #: its kind T
    uncertainty: str = ""    #: a CSV `pressure_mbar, temperature_uncertainty_K`, where its kind T carries none


@dataclass(frozen=True)
class Target(Profile):
    run: str = ""            #: the base run, whose inputs and namelist the transfer copies


@dataclass(frozen=True)
class Representativeness:
    """`sigma_rep(p) = sigma_0 (p_0 / p)^exponent`, capped, with its correlation in `ln p`."""
    sigma_K: float
    reference_Pa: float
    exponent: float
    cap_K: float
    length: float
    form: str

    def sigma(self, p):
        return np.minimum(self.sigma_K * (self.reference_Pa / np.asarray(p, dtype="float64")) ** self.exponent, self.cap_K)


def read_profile(profile: Profile):
    """`(pressure_Pa, temperature_K, uncertainty_K)` of a profile, increasing pressure."""
    with xr.open_dataset(profile.thermo) as th:
        p = np.asarray(th["pressure_Pa"].values, dtype="float64")
        T = np.asarray(th["temperature_K"].values, dtype="float64")
        s = np.asarray(th["temperature_uncertainty_K"].values, dtype="float64")
    if profile.uncertainty:
        with open(profile.uncertainty, encoding="utf-8") as handle:
            table = np.genfromtxt([line for line in handle if not line.lstrip().startswith("#")], delimiter=",", names=True)
        s = np.interp(np.log(p), np.log(table["pressure_mbar"] * 100.0), table["temperature_uncertainty_K"])
    if not np.all(np.isfinite(s)):
        raise ValueError(f"{profile.name}: its kind T carries no uncertainty at some levels and no file names one")
    order = np.argsort(p)
    return p[order], T[order], s[order]


def correlation(a, b, length, form="matern32"):
    """The correlation between the pressures `a` and `b` over `|d ln p|`, `(a, b)`."""
    return prior_correlation(np.log(a)[:, None] - np.log(b)[None, :], length, form)


def interp_ln_p(p, values, at):
    p, values = np.asarray(p, dtype="float64"), np.asarray(values, dtype="float64")
    order = np.argsort(p)
    return np.interp(np.log(at), np.log(p[order]), values[order], left=np.nan, right=np.nan)


def layer_mean(x, edges):
    """The mean over each layer `[edges[k], edges[k + 1]]` of the piecewise-linear curve through samples at `x`,
    both increasing: `(W, covered)`, `W` `(layers, samples)` with `W @ f` the layer means of `f`, and `covered`
    the layers the samples span entirely. Rows of uncovered layers are zero."""
    x, edges = np.asarray(x, dtype="float64"), np.asarray(edges, dtype="float64")
    W = np.zeros((edges.size - 1, x.size))
    covered = (edges[:-1] >= x[0]) & (edges[1:] <= x[-1])
    for k in np.flatnonzero(covered):
        a, b = edges[k], edges[k + 1]
        for i in range(max(int(np.searchsorted(x, a, side="right")) - 1, 0), min(int(np.searchsorted(x, b)), x.size - 1)):
            lo, hi = max(a, x[i]), min(b, x[i + 1])
            if hi <= lo:
                continue
            h = x[i + 1] - x[i]
            W[k, i] += (x[i + 1] * (hi - lo) - 0.5 * (hi ** 2 - lo ** 2)) / h
            W[k, i + 1] += (0.5 * (hi ** 2 - lo ** 2) - x[i] * (hi - lo)) / h
        W[k] /= b - a
    return W, covered


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
    # The isobars carry `ln N` under the run's anchor's slug (`forward.transfer`), one per run; not always Lindal's.
    names = [n for n in isobars.data_vars if str(n).startswith("ln_refractivity_")]
    if len(names) != 1:
        raise ValueError(f"the run's isobars carry {len(names)} ln_refractivity_ variables, not one: {names}")
    T = labels[:, None] * r_bar[:, None] / (BOLTZMANN_CONSTANT * np.exp(np.asarray(isobars[names[0]].values)))
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

    def __init__(self, reference: Profile, targets, correlation_x, floor_K, representativeness: Representativeness,
                 correlation_form="matern32", floor_correlation_x=0.0, floor_correlation_form="exponential",
                 average_layers=False):
        self.reference = reference
        self.targets = tuple(targets)
        self.anchor_latitude_deg = float(reference.latitude_deg)
        self.correlation_x = float(correlation_x)
        self.correlation_form = correlation_form
        self.floor_K = float(floor_K)
        self.floor_correlation_x = float(floor_correlation_x)
        self.floor_correlation_form = floor_correlation_form
        self.average_layers = bool(average_layers)
        self.shared = None   #: the shared anchor term of the last block, its rows the block's
        self.rep = representativeness
        self.profiles = {p.name: read_profile(p) for p in (reference, *self.targets)}
        self.observed = {name: {"pressure_mbar": v[0] / 100.0, "temperature_K": v[1], "temperature_uncertainty_K": v[2]}
                         for name, v in self.profiles.items()}

    def stated(self, name, p):
        """A profile's stated uncertainty at `p` (K), linear in `ln p`, its edges held."""
        pp, _, s = self.profiles[name]
        return np.interp(np.log(p), np.log(pp), s)

    def delivered(self, basis, products, unit_weight, lat_wind, nodes_Pa=None) -> dict:
        """Each target's delivered profile at one linearization, every level: `{name: (p, ln T, H, labels,
        H_forward)}`, with `H` the rows of `ln T` against the state, `int w Delta s dphi` along the path (SPEC_17).
        Increasing pressure. With the wind grid's pressure nodes, `H_forward` is the same rows with the shear as the
        forward run sees it, linear wind in `ln p` between the nodes (SPEC_21 v0.5 section 1 item 1b); else None."""
        out = {}
        for target in self.targets:
            product = products[target.name]
            root = product.to_dataset(inherit=False)
            p = np.asarray(root["pressure_Pa"].values, dtype="float64")
            T = np.asarray(root["temperature_K"].values, dtype="float64")
            labels, T_wind = temperature_on_wind_grid(root, product["isobars"].to_dataset(inherit=False), lat_wind)
            weight = unit_weight / T_wind
            lat_part = path_integral(lat_wind, weight, basis, self.anchor_latitude_deg, target.latitude_deg)
            H = basis.row(lat_part, basis.x_matrix(labels))
            order = np.argsort(p)
            H_forward = None if nodes_Pa is None else basis.row(lat_part, basis.x_matrix_forward(labels, nodes_Pa))[order]
            out[target.name] = (p[order], np.log(T[order]), H[order], labels[order], H_forward)
        return out

    def block(self, basis, products, unit_weight, lat_wind, delivered=None, nodes_Pa=None) -> Blocks:
        """The term at one linearization. `products` maps a target's name to its run's product;
        `unit_weight` is `balance_weight` with T = 1 on the run's labels. With `average_layers`, `nodes_Pa`, the
        forward wind grid's nodes, bound the layers."""
        delivered = delivered or self.delivered(basis, products, unit_weight, lat_wind)
        if self.average_layers:
            if nodes_Pa is None:
                raise ValueError("average_layers needs the forward wind grid's nodes")
            return self.layer_block(delivered, nodes_Pa)
        ys, hs, rows, groups, sets = [], [], [], [], []
        for target in self.targets:
            p, lnT, H, labels = delivered[target.name][:4]
            pp, TT, _ = self.profiles[target.name]
            T_obs = interp_ln_p(pp, TT, p)
            keep = np.isfinite(T_obs)
            ap, aT, _ = self.profiles[self.reference.name]
            T_a = interp_ln_p(ap, aT, labels)
            pk, lk = p[keep], labels[keep]
            sets.append({"p": pk, "labels": lk,
                         "own": self.stated(target.name, pk) / T_obs[keep],
                         "rep": self.rep.sigma(pk) / T_obs[keep],
                         "anchor": self.stated(self.reference.name, lk) / T_a[keep],
                         "anchor_rep": self.rep.sigma(lk) / T_a[keep],
                         "floor": self.floor_K / T_obs[keep]})
            ys.append(np.log(T_obs[keep]))
            hs.append(lnT[keep])
            rows.append(H[keep])
            groups += [target.name] * int(keep.sum())
        size = sum(s["p"].size for s in sets)
        R = np.zeros((size, size))
        shared = np.zeros((size, size))
        rho = lambda a, b: correlation(a, b, self.correlation_x, self.correlation_form)
        rho_rep = lambda a, b: correlation(a, b, self.rep.length, self.rep.form)

        def floor(a):
            # The floor over the profile's own levels: correlated over a fixed length in ln p, so that its weight per
            # unit ln p is the same on 58 levels as on 307; uncorrelated where no length is given.
            if self.floor_correlation_x <= 0.0:
                return np.diag(a["floor"] ** 2)
            return np.outer(a["floor"], a["floor"]) * correlation(a["p"], a["p"], self.floor_correlation_x,
                                                                  self.floor_correlation_form)

        start_i = 0
        for i, a in enumerate(sets):
            start_j = 0
            for j, b in enumerate(sets):
                # The anchor's error and its representativeness, shared by every target, on the anchor's labels.
                block = (np.outer(a["anchor"], b["anchor"]) * rho(a["labels"], b["labels"])
                         + np.outer(a["anchor_rep"], b["anchor_rep"]) * rho_rep(a["labels"], b["labels"]))
                shared[start_i:start_i + a["p"].size, start_j:start_j + b["p"].size] = block
                if i == j:
                    block = block + (np.outer(a["own"], a["own"]) * rho(a["p"], a["p"])
                                     + np.outer(a["rep"], a["rep"]) * rho_rep(a["p"], a["p"]) + floor(a))
                R[start_i:start_i + a["p"].size, start_j:start_j + b["p"].size] = block
                start_j += b["p"].size
            start_i += a["p"].size
        self.shared = shared
        return Blocks(np.concatenate(ys), np.concatenate(hs), np.vstack(rows), R, np.array(groups))

    def layer_block(self, delivered, nodes_Pa) -> Blocks:
        """The term at the forward model's resolution: each profile averaged in `ln p` into the layers between
        the wind grid's nodes, the observed one on its own levels. Every error part is built on the delivered rows,
        as `block`, and averaged with them; the floor is added once per layer."""
        edges = np.log(np.sort(np.asarray(nodes_Pa, dtype="float64")))
        centers = np.exp(0.5 * (edges[:-1] + edges[1:]))
        held = lambda pp, vv, at: np.interp(np.log(at), np.log(pp), vv)
        ys, hs, rows, groups, sets = [], [], [], [], []
        ap, aT, _ = self.profiles[self.reference.name]
        for target in self.targets:
            p, lnT, H, labels = delivered[target.name][:4]
            pp, TT, _ = self.profiles[target.name]
            W_obs, cov_obs = layer_mean(np.log(pp), edges)
            W, cov_del = layer_mean(np.log(p), edges)
            keep = cov_obs & cov_del
            W_obs, W = W_obs[keep], W[keep]
            y = W_obs @ np.log(TT)
            T_rows, T_a = held(pp, TT, p), held(ap, aT, labels)
            sets.append({"W": W, "p": p, "labels": labels, "center": centers[keep],
                         "own": self.stated(target.name, p) / T_rows,
                         "rep": self.rep.sigma(p) / T_rows,
                         "anchor": self.stated(self.reference.name, labels) / T_a,
                         "anchor_rep": self.rep.sigma(labels) / T_a,
                         "floor": self.floor_K / np.exp(y)})
            ys.append(y)
            hs.append(W @ lnT)
            rows.append(W @ H)
            groups += [target.name] * int(keep.sum())
        size = sum(a["W"].shape[0] for a in sets)
        R = np.zeros((size, size))
        shared = np.zeros((size, size))
        rho = lambda a, b: correlation(a, b, self.correlation_x, self.correlation_form)
        rho_rep = lambda a, b: correlation(a, b, self.rep.length, self.rep.form)

        def floor(a):
            if self.floor_correlation_x <= 0.0:
                return np.diag(a["floor"] ** 2)
            return np.outer(a["floor"], a["floor"]) * correlation(a["center"], a["center"], self.floor_correlation_x,
                                                                  self.floor_correlation_form)

        start_i = 0
        for i, a in enumerate(sets):
            start_j = 0
            for j, b in enumerate(sets):
                # The anchor's error and its representativeness, shared by every target, on the anchor's labels.
                part = (np.outer(a["anchor"], b["anchor"]) * rho(a["labels"], b["labels"])
                        + np.outer(a["anchor_rep"], b["anchor_rep"]) * rho_rep(a["labels"], b["labels"]))
                if i == j:
                    part = part + (np.outer(a["own"], a["own"]) * rho(a["p"], a["p"])
                                   + np.outer(a["rep"], a["rep"]) * rho_rep(a["p"], a["p"]))
                block = a["W"] @ part @ b["W"].T
                if i != j:
                    shared[start_i:start_i + a["W"].shape[0], start_j:start_j + b["W"].shape[0]] = block
                else:
                    shared[start_i:start_i + a["W"].shape[0], start_j:start_j + b["W"].shape[0]] = a["W"] @ (
                        np.outer(a["anchor"], a["anchor"]) * rho(a["labels"], a["labels"])
                        + np.outer(a["anchor_rep"], a["anchor_rep"]) * rho_rep(a["labels"], a["labels"])) @ a["W"].T
                    block = block + floor(a)
                R[start_i:start_i + a["W"].shape[0], start_j:start_j + b["W"].shape[0]] = block
                start_j += b["W"].shape[0]
            start_i += a["W"].shape[0]
        self.shared = shared
        return Blocks(np.concatenate(ys), np.concatenate(hs), np.vstack(rows), R, np.array(groups))

    def shape(self, blocks: Blocks, group, products):
        """The correlation shape a systematic term in this group would take."""
        root = products[group].to_dataset(inherit=False)
        p = np.asarray(root["pressure_Pa"].values, dtype="float64")
        pp, TT, _ = self.profiles[group]
        keep = np.isfinite(interp_ln_p(pp, TT, p))
        return correlation(p[keep], p[keep], self.correlation_x, self.correlation_form)
