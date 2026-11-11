"""Observation type: temperature gradients along isobars. SPEC_19 section 2 item 3, section 3 deliverable 2.

**Data.** `d ln T / dphi_c` from `iris_temperatures.fit` of IRIS-like data at a level, the derivative of
the fitted value with the smooth window as the base run's shear uses it (SPEC_13), sampled every 4 deg. A row
where no data point lies within `gap_deg` is not made: there the fit only bridges, and a bridge is not a
measurement (REVIEW_21_step0 ruling 7). Nothing else is omitted: near the equator the operator's own weight,
which `w` drives toward zero, says how little a row carries (SPEC_21 section 1 item 9).
**Operator.** The model's balance, `(d ln T / dphi)_p = w s + a`, with
`a = -(2 Omega_abs cos(phi) / g) (du/dphi)_p` (`lindal_wind` module docstring), averaged uniformly in `ln p`
over a layer centered on the level: IRIS's values "represent mean values for atmospheric layers
approximately one-half scale height thick" (Conrath and Pirraglia 1983, Fig. 1; SPEC_19 v0.3). Linear in the
shear with the wind held, as `balance_weight` is; `T` is the fit's.
**Error.** The standard error of the fit's gradient (SPEC_10), independent row to row; and, as the general rules
of the design note v0.6 section 4 make them, two correlated parts:
* **An uncertain level.** Where a level's placement is uncertain (IRIS's top level, 110 mbar to Lindal and 150 to
  Conrath and Pirraglia), one label serves every row at that level, so the error is rank one: the vector `d` of
  the operator's change across the range, scaled for a uniform range, `d d^T / 12`, recomputed each pass
  (REVIEW_21_step0 ruling 8). The range is the data's property, read at ingest.
* **A systematic part,** correlated over a long distance in latitude, sized by the consistency check after a
  first fit (SPEC_21 section 1 item 7); its size, length and form are the fit's.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from casspian.lib import geoid as gd
from casspian.lib.gravity import g_eff_radial, omega_abs
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.wind.estimation.prior import correlation
from casspian.tools.wind.estimation.solver import Blocks

#: The layer's thickness in `ln p`, half a scale height (Conrath and Pirraglia 1983, Fig. 1).
LAYER_THICKNESS = 0.5
#: Gauss-Legendre in `ln p`: the shear is a cubic spline with knots about as far apart as the layer is thick.
_NODES, _QUADRATURE = np.polynomial.legendre.leggauss(5)
LAYER_WEIGHTS = _QUADRATURE / 2.0


def layer_pressures(pressure_Pa) -> np.ndarray:
    """The pressures at which the layer centered on `pressure_Pa` is sampled, with `LAYER_WEIGHTS`."""
    return float(pressure_Pa) * np.exp(0.5 * LAYER_THICKNESS * _NODES)


@dataclass(frozen=True)
class Level:
    name: str                 #: the group label
    pressure_Pa: float
    indices: np.ndarray       #: the points, as indices of the wind's latitude grid
    y: np.ndarray             #: d ln T / dphi_c, per radian
    sigma: np.ndarray
    temperature_K: np.ndarray
    range_Pa: tuple = None    #: the pressures its placement may take, `(low, high)`, or None


def read_levels(table, levels, lat_c, lat_g, sample_deg, gap_deg, ranges=None) -> list:
    """`levels` maps a printed level (mbar) to its pressure (Pa). `lat_c`, `lat_g` are the wind's grid. `ranges`
    maps a printed level to the `(low, high)` pressures (Pa) its placement may take."""
    per_rad = np.gradient(lat_g, lat_c) * (180.0 / np.pi)
    ranges = ranges or {}
    out = []
    on_sample = np.abs(lat_c / sample_deg - np.round(lat_c / sample_deg)) < 1e-9
    for printed, pressure in levels.items():
        rows = table["pressure_mbar"] == printed
        data = table["latitude_planetographic_deg"][rows]
        result = iris.fit(data, table["temperature_K"][rows], grid=lat_g, window="smooth")
        slope = lw.derivative_of_value(result.value_K, lat_g)
        near = np.min(np.abs(lat_g[:, None] - data[None, :]), axis=1) <= gap_deg
        use = on_sample & near & np.isfinite(slope) & (np.abs(lat_c) < 90.0)
        T = result.value_K[use]
        span = ranges.get(printed)
        out.append(Level(f"IRIS {printed} mbar", float(pressure), np.flatnonzero(use),
                         slope[use] * per_rad[use] / T,
                         result.gradient_error_K_per_deg[use] * per_rad[use] / T, T,
                         None if span is None else (float(span[0]), float(span[1]))))
    return out


def wind_term(wind, gravity, rotation, pressure_Pa) -> np.ndarray:
    """`a` on the wind's latitude grid at one level: the balance's part that the wind itself carries."""
    lat = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    phi = np.radians(lat)
    grid = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    rising = np.argsort(np.log(grid))
    u_grid = np.asarray(wind["u_total_ms"].values, dtype="float64")
    u = np.array([np.interp(np.log(pressure_Pa), np.log(grid[rising]), row[rising]) for row in u_grid])
    Omega = float(rotation["angular_rate_rad_s"])
    GM, R_norm = float(gravity["GM_m3s2"]), float(gravity["normalization_radius_m"])
    J, degrees = np.asarray(gravity["J"].values, dtype="float64"), np.asarray(gravity["degree"].values)
    r = np.asarray(gd.wind_geoid(phi, R_norm, "equatorial_radius", lambda x, u=u: np.interp(x, phi, u),
                                 Omega, GM, J, degrees, R_norm).radius, dtype="float64")
    g = g_eff_radial(u, r, phi, Omega, GM, J, degrees, R_norm)
    return -2.0 * omega_abs(u, r, phi, Omega) * np.cos(phi) / g * np.gradient(u, phi)


@dataclass(frozen=True)
class Systematic:
    """A level's systematic part: its size in `d ln T / dphi` per radian, and its correlation in latitude."""
    sigma: float
    length_deg: float
    form: str


class TemperatureGradients:
    def __init__(self, levels, systematic=None):
        self.levels = list(levels)
        self.systematic = dict(systematic or {})

    def placements(self):
        """Every pressure a layer is centered on: each level's own, and its range's ends."""
        out = []
        for level in self.levels:
            out.append(level.pressure_Pa)
            if level.range_Pa is not None:
                out += list(level.range_Pa)
        return sorted(set(out))

    @staticmethod
    def _operator(basis, state, level, lat_part, at):
        """`H` and `h` of one level's rows with the layer centered at the pressure `at` describes."""
        unit, terms, pressure = at
        x_parts = basis.x_matrix(layer_pressures(pressure))
        n = level.indices.size
        H = np.zeros((n, basis.size))
        term = np.zeros(n)
        for q, quadrature in enumerate(LAYER_WEIGHTS):
            w = unit[q, level.indices] / level.temperature_K
            H += quadrature * basis.row(lat_part * w[:, None], np.repeat(x_parts[q:q + 1], n, axis=0))
            term += quadrature * terms[q, level.indices]
        return H, H @ state + term

    def block(self, basis, state, at, lat_wind) -> Blocks:
        """The term at one linearization. `at` maps a layer's center (Pa) to `(unit_weight, wind_terms, pressure)`,
        `balance_weight` with T = 1 and the matching `a` at its `layer_pressures`, each `(node, latitude)`."""
        ys, hs, rows, groups, Rs = [], [], [], [], []
        for level in self.levels:
            lat = lat_wind[level.indices]
            lat_part = basis.lat_matrix(lat)
            H, h = self._operator(basis, state, level, lat_part, at[level.pressure_Pa])
            R = np.diag(level.sigma ** 2)
            if level.range_Pa is not None:
                _, low = self._operator(basis, state, level, lat_part, at[level.range_Pa[0]])
                _, high = self._operator(basis, state, level, lat_part, at[level.range_Pa[1]])
                d = high - low
                R = R + np.outer(d, d) / 12.0
            part = self.systematic.get(level.name)
            if part is not None and part.sigma > 0.0:
                R = R + part.sigma ** 2 * correlation(lat[:, None] - lat[None, :], part.length_deg, part.form)
            ys.append(level.y)
            hs.append(h)
            rows.append(H)
            groups += [level.name] * lat.size
            Rs.append(R)
        n = sum(r.shape[0] for r in Rs)
        R = np.zeros((n, n))
        start = 0
        for r in Rs:
            R[start:start + r.shape[0], start:start + r.shape[0]] = r
            start += r.shape[0]
        return Blocks(np.concatenate(ys), np.concatenate(hs), np.vstack(rows), R, np.array(groups))

    def shape(self, group, lat_wind, length_deg, form):
        """The latitude correlation a systematic part in this group would take."""
        level = next(v for v in self.levels if v.name == group)
        lat = lat_wind[level.indices]
        return correlation(lat[:, None] - lat[None, :], length_deg, form)
