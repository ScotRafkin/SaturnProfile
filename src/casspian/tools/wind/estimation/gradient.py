"""Observation type: temperature gradients along isobars. SPEC_19 section 2 item 3, section 3 deliverable 2.

**Data.** `d ln T / dphi_c` from `iris_temperatures.fit` of IRIS-like data at a level, the derivative of
the fitted value with the smooth window as the base run's shear uses it (SPEC_13), sampled every 4 deg.
Omitted within 5 deg of the equator and where no data point lies within 2 deg.
**Operator.** The model's balance, `(d ln T / dphi)_p = w s + a`, with
`a = -(2 Omega_abs cos(phi) / g) (du/dphi)_p` (`lindal_wind` module docstring), averaged uniformly in `ln p`
over a layer centered on the level: IRIS's values "represent mean values for atmospheric layers
approximately one-half scale height thick" (Conrath and Pirraglia 1983, Fig. 1; SPEC_19 v0.3). Linear in the
shear with the wind held, as `balance_weight` is; `T` is the fit's.
**Error.** Random and independent: the standard error of the fit's gradient (SPEC_10).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from casspian.lib import geoid as gd
from casspian.lib.gravity import g_eff_radial, omega_abs
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.lindal import lindal_wind as lw
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


def read_levels(table, levels, lat_c, lat_g, sample_deg, equatorial_deg, gap_deg) -> list:
    """`levels` maps a printed level (mbar) to its pressure (Pa). `lat_c`, `lat_g` are the wind's grid."""
    per_rad = np.gradient(lat_g, lat_c) * (180.0 / np.pi)
    out = []
    on_sample = np.abs(lat_c / sample_deg - np.round(lat_c / sample_deg)) < 1e-9
    for printed, pressure in levels.items():
        rows = table["pressure_mbar"] == printed
        data = table["latitude_planetographic_deg"][rows]
        result = iris.fit(data, table["temperature_K"][rows], grid=lat_g, window="smooth")
        slope = lw.derivative_of_value(result.value_K, lat_g)
        near = np.min(np.abs(lat_g[:, None] - data[None, :]), axis=1) <= gap_deg
        use = on_sample & (np.abs(lat_c) >= equatorial_deg) & near & np.isfinite(slope) & (np.abs(lat_c) < 90.0)
        T = result.value_K[use]
        out.append(Level(f"IRIS {printed} mbar", float(pressure), np.flatnonzero(use),
                         slope[use] * per_rad[use] / T,
                         result.gradient_error_K_per_deg[use] * per_rad[use] / T, T))
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


class TemperatureGradients:
    def __init__(self, levels):
        self.levels = list(levels)

    def block(self, basis, state, unit_weight, wind_terms, lat_wind) -> Blocks:
        """The term at one linearization. `unit_weight` is `balance_weight` with T = 1 at each level's
        `layer_pressures`, `(level, node, latitude)`, and `wind_terms` the matching `a`."""
        ys, hs, rows, groups, sigmas = [], [], [], [], []
        for k, level in enumerate(self.levels):
            lat = lat_wind[level.indices]
            lat_part = basis.lat_matrix(lat)
            x_parts = basis.x_matrix(layer_pressures(level.pressure_Pa))
            H = np.zeros((lat.size, basis.size))
            term = np.zeros(lat.size)
            for q, quadrature in enumerate(LAYER_WEIGHTS):
                w = unit_weight[k, q, level.indices] / level.temperature_K
                H += quadrature * basis.row(lat_part * w[:, None], np.repeat(x_parts[q:q + 1], lat.size, axis=0))
                term += quadrature * wind_terms[k, q, level.indices]
            ys.append(level.y)
            hs.append(H @ state + term)
            rows.append(H)
            groups += [level.name] * lat.size
            sigmas.append(level.sigma)
        return Blocks(np.concatenate(ys), np.concatenate(hs), np.vstack(rows),
                      np.diag(np.concatenate(sigmas) ** 2), np.array(groups))
