"""The state: the shear `s(phi, p)` at one epoch on a smooth basis. SPEC_19 section 1.

Cubic B-splines in planetocentric latitude and in `ln p`, clamped at both ends. Above and below the
`ln p` knots the shear is held at its edge value. At the poles it is zero: the two polar latitude
functions are left out, since a kind W file's wind is exactly zero there (SPEC_00 section 6.6) and the
reference wind already is. The state is a coefficient array `(latitude, ln p)`, flattened when the
solver sees it. The wind follows as `u = u_ref + int s dln p` from the reference level.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import BSpline

DEGREE = 3


def _clamped(knots):
    knots = np.asarray(knots, dtype="float64")
    return np.concatenate([[knots[0]] * DEGREE, knots, [knots[-1]] * DEGREE])


def _design(t, x):
    """Every basis function of the clamped knot vector `t` at `x`, `(point, function)`."""
    x = np.clip(np.atleast_1d(np.asarray(x, dtype="float64")), t[0], t[-1])
    return BSpline.design_matrix(x, t, DEGREE, extrapolate=False).toarray()


@dataclass(frozen=True)
class Basis:
    """The two one-dimensional bases, and the epoch the state belongs to."""

    epoch: str
    t_lat: np.ndarray
    t_x: np.ndarray
    reference_Pa: float

    @classmethod
    def build(cls, epoch, lat_step_deg, x_step, bottom_Pa, top_Pa, reference_Pa):
        lat = np.linspace(-90.0, 90.0, int(round(180.0 / lat_step_deg)) + 1)
        span = np.log(bottom_Pa / top_Pa)
        x = np.linspace(np.log(top_Pa), np.log(bottom_Pa), int(round(span / x_step)) + 1)
        return cls(str(epoch), _clamped(lat), _clamped(x), float(reference_Pa))

    @property
    def n_lat(self):
        return self.t_lat.size - DEGREE - 1 - 2

    @property
    def n_x(self):
        return self.t_x.size - DEGREE - 1

    @property
    def size(self):
        return self.n_lat * self.n_x

    def lat_matrix(self, latitude_deg):
        """The free latitude functions at `latitude_deg`, `(point, n_lat)`."""
        return _design(self.t_lat, latitude_deg)[:, 1:-1]

    def x_matrix(self, pressure_Pa):
        """The `ln p` functions at `pressure_Pa`, held at the edges, `(point, n_x)`."""
        return _design(self.t_x, np.log(pressure_Pa))

    def x_integral(self, pressure_Pa):
        """`int_{ln p_ref}^{ln p}` of each `ln p` function, the held edges included, `(point, n_x)`."""
        x = np.log(np.atleast_1d(np.asarray(pressure_Pa, dtype="float64")))
        low, high = self.t_x[0], self.t_x[-1]
        inside = np.clip(x, low, high)
        ref = np.clip(np.log(self.reference_Pa), low, high)
        out = np.empty((x.size, self.n_x))
        edge_low, edge_high = _design(self.t_x, low)[0], _design(self.t_x, high)[0]
        for j in range(self.n_x):
            c = np.zeros(self.n_x)
            c[j] = 1.0
            anti = BSpline(self.t_x, c, DEGREE, extrapolate=False).antiderivative()
            out[:, j] = (anti(inside) - anti(ref)
                         + edge_low[j] * np.minimum(x - low, 0.0) + edge_high[j] * np.maximum(x - high, 0.0))
        return out

    def x_matrix_forward(self, pressure_Pa, nodes_Pa):
        """The shear as the forward run sees it, `(point, n_x)`: the wind on the wind grid's pressure nodes, linear
        in `ln p` between them, so `du / dln p` is the node difference over the cell (`lib.windfield`, decision L).
        A point on a node takes the cell to higher pressure, as `WindField._cell` does. SPEC_21 v0.5 section 1
        item 1b."""
        nodes = np.sort(np.asarray(nodes_Pa, dtype="float64"))
        x_nodes = np.log(nodes)
        k = np.clip(np.searchsorted(x_nodes, np.log(np.atleast_1d(pressure_Pa)), side="right") - 1, 0, nodes.size - 2)
        integral = self.x_integral(nodes)
        return (integral[k + 1] - integral[k]) / (x_nodes[k + 1] - x_nodes[k])[:, None]

    def greville(self):
        """The coefficients' abscissae, where the prior's distances are taken: (latitude, ln p)."""
        lat = np.array([self.t_lat[i + 1:i + DEGREE + 1].mean() for i in range(self.t_lat.size - DEGREE - 1)])
        x = np.array([self.t_x[i + 1:i + DEGREE + 1].mean() for i in range(self.n_x)])
        return lat[1:-1], x

    def row(self, latitude_matrix, x_matrix):
        """Rows of the operator from `(m, n_lat)` and `(m, n_x)`, one point each: `(m, size)`."""
        return np.einsum("mi,mj->mij", latitude_matrix, x_matrix).reshape(latitude_matrix.shape[0], -1)

    def shear(self, state, latitude_deg, pressure_Pa):
        """`s` at the grid `(pressure, latitude)`."""
        c = np.asarray(state, dtype="float64").reshape(self.n_lat, self.n_x)
        return self.x_matrix(pressure_Pa) @ c.T @ self.lat_matrix(latitude_deg).T

    def wind_change(self, state, latitude_deg, pressure_Pa):
        """`int s dln p` from the reference level at the grid `(latitude, pressure)`."""
        c = np.asarray(state, dtype="float64").reshape(self.n_lat, self.n_x)
        return self.lat_matrix(latitude_deg) @ c @ self.x_integral(pressure_Pa).T

    def project(self, shear, latitude_deg, pressure_Pa, ridge=1e-8):
        """The state nearest a shear `(pressure, latitude)` on a grid, by least squares."""
        a = self.lat_matrix(latitude_deg)
        b = self.x_matrix(pressure_Pa)
        left = np.linalg.solve(a.T @ a + ridge * np.eye(a.shape[1]), a.T)
        right = np.linalg.solve(b.T @ b + ridge * np.eye(b.shape[1]), b.T)
        return (left @ np.asarray(shear).T @ right.T).ravel()
