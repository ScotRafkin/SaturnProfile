"""The prior: a mean and a covariance. Design note v0.6 section 5, SPEC_21 section 1 item 1.

The mean is zero shear here and carries no observation; a climatology for the date enters as the mean.

**Stated in temperature (SPEC_21, REVIEW_21_step0).** The input is a temperature-anomaly covariance on
isobars: an amplitude `sigma_T(p)` and correlations in latitude and `ln p`, separable. The model's own
linearized balance, `d ln T / dphi = w s` on an isobar (`lindal_wind.balance_weight`, SPEC_17), carries it to
the shear:

1. The gradient `g = (1 / T_bar) dT'/dphi` has, in latitude, the covariance `-rho''` of the temperature's
   correlation `rho`; for Matern 5/2 of length `l`, with `q = sqrt(5) |dphi| / l` in radians,
   `(5 / 3 l^2)(1 + q - q^2) exp(-q)`.
2. It is conditioned on zero at the equator, where the balance requires it for any finite shear, and at both
   poles, where the basis holds the shear at zero and a smooth field has no gradient (rulings 2).
3. `s = g / w`, with `w` held at the base run's wind and the anchor's `T_bar(p)`, uniform in latitude
   (ruling 3), so the prior does not move with the iterate.
4. It is formed on a fine grid that does not hold the equator or the poles as nodes, where `w` or the
   conditioned covariance vanishes: the nearest nodes, a quarter step away, carry the finite limit
   `s(0) = g'(0) / w'(0)` to within the grid's resolution. It is carried onto the basis by the basis's own
   least-squares projection, `B = P B_s P^T`, which is dense, since `w` breaks the separable form.

**A deformation-radius floor (SPEC_21 v0.4 section 1 item 1a).** A balanced structure cannot be narrower than
the deformation radius; narrower ones are waves. The latitude length at each latitude and level is the larger of
the input length and the local deformation radius `L_R = c / sqrt(f^2 + beta c)`, `f = 2 Omega sin(phi)`,
`beta = 2 Omega / a`: `c / |f|` away from the equator and `sqrt(c a / 2 Omega)` at it, the standard form joining
the two. `c = N h` with `h` the prior's vertical length in height, from the anchor's profile (the driver's).
A varying length is made by warping latitude: `psi(phi) = int_0^phi l_0 / l(phi') dphi'`, with the stationary
correlation of length `l_0` in `psi`. A warped positive-definite kernel is positive definite, also when the warp
changes with level, which is why this form and not a varying length inside the stationary formula. The
gradient's covariance stays exact: `-rho''(psi - psi') psi_phi psi'_phi'`. The field is no longer separable, so
the conditioning on zero gradient at the equator and the poles is made on every level at once.

The first form, an amplitude times Gaussian correlations in shear (SPEC_19), is kept for SPEC_19's control files.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from casspian.tools.wind.estimation.state import Basis

#: The fine grid the shear covariance is formed on before projection: latitude nodes between the knots, half a
#: degree apart and offset from the equator and the poles by a quarter degree; `ln p` a quarter scale height.
FINE_LAT_STEP_DEG = 0.5
FINE_X_STEP = 0.25


def _gaussian(x, length):
    d = x[:, None] - x[None, :]
    return np.exp(-0.5 * (d / length) ** 2)


def correlation(d, length, form):
    """A correlation of the distance `d` (the same units as `length`): exponential (Matern 1/2), Matern 3/2 or 5/2,
    or Gaussian."""
    q = np.abs(np.asarray(d, dtype="float64")) / float(length)
    if form == "exponential":
        return np.exp(-q)
    if form == "matern32":
        a = np.sqrt(3.0) * q
        return (1.0 + a) * np.exp(-a)
    if form == "matern52":
        a = np.sqrt(5.0) * q
        return (1.0 + a + a ** 2 / 3.0) * np.exp(-a)
    if form == "gaussian":
        return np.exp(-0.5 * q ** 2)
    raise ValueError(f"no correlation form {form!r}; exponential, matern32, matern52 or gaussian")


def gradient_covariance(a, b, length, form):
    """The covariance in latitude of `d/dphi` of a unit field of correlation `form`, `-rho''(a - b)`, radians."""
    if form != "matern52":
        raise ValueError("the temperature's correlation must have a curvature: matern52 (SPEC_21 section 2)")
    q = np.sqrt(5.0) * np.abs(np.asarray(a)[:, None] - np.asarray(b)[None, :]) / length
    return (5.0 / (3.0 * length ** 2)) * (1.0 + q - q ** 2) * np.exp(-q)


@dataclass(frozen=True)
class Prior:
    """SPEC_19's: an amplitude times Gaussian correlations in shear between the coefficients' abscissae."""

    mean: np.ndarray
    amplitude: float
    c_lat: np.ndarray
    c_x: np.ndarray

    @classmethod
    def build(cls, basis: Basis, amplitude, lat_length_deg, x_length):
        lat, x = basis.greville()
        return cls(np.zeros(basis.size), float(amplitude), _gaussian(lat, lat_length_deg), _gaussian(x, x_length))

    def apply(self, matrix):
        """`matrix B`, for `matrix` `(m, n_lat n_x)`."""
        m = matrix.shape[0]
        stacked = matrix.reshape(m, self.c_lat.shape[0], self.c_x.shape[0])
        return (self.amplitude ** 2 * np.einsum("mij,ik,jl->mkl", stacked, self.c_lat, self.c_x,
                                                optimize=True)).reshape(m, -1)

    def covariance(self):
        return self.apply(np.eye(self.mean.size))


@dataclass(frozen=True)
class DensePrior:
    """A prior whose covariance is formed: the temperature-stated one."""

    mean: np.ndarray
    B: np.ndarray
    lengths_deg: np.ndarray = None   #: the latitude length used, `(x, lat)` on the fine grid

    def apply(self, matrix):
        return np.asarray(matrix) @ self.B

    def covariance(self):
        return self.B

    def draw(self, rng, count=1):
        """`count` states drawn from the prior, `(count, state)`."""
        jitter = 1e-10 * float(np.max(np.diag(self.B)))
        L = np.linalg.cholesky(self.B + jitter * np.eye(self.B.shape[0]))
        return self.mean + (L @ rng.standard_normal((self.B.shape[0], count))).T


def fine_grid(basis: Basis):
    """The latitudes (deg) and `ln p` of the fine grid the shear covariance is formed on."""
    lat = np.arange(-90.0 + FINE_LAT_STEP_DEG / 2, 90.0, FINE_LAT_STEP_DEG)
    x = np.linspace(basis.t_x[0], basis.t_x[-1], int(round((basis.t_x[-1] - basis.t_x[0]) / FINE_X_STEP)) + 1)
    return lat, x


def deformation_length(phi, speed, radius_m, omega):
    """The deformation radius in radians of latitude at `phi` (radians) for the speed `c` (m/s)."""
    f = 2.0 * omega * np.sin(phi)
    beta = 2.0 * omega / radius_m
    c = np.maximum(np.asarray(speed, dtype="float64"), 0.0)
    return np.where(c > 0.0, c / np.sqrt(f ** 2 + beta * c + 1e-300), 0.0) / radius_m


#: The dense grid the warp is integrated on, radians: a twentieth of a degree.
_WARP_STEP = np.radians(0.05)


def warp(phi, length, speed, radius_m, omega):
    """`(psi, psi_phi, l)` at `phi` (radians) for one level: the warped latitude, its derivative `l_0 / l`, and
    the length used, `max(l_0, L_R)`."""
    dense = np.linspace(-np.pi / 2, np.pi / 2, int(round(np.pi / _WARP_STEP)) + 1)
    used = np.maximum(length, deformation_length(dense, speed, radius_m, omega))
    rate = length / used
    psi = np.concatenate([[0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(dense))])
    psi -= np.interp(0.0, dense, psi)
    at_used = np.maximum(length, deformation_length(phi, speed, radius_m, omega))
    return np.interp(phi, dense, psi), length / at_used, at_used


def temperature_prior(basis: Basis, w, T_bar, sigma_T, lat_length_deg, lat_form, x_length, x_form,
                      speed=None, radius_m=None, omega=None, ridge=1e-8, jitter=1e-10):
    """The shear prior from the temperature covariance. `w` is `balance_weight` on the fine grid `(x, lat)` with
    `T_bar`, the anchor's temperature on the fine `x`; `sigma_T` the amplitude there (K). `speed`, `c` on the
    fine `x` (m/s), with the planet's `radius_m` and `omega`, sets the deformation-radius floor; None leaves the
    input length everywhere. The prior carries the length used on the fine grid."""
    lat, x = fine_grid(basis)
    phi = np.radians(lat)
    ends = np.radians([-90.0, 0.0, 90.0])
    points = np.concatenate([phi, ends])
    n_f, n_x = phi.size, x.size
    length = np.radians(float(lat_length_deg))
    if lat_form != "matern52":
        raise ValueError("the temperature's correlation must have a curvature: matern52 (SPEC_21 section 2)")
    speed = np.zeros(n_x) if speed is None else np.asarray(speed, dtype="float64")
    psi, rate, used = (np.stack(v) for v in zip(*(warp(points, length, speed[k], radius_m or 1.0, omega or 0.0)
                                                  for k in range(n_x))))
    k_x = correlation(x[:, None] - x[None, :], x_length, x_form)

    def block(k, l, rows, cols):
        """The gradient covariance between `points[rows]` on level `k` and `points[cols]` on level `l`."""
        return (k_x[k, l] * gradient_covariance(psi[k, rows], psi[l, cols], length, lat_form)
                * rate[k, rows][:, None] * rate[l, cols][None, :])

    fine, end = np.arange(n_f), np.arange(n_f, n_f + 3)
    # The conditioning on zero gradient at the equator and both poles, every level at once.
    K_EE = np.block([[block(k, l, end, end) for l in range(n_x)] for k in range(n_x)])
    K_EE = K_EE + jitter * float(np.max(np.diag(K_EE))) * np.eye(K_EE.shape[0])
    L_EE = np.linalg.cholesky(K_EE)
    K_FE = [np.hstack([block(k, l, fine, end) for l in range(n_x)]) for k in range(n_x)]
    scale = (np.asarray(sigma_T) / np.asarray(T_bar))[:, None] / np.asarray(w)        # (x, lat)
    a = basis.lat_matrix(lat)
    b = basis.x_matrix(np.exp(x))
    p_lat = np.linalg.solve(a.T @ a + ridge * np.eye(a.shape[1]), a.T)                 # (n_lat, lat)
    p_x = np.linalg.solve(b.T @ b + ridge * np.eye(b.shape[1]), b.T)                   # (n_x, x)
    A = [p_lat * scale[k][None, :] for k in range(n_x)]                                # (n_lat, lat) per level
    U = [np.linalg.solve(L_EE, (A[k] @ K_FE[k]).T).T for k in range(n_x)]              # (n_lat, 3 n_x)
    n_lat = p_lat.shape[0]
    M = np.empty((n_x, n_x, n_lat, n_lat))
    for k in range(n_x):
        for l in range(k, n_x):
            m = A[k] @ block(k, l, fine, fine) @ A[l].T - U[k] @ U[l].T
            M[k, l], M[l, k] = m, m.T
    B = np.einsum("klij,ak,bl->iajb", M, p_x, p_x, optimize=True).reshape(basis.size, basis.size)
    return DensePrior(np.zeros(basis.size), 0.5 * (B + B.T), np.degrees(used[:, :n_f]))
