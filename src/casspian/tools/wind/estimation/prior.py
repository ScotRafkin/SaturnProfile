"""The prior: a mean and a covariance. SPEC_19 section 1 item 4, design note v0.5 section 5.

The mean is zero shear here and carries no observation; an imposed seasonal cycle enters as the mean.
The covariance is an amplitude times Gaussian correlations in latitude and `ln p`, between the
coefficients' abscissae, so it is a Kronecker product and is applied without being formed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from casspian.tools.wind.estimation.state import Basis


def _gaussian(x, length):
    d = x[:, None] - x[None, :]
    return np.exp(-0.5 * (d / length) ** 2)


@dataclass(frozen=True)
class Prior:
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
