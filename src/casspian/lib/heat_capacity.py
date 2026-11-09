"""Heat capacity of a run's mixture, and the dry adiabat. SPEC_20 v0.5, REVIEW_20_step1 ruling 1.

No run specifies a heat capacity. The mixture's molar `c_p(T)` is formed from the run's own kind C,
`c_p = sum x_i c_p,i(T)`, each species' `c_p,i(T)` from physics with its constants in
`data_static/heat_capacity.toml`:

* **H2:** `5/2 R` of translation plus its rotation, from the rotational partition function on rigid
  rotor levels `E_J = k theta J(J+1)`, para on even `J` (spin weight 1) and ortho on odd `J` (weight 3).
  The rotational heat capacity of an ensemble is the variance of its energy, `c / R = <e^2> - <e>^2` with
  `e = E / kT`. **Normal** hydrogen holds the ortho fraction frozen at 3/4, so it is `1/4 c_para + 3/4
  c_ortho`, each over its own levels. **Equilibrium** hydrogen sums over every `J` with the weights, so
  the heat of conversion is in it.
* **He:** `5/2 R`.
* **Minor species:** classical rigid rotors, `7/2 R` if linear and `4 R` if not.

The dry adiabat is `d ln T / d ln p = R / c_p(T)` with `c_p` molar; its lapse rate along the local
vertical is `g m_bar / c_p`.

Manuscript equations implemented: none. The constants are in `data_static/heat_capacity.toml`.
"""

from __future__ import annotations

import tomllib
from functools import lru_cache
from pathlib import Path

import numpy as np

from casspian.lib.constants import MOLAR_GAS_CONSTANT

DATA = Path(__file__).resolve().parents[3] / "data_static" / "heat_capacity.toml"
#: Levels summed; at 300 K the J = 20 level holds e^-480 of the population.
J_MAX = 20
MODES = ("normal", "equilibrium")


@lru_cache(maxsize=1)
def constants(path=DATA) -> dict:
    with open(path, "rb") as handle:
        return tomllib.load(handle)


def _rotational(T, J, weights, theta):
    """`c_rot / R` of one ensemble of levels."""
    T = np.atleast_1d(np.asarray(T, dtype="float64"))
    e = theta * (J * (J + 1))[None, :] / T[:, None]
    w = (weights * (2 * J + 1))[None, :] * np.exp(-(e - e.min(axis=1, keepdims=True)))
    z = w.sum(axis=1)
    mean = (w * e).sum(axis=1) / z
    return (w * e ** 2).sum(axis=1) / z - mean ** 2


def h2_rotational(T, mode: str) -> np.ndarray:
    """H2's rotational `c / R` at `T`, normal or equilibrium."""
    c = constants()["H2"]
    theta = float(c["rotational_temperature_K"])
    J = np.arange(J_MAX + 1)
    even, odd = J[J % 2 == 0], J[J % 2 == 1]
    if mode == "normal":
        f = float(c["normal_ortho_fraction"])
        return (1.0 - f) * _rotational(T, even, np.ones(even.size), theta) + f * _rotational(T, odd, np.ones(odd.size), theta)
    if mode == "equilibrium":
        weights = np.where(J % 2 == 0, float(c["spin_weight_even"]), float(c["spin_weight_odd"]))
        return _rotational(T, J, weights, theta)
    raise ValueError(f"mode {mode!r} is not one of {MODES}")


def species_cp(name: str, T, mode: str) -> np.ndarray:
    """Molar `c_p / R` of one species at `T`."""
    T = np.atleast_1d(np.asarray(T, dtype="float64"))
    if name == "H2":
        return 2.5 + h2_rotational(T, mode)
    if name == "He":
        return np.full(T.shape, 2.5)
    linear = constants()["rigid_rotor"].get(name)
    if linear is None:
        raise ValueError(f"no rigid rotor geometry for {name} in {DATA.name}")
    return np.full(T.shape, 3.5 if linear else 4.0)


def mixture_cp(composition, pressure_Pa, latitude_deg, T, mode: str) -> np.ndarray:
    """Molar `c_p` (J/mol/K) of the kind C mixture at `pressure_Pa` and one latitude, at `T`."""
    root = composition.to_dataset(inherit=False)
    names = [str(n) for n in composition["species"].to_dataset(inherit=False)["species_name"].values]
    levels = np.asarray(root["pressure_Pa"].values, dtype="float64")
    lat = np.asarray(root["latitude_planetocentric_deg"].values, dtype="float64")
    rising, north = np.argsort(np.log(levels)), np.argsort(lat)
    p = np.atleast_1d(np.asarray(pressure_Pa, dtype="float64"))
    total = np.zeros(p.shape)
    for name in names:
        x = np.asarray(root[f"x_{name}"].values, dtype="float64")[np.ix_(rising, north)]
        column = np.array([np.interp(latitude_deg, lat[north], row) for row in x])
        total = total + np.interp(np.log(p), np.log(levels[rising]), column) * species_cp(name, T, mode)
    return MOLAR_GAS_CONSTANT * total


def adiabat(T0, p0_Pa, p_top_Pa, composition, latitude_deg, mode: str, steps=400):
    """The dry adiabat from `(p0, T0)` up to `p_top`: `(pressure_Pa, temperature_K)`, RK4 in `ln p`."""
    x = np.linspace(np.log(p0_Pa), np.log(p_top_Pa), steps + 1)
    T = np.empty(x.size)
    T[0] = T0

    def rate(xi, Ti):
        return Ti * MOLAR_GAS_CONSTANT / mixture_cp(composition, np.exp(xi), latitude_deg, Ti, mode)[0]

    for i in range(steps):
        h = x[i + 1] - x[i]
        k1 = rate(x[i], T[i])
        k2 = rate(x[i] + h / 2, T[i] + h * k1 / 2)
        k3 = rate(x[i] + h / 2, T[i] + h * k2 / 2)
        k4 = rate(x[i + 1], T[i] + h * k3)
        T[i + 1] = T[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return np.exp(x), T
