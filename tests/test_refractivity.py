"""Tests for refractivity.py."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from occultation_recon.constants import BOLTZMANN, LOSCHMIDT, SPECIES, per_particle_refractivity
from occultation_recon.composition import lindal_composition
from occultation_recon.composition_utils import mean_refractivity_per_particle
from occultation_recon.io_standard import read_profile
from occultation_recon.refractivity import (
    refractivity_from_number_density,
    number_density_from_refractivity,
    generate_refractivity_product,
)


def test_roundtrip():
    """Forward then inverse should recover the original number density."""
    n_levels = 10
    n_original = np.linspace(1e25, 1e23, n_levels)
    comp = lindal_composition(n_levels)

    refrac = refractivity_from_number_density(n_original, comp)
    n_recovered = number_density_from_refractivity(refrac, comp)

    assert np.allclose(n_recovered, n_original, rtol=1e-12)


def test_lindal_refractivity_plausible():
    """For the Lindal profile, refractivity should be positive and plausible."""
    nc_path = (
        Path(__file__).resolve().parent.parent
        / "data"
        / "lindal1985_voyager2.nc"
    )
    if not nc_path.exists():
        pytest.skip("lindal1985_voyager2.nc not built yet")

    profile = read_profile(str(nc_path))
    n_levels = len(profile.height)
    comp = lindal_composition(n_levels)

    p = profile.variables["pressure"]
    t = profile.variables["temperature"]
    n = p / (BOLTZMANN * t)

    refrac = refractivity_from_number_density(n, comp)

    assert np.all(refrac > 0), "Refractivity should be positive everywhere"

    # Spot check: at the 1-bar level (1000 mbar, 100000 Pa, 134.8 K)
    # n = p / (kB T) ~ 5.37e25 m^-3
    # mean per-particle refractivity for 0.94 H2 / 0.06 He:
    #   0.94 * 136/LOSCHMIDT + 0.06 * 35/LOSCHMIDT
    r_mean = mean_refractivity_per_particle({"H2": 0.94, "He": 0.06})
    idx_1bar = np.argmin(np.abs(p - 100000.0))
    n_1bar = p[idx_1bar] / (BOLTZMANN * t[idx_1bar])
    expected_N = n_1bar * r_mean
    assert abs(refrac[idx_1bar] - expected_N) / expected_N < 1e-10


def test_composition_sensitivity():
    """Different composition changes N by the ratio of mean refractivities."""
    n_levels = 5
    n_density = np.full(n_levels, 1e25)

    comp_a = {"H2": np.full(n_levels, 0.94), "He": np.full(n_levels, 0.06)}
    comp_b = {"H2": np.full(n_levels, 0.50), "He": np.full(n_levels, 0.50)}

    N_a = refractivity_from_number_density(n_density, comp_a)
    N_b = refractivity_from_number_density(n_density, comp_b)

    r_a = mean_refractivity_per_particle({"H2": 0.94, "He": 0.06})
    r_b = mean_refractivity_per_particle({"H2": 0.50, "He": 0.50})

    expected_ratio = r_a / r_b
    actual_ratio = N_a / N_b

    assert np.allclose(actual_ratio, expected_ratio, rtol=1e-12)


def test_generate_product():
    """End-to-end: generate refractivity product from Lindal netCDF."""
    nc_path = (
        Path(__file__).resolve().parent.parent
        / "data"
        / "lindal1985_voyager2.nc"
    )
    if not nc_path.exists():
        pytest.skip("lindal1985_voyager2.nc not built yet")

    profile = read_profile(str(nc_path))
    n_levels = len(profile.height)
    comp = lindal_composition(n_levels)

    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
        out_path = f.name
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        plot_path = f.name

    generate_refractivity_product(str(nc_path), comp, out_path, plot_path)

    # Read back and check
    product = read_profile(out_path)
    assert product.has("number_density")
    assert product.has("refractivity")
    assert product.has("x_H2")
    assert product.has("x_He")

    # Number density should decrease with altitude (increase with pressure)
    # Height is ordered by increasing pressure = decreasing altitude in
    # the Lindal file, so number density should increase along the array
    n = product.variables["number_density"]
    assert n[-1] > n[0], "Number density should increase with depth"

    # Refractivity should span several decades
    refrac = product.variables["refractivity"]
    assert refrac.max() / refrac.min() > 1000

    Path(out_path).unlink()
    Path(plot_path).unlink()
