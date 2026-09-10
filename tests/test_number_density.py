"""Tests for number_density.py."""

import numpy as np
import pytest

from occultation_recon.constants import BOLTZMANN
from occultation_recon.io_standard import ProfileData
from occultation_recon.number_density import ensure_number_density


def _make_profile(**variables):
    """Helper to build a minimal ProfileData with given variables."""
    height = np.array([0.0, 1000.0, 2000.0])
    return ProfileData(
        height=height,
        variables=variables,
    )


def test_case_a_passthrough():
    """If number_density is already present, return it directly."""
    n_expected = np.array([1e25, 5e24, 1e24])
    profile = _make_profile(number_density=n_expected)
    n = ensure_number_density(profile)
    assert np.allclose(n, n_expected)


def test_case_b_from_p_and_t():
    """Compute n = p / (k_B * T) from pressure and temperature."""
    p = np.array([1e5, 5e4, 2e4])
    t = np.array([300.0, 250.0, 200.0])
    profile = _make_profile(pressure=p, temperature=t)
    n = ensure_number_density(profile)
    expected = p / (BOLTZMANN * t)
    assert np.allclose(n, expected)


def test_missing_variables_raises():
    """Profile with neither n, nor p+T, nor refractivity raises ValueError."""
    profile = _make_profile(temperature=np.array([300.0, 250.0, 200.0]))
    with pytest.raises(ValueError, match="Cannot compute number density"):
        ensure_number_density(profile)


def test_case_c_stub_raises():
    """Case C (refractivity) raises NotImplementedError for now."""
    profile = _make_profile(refractivity=np.array([100.0, 50.0, 10.0]))
    with pytest.raises(NotImplementedError, match="Case C"):
        ensure_number_density(profile)


def test_lindal_1bar_sanity():
    """At 1000 mbar / 134.8 K, n should be p / (k_B * T)."""
    p_1bar = 1000.0 * 100.0  # 1000 mbar in Pa
    t_1bar = 134.8            # K
    expected_n = p_1bar / (BOLTZMANN * t_1bar)

    profile = _make_profile(
        pressure=np.array([p_1bar]),
        temperature=np.array([t_1bar]),
    )
    # Adjust height to match single level
    profile.height = np.array([0.0])

    n = ensure_number_density(profile)
    assert abs(n[0] - expected_n) / expected_n < 1e-10
    # Sanity: should be approximately 5.37e25 m^-3
    assert 5.0e25 < n[0] < 6.0e25
