"""Tests for composition.py."""

from pathlib import Path

import numpy as np
import pytest

from occultation_recon.composition import (
    constant_composition,
    lindal_composition,
    lindal_composition_with_nh3,
)
from occultation_recon.composition_utils import (
    mean_molecular_weight_amu,
    mean_molar_mass_profile,
)
from occultation_recon.io_standard import read_profile


def test_constant_composition_sum_and_length():
    comp = constant_composition({"H2": 0.94, "He": 0.06}, 10)
    assert comp["H2"].shape == (10,)
    assert comp["He"].shape == (10,)
    total = comp["H2"] + comp["He"]
    assert np.allclose(total, 1.0)


def test_constant_composition_bad_sum():
    with pytest.raises(ValueError, match="sum"):
        constant_composition({"H2": 0.5, "He": 0.3}, 5)


def test_lindal_composition_values():
    comp = lindal_composition(5)
    assert np.allclose(comp["H2"], 0.94)
    assert np.allclose(comp["He"], 0.06)


def test_lindal_composition_mu():
    comp = lindal_composition(1)
    scalar = {sp: arr[0] for sp, arr in comp.items()}
    mu = mean_molecular_weight_amu(scalar)
    assert abs(mu - 2.135) < 1e-3, f"Expected ~2.135 amu, got {mu}"


def test_lindal_with_nh3_filling():
    """Test NH3 gap-filling and renormalization against the Lindal netCDF."""
    nc_path = (
        Path(__file__).resolve().parent.parent
        / "data"
        / "lindal1985_voyager2.nc"
    )
    if not nc_path.exists():
        pytest.skip("lindal1985_voyager2.nc not built yet")

    profile = read_profile(str(nc_path))
    n_levels = len(profile.height)
    comp = lindal_composition_with_nh3(profile, n_levels)

    # All three species present
    assert "H2" in comp and "He" in comp and "NH3" in comp

    # Sum to 1 at every level
    total = comp["H2"] + comp["He"] + comp["NH3"]
    assert np.allclose(total, 1.0), f"Sum not 1: min={total.min()}, max={total.max()}"

    # No NaN in any species
    for sp, arr in comp.items():
        assert not np.any(np.isnan(arr)), f"NaN found in {sp}"

    # NH3 should be ~0 above 794 mbar (79400 Pa)
    p = profile.variables["pressure"]
    above_794 = p < 79400.0
    nh3_ppm_above = comp["NH3"][above_794] * 1e6
    # Allow a small transitional value near the boundary
    assert np.all(nh3_ppm_above < 1.0), (
        f"NH3 above 794 mbar should be ~0 ppm, max is {nh3_ppm_above.max()}"
    )

    # Interior gap at 1047.13 mbar (104713 Pa): expect ~15.9 ppm
    idx_gap = np.argmin(np.abs(p - 104713.0))
    nh3_gap_ppm = comp["NH3"][idx_gap] * 1e6
    assert abs(nh3_gap_ppm - 15.9) < 1.0, (
        f"NH3 at 1047 mbar gap: expected ~15.9 ppm, got {nh3_gap_ppm:.1f}"
    )

    # Deepest level at 1298.48 mbar (129848 Pa): expect ~79.3 ppm (extrapolated)
    idx_deep = np.argmin(np.abs(p - 129848.0))
    nh3_deep_ppm = comp["NH3"][idx_deep] * 1e6
    assert abs(nh3_deep_ppm - 79.3) < 5.0, (
        f"NH3 at 1298 mbar: expected ~79.3 ppm, got {nh3_deep_ppm:.1f}"
    )


def test_lindal_with_nh3_mu_perturbation():
    """NH3 trace species should perturb mean molecular weight by < 0.1%."""
    nc_path = (
        Path(__file__).resolve().parent.parent
        / "data"
        / "lindal1985_voyager2.nc"
    )
    if not nc_path.exists():
        pytest.skip("lindal1985_voyager2.nc not built yet")

    profile = read_profile(str(nc_path))
    n_levels = len(profile.height)

    comp_pure = lindal_composition(n_levels)
    comp_nh3 = lindal_composition_with_nh3(profile, n_levels)

    m_pure = mean_molar_mass_profile(comp_pure)
    m_nh3 = mean_molar_mass_profile(comp_nh3)

    rel_diff = np.abs(m_nh3 - m_pure) / m_pure
    assert np.all(rel_diff < 0.001), (
        f"NH3 perturbation to mean molar mass exceeds 0.1%: max={rel_diff.max():.4%}"
    )


def test_nan_in_composition_triggers_validation():
    """A composition with NaN should fail the sum-to-one check."""
    from occultation_recon.composition import _validate_fractions_sum
    bad_comp = {
        "H2": np.array([0.94, np.nan]),
        "He": np.array([0.06, 0.06]),
    }
    with pytest.raises(ValueError, match="sum"):
        _validate_fractions_sum(bad_comp)
