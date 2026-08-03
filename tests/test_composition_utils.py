"""Tests for composition_utils.py."""

import numpy as np
import pytest

from occultation_recon.constants import SPECIES, per_particle_refractivity
from occultation_recon.composition_utils import (
    mean_molar_mass,
    mean_molecular_weight_amu,
    mean_refractivity_per_particle,
    mean_molar_mass_profile,
    mean_refractivity_per_particle_profile,
)


def test_pure_h2_molar_mass():
    result = mean_molar_mass({"H2": 1.0})
    assert abs(result - SPECIES["H2"]["molar_mass_kg_per_mol"]) < 1e-12


def test_pure_h2_refractivity():
    result = mean_refractivity_per_particle({"H2": 1.0})
    assert abs(result - per_particle_refractivity("H2")) < 1e-20


def test_5050_h2_he_molar_mass():
    result = mean_molar_mass({"H2": 0.5, "He": 0.5})
    expected = 0.5 * (
        SPECIES["H2"]["molar_mass_kg_per_mol"]
        + SPECIES["He"]["molar_mass_kg_per_mol"]
    )
    assert abs(result - expected) < 1e-12


def test_5050_h2_he_refractivity():
    result = mean_refractivity_per_particle({"H2": 0.5, "He": 0.5})
    expected = 0.5 * (
        per_particle_refractivity("H2") + per_particle_refractivity("He")
    )
    assert abs(result - expected) < 1e-30


def test_bad_sum_raises():
    with pytest.raises(ValueError, match="sum"):
        mean_molar_mass({"H2": 0.5, "He": 0.3})


def test_lindal_mu():
    comp = {"H2": 0.94, "He": 0.06}
    mu = mean_molecular_weight_amu(comp)
    assert abs(mu - 2.135) < 0.01, f"Lindal mu should be ~2.135, got {mu}"


def test_profile_agrees_with_scalar():
    n_levels = 5
    comp_scalar = {"H2": 0.94, "He": 0.06}
    comp_profile = {
        "H2": np.full(n_levels, 0.94),
        "He": np.full(n_levels, 0.06),
    }
    m_scalar = mean_molar_mass(comp_scalar)
    m_profile = mean_molar_mass_profile(comp_profile)
    assert np.allclose(m_profile, m_scalar)

    r_scalar = mean_refractivity_per_particle(comp_scalar)
    r_profile = mean_refractivity_per_particle_profile(comp_profile)
    assert np.allclose(r_profile, r_scalar)
