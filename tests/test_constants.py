"""Tests for constants.py."""

from occultation_recon.constants import (
    LOSCHMIDT,
    SPECIES,
    _REQUIRED_SPECIES_KEYS,
    per_particle_refractivity,
)


def test_all_species_have_required_keys():
    for name, props in SPECIES.items():
        missing = _REQUIRED_SPECIES_KEYS - set(props.keys())
        assert not missing, f"Species {name} missing keys: {missing}"


def test_molar_masses_positive_and_plausible():
    for name, props in SPECIES.items():
        m = props["molar_mass_kg_per_mol"]
        assert m > 0, f"{name} molar mass not positive"
        # All species here are lighter than 100 g/mol
        assert m < 0.1, f"{name} molar mass implausibly large: {m}"


def test_per_particle_refractivity_roundtrip():
    """per_particle_refractivity(H2) * LOSCHMIDT should recover refractivity_stp."""
    alpha = per_particle_refractivity("H2")
    recovered = alpha * LOSCHMIDT
    assert abs(recovered - SPECIES["H2"]["refractivity_stp"]) < 1e-10
