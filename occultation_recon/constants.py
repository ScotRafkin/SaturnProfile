"""Physical constants and species data for occultation reconstruction.

All physical constants are SI, CODATA values.
All species properties are collected here as the single source of truth.
No constant should be hard-coded anywhere else in this package.
"""


# ---------------------------------------------------------------------------
# Physical constants (SI, CODATA)
# ---------------------------------------------------------------------------

BOLTZMANN = 1.380649e-23          # J/K, exact (2019 redefinition)
ATOMIC_MASS_UNIT = 1.66053906660e-27  # kg
AVOGADRO = 6.02214076e23          # 1/mol, exact (2019 redefinition)
LOSCHMIDT = 2.6867811e25          # m^-3, number density at 273.15 K, 101325 Pa
GAS_CONSTANT = 8.314462618        # J/(mol K)


# ---------------------------------------------------------------------------
# Species data
# ---------------------------------------------------------------------------
#
# IMPORTANT: The refractivity_stp values below are approximate placeholders
# suitable for initial development and cross-checks against Lindal et al.
# (1985). Before scientific use, replace each value with an authoritative
# radio-frequency refractivity measurement from a spectroscopic reference.
# In particular:
#   - CH4: the radio-frequency value may differ from the optical value.
#   - NH3: refractivity is frequency-dependent (polar molecule); the value
#     here is a rough estimate and should be replaced with a value at the
#     relevant frequency (2.3 GHz S-band or 8.4 GHz X-band).

SPECIES = {
    "H2": {
        "molar_mass_kg_per_mol": 2.01588e-3,
        "refractivity_stp": 136.0,
        "refractivity_reference": "Essen 1953, Orcutt & Cole 1967 (normal H2, radio freq)",
        "is_polar": False,
    },
    "He": {
        "molar_mass_kg_per_mol": 4.002602e-3,
        "refractivity_stp": 35.0,
        "refractivity_reference": "Essen 1953 (radio freq)",
        "is_polar": False,
    },
    "CH4": {
        "molar_mass_kg_per_mol": 16.0425e-3,
        "refractivity_stp": 430.0,
        # TODO: replace with authoritative radio-frequency value
        "refractivity_reference": "approximate, see Orcutt & Cole 1967; needs radio-freq verification",
        "is_polar": False,
    },
    "NH3": {
        "molar_mass_kg_per_mol": 17.0305e-3,
        "refractivity_stp": 375.0,
        # TODO: replace with frequency-specific value (polar molecule)
        "refractivity_reference": "approximate; NH3 is polar, refractivity is frequency-dependent",
        "is_polar": True,
    },
}

# Required keys that every species entry must contain
_REQUIRED_SPECIES_KEYS = {
    "molar_mass_kg_per_mol",
    "refractivity_stp",
    "refractivity_reference",
    "is_polar",
}


def per_particle_refractivity(species_name: str) -> float:
    """Return refractivity contribution per unit number density for a species.

    The per-particle refractivity alpha_i is defined such that for a pure gas
    of species i at number density n:
        N = n * alpha_i
    where N is the gas refractivity (refractive index minus 1, times 1e6).

    Computed as refractivity_stp / LOSCHMIDT, since refractivity_stp is N
    at STP where the number density equals LOSCHMIDT.

    Parameters
    ----------
    species_name : str
        Key into the SPECIES dictionary (e.g. "H2", "He").

    Returns
    -------
    float
        Per-particle refractivity in units such that N = n * alpha
        with n in m^-3 and N dimensionless (scaled by 1e6).
    """
    return SPECIES[species_name]["refractivity_stp"] / LOSCHMIDT
