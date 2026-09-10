"""Derived quantities from mole fractions.

Pure functions that compute mean molecular weight and mean per-particle
refractivity from a set of species mole fractions. These guarantee that
the mean molecular weight and mean refractivity always come from the same
composition, preventing inconsistencies.
"""

import numpy as np

from .constants import SPECIES, per_particle_refractivity


_SUM_TOLERANCE = 1e-6


def _validate_sum(mole_fractions: dict, label: str = "") -> None:
    """Check that mole fractions sum to 1 within tolerance."""
    total = sum(mole_fractions.values())
    if isinstance(total, np.ndarray):
        if not np.allclose(total, 1.0, atol=_SUM_TOLERANCE):
            bad = np.where(np.abs(total - 1.0) > _SUM_TOLERANCE)[0]
            raise ValueError(
                f"Mole fractions {label}do not sum to 1 at levels {bad}; "
                f"sums range from {total.min():.8f} to {total.max():.8f}"
            )
    else:
        if abs(total - 1.0) > _SUM_TOLERANCE:
            raise ValueError(
                f"Mole fractions {label}sum to {total:.8f}, not 1"
            )


def mean_molar_mass(mole_fractions: dict) -> float:
    """Composition-weighted mean molar mass.

    Parameters
    ----------
    mole_fractions : dict
        Mapping species name -> scalar mole fraction. Must sum to 1.

    Returns
    -------
    float
        Mean molar mass in kg/mol.
    """
    _validate_sum(mole_fractions)
    return sum(
        x * SPECIES[sp]["molar_mass_kg_per_mol"]
        for sp, x in mole_fractions.items()
    )


def mean_molecular_weight_amu(mole_fractions: dict) -> float:
    """Mean molecular weight in atomic mass units (dimensionless mu).

    Parameters
    ----------
    mole_fractions : dict
        Mapping species name -> scalar mole fraction. Must sum to 1.

    Returns
    -------
    float
        Mean molecular weight in amu.
    """
    from .constants import ATOMIC_MASS_UNIT, AVOGADRO
    return mean_molar_mass(mole_fractions) / (ATOMIC_MASS_UNIT * AVOGADRO)


def mean_refractivity_per_particle(mole_fractions: dict) -> float:
    """Composition-weighted per-particle refractivity.

    Returns R_mean such that total refractivity N = number_density * R_mean.

    Parameters
    ----------
    mole_fractions : dict
        Mapping species name -> scalar mole fraction. Must sum to 1.

    Returns
    -------
    float
        Mean per-particle refractivity (same units as per_particle_refractivity).
    """
    _validate_sum(mole_fractions)
    return sum(
        x * per_particle_refractivity(sp)
        for sp, x in mole_fractions.items()
    )


# ---- Vectorized (profile) variants ----

def mean_molar_mass_profile(mole_fraction_profiles: dict) -> np.ndarray:
    """Composition-weighted mean molar mass at each level.

    Parameters
    ----------
    mole_fraction_profiles : dict
        Mapping species name -> np.ndarray of mole fractions (one per level).

    Returns
    -------
    np.ndarray
        Mean molar mass in kg/mol at each level.
    """
    _validate_sum(mole_fraction_profiles, label="(profile) ")
    result = np.zeros_like(next(iter(mole_fraction_profiles.values())),
                           dtype=np.float64)
    for sp, x_arr in mole_fraction_profiles.items():
        result += x_arr * SPECIES[sp]["molar_mass_kg_per_mol"]
    return result


def mean_refractivity_per_particle_profile(
    mole_fraction_profiles: dict,
) -> np.ndarray:
    """Composition-weighted per-particle refractivity at each level.

    Parameters
    ----------
    mole_fraction_profiles : dict
        Mapping species name -> np.ndarray of mole fractions (one per level).

    Returns
    -------
    np.ndarray
        Mean per-particle refractivity at each level.
    """
    _validate_sum(mole_fraction_profiles, label="(profile) ")
    result = np.zeros_like(next(iter(mole_fraction_profiles.values())),
                           dtype=np.float64)
    for sp, x_arr in mole_fraction_profiles.items():
        result += x_arr * per_particle_refractivity(sp)
    return result
