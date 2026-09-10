"""Ensure total number density n(h) is available from a profile.

Dispatches on which variables the profile contains to compute number density
from the available data using local algebraic relations. No integration is
performed in this module.
"""

import numpy as np

from .constants import BOLTZMANN
from .io_standard import ProfileData


def ensure_number_density(profile: ProfileData, composition=None) -> np.ndarray:
    """Return total number density in m^-3 at each level of the profile.

    Dispatches on the variables present in the profile:

    Case A: profile has 'number_density' already -- return it directly.
    Case B: profile has 'pressure' and 'temperature' -- compute
            n = pressure / (k_B * temperature). This is the exact algebraic
            identity for an ideal gas and carries no composition assumption.
            This is the Lindal case.
    Case C: profile has 'refractivity' and composition is provided -- compute
            n = refractivity / mean_refractivity_per_particle(composition).
            Not yet implemented.

    Parameters
    ----------
    profile : ProfileData
        A standard-format profile dataset.
    composition : dict, optional
        Mole fraction profiles (species -> array), required only for Case C.

    Returns
    -------
    np.ndarray
        Total number density in m^-3 at each level.

    Raises
    ------
    ValueError
        If the profile does not contain sufficient variables for any case.
    NotImplementedError
        If Case C is triggered (refractivity + composition).
    """
    # Case A: number_density already present
    if profile.has("number_density"):
        return profile.variables["number_density"].copy()

    # Case B: compute from pressure and temperature (ideal gas law)
    if profile.has("pressure") and profile.has("temperature"):
        p = profile.variables["pressure"]
        t = profile.variables["temperature"]
        return p / (BOLTZMANN * t)

    # Case C: from refractivity and composition (stub)
    if profile.has("refractivity"):
        raise NotImplementedError(
            "Case C (number density from refractivity and composition) is not "
            "yet implemented. This will be added when refractivity-as-input is "
            "needed. Requires a composition argument providing mole fraction "
            "profiles to compute mean per-particle refractivity."
        )

    raise ValueError(
        "Cannot compute number density: profile must contain either "
        "'number_density', or 'pressure' and 'temperature', or "
        "'refractivity' (with a composition argument). "
        f"Available variables: {list(profile.variables.keys())}"
    )
