"""Composition profiles as mole fractions versus height.

Multiple constructors produce the same representation: a dict mapping
species name to a numpy array of mole fractions (one per level), summing
to 1 at every level. This representation is directly consumable by
composition_utils functions.
"""

import numpy as np
from scipy.interpolate import interp1d

from .io_standard import ProfileData


_SUM_TOLERANCE = 1e-6


def _validate_fractions_sum(comp: dict, label: str = "") -> None:
    """Check that mole fractions sum to 1 at every level."""
    total = sum(comp.values())
    if isinstance(total, np.ndarray):
        if not np.allclose(total, 1.0, atol=_SUM_TOLERANCE):
            bad = np.where(np.abs(total - 1.0) > _SUM_TOLERANCE)[0]
            raise ValueError(
                f"Mole fractions {label}do not sum to 1 at levels {bad}"
            )
    else:
        if abs(total - 1.0) > _SUM_TOLERANCE:
            raise ValueError(
                f"Mole fractions {label}sum to {total}, not 1"
            )


def constant_composition(species_fractions: dict, n_levels: int) -> dict:
    """Constant mole fractions at every level.

    Parameters
    ----------
    species_fractions : dict
        Mapping species name -> scalar mole fraction. Must sum to 1.
    n_levels : int
        Number of vertical levels.

    Returns
    -------
    dict
        Mapping species name -> np.ndarray of length n_levels.
    """
    total = sum(species_fractions.values())
    if abs(total - 1.0) > _SUM_TOLERANCE:
        raise ValueError(
            f"species_fractions sum to {total}, not 1"
        )
    return {sp: np.full(n_levels, x) for sp, x in species_fractions.items()}


def lindal_composition(n_levels: int) -> dict:
    """Lindal et al. 1985 assumed composition: 0.94 H2, 0.06 He by number.

    Constant with height. Reproduces mean molecular weight 2.135 amu.
    This is the composition that recovers the measured refractivity when
    combined with the recovered number density.

    Parameters
    ----------
    n_levels : int
        Number of vertical levels.

    Returns
    -------
    dict
        Mapping species name -> np.ndarray of length n_levels.
    """
    return constant_composition({"H2": 0.94, "He": 0.06}, n_levels)


def lindal_composition_with_nh3(profile: ProfileData, n_levels: int) -> dict:
    """Lindal H2/He composition with NH3 included as a trace species.

    Reads nh3_mixing_ratio from the profile, fills NaN values by linear
    interpolation and extrapolation in pressure, clamps to non-negative,
    then renormalizes H2 and He so all three species sum to 1 at every level.
    The H2:He number ratio is held at 0.94:0.06.

    Parameters
    ----------
    profile : ProfileData
        A standard-format profile containing 'pressure' and 'nh3_mixing_ratio'.
    n_levels : int
        Number of vertical levels (must match profile dimensions).

    Returns
    -------
    dict
        Mapping species name -> np.ndarray of length n_levels.
    """
    if not profile.has("nh3_mixing_ratio"):
        raise ValueError("Profile must contain 'nh3_mixing_ratio'")
    if not profile.has("pressure"):
        raise ValueError("Profile must contain 'pressure' for NH3 interpolation")

    pressure = profile.variables["pressure"]
    nh3_raw = profile.variables["nh3_mixing_ratio"].copy()

    # Fill NaN values by linear interpolation/extrapolation in pressure.
    # Use scipy interp1d with fill_value="extrapolate" so that both ends
    # are linearly extrapolated (not clamped to the endpoint value).
    valid = ~np.isnan(nh3_raw)
    if valid.sum() == 0:
        # No NH3 data at all: set to zero everywhere
        nh3_filled = np.zeros(n_levels)
    elif valid.sum() == 1:
        # Single point: cannot interpolate, use that value everywhere it
        # would be positive, zero elsewhere
        nh3_filled = np.full(n_levels, nh3_raw[valid][0])
    else:
        interp_func = interp1d(
            pressure[valid],
            nh3_raw[valid],
            kind="linear",
            fill_value="extrapolate",
            bounds_error=False,
        )
        nh3_filled = interp_func(pressure)

    # Clamp to non-negative. This is load-bearing: without it, the upper
    # atmosphere receives large negative mixing ratios from extrapolation
    # of the rising deep-troposphere trend, which are unphysical.
    nh3_filled = np.maximum(nh3_filled, 0.0)

    # NOTE: The deepest-level NH3 value is an unconstrained linear
    # extrapolation of a trace species, adequate because NH3 has negligible
    # effect on mean molecular weight and refractivity, but not a measured
    # or physically modeled value and should not be over-interpreted.
    # Setting NH3 to zero above the cloud is an approximation to the true
    # sub-saturation vapor profile, adequate here because the vapor is
    # negligible at those levels.

    # Renormalize: assign trace species first, distribute remainder among
    # major species in their fixed ratio (0.94 H2 : 0.06 He).
    r_h2 = 0.94
    r_he = 0.06
    r_total = r_h2 + r_he

    remaining = 1.0 - nh3_filled
    x_h2 = remaining * r_h2 / r_total
    x_he = remaining * r_he / r_total

    comp = {"H2": x_h2, "He": x_he, "NH3": nh3_filled}
    _validate_fractions_sum(comp, label="(lindal_with_nh3) ")
    return comp


def analytical_composition(height: np.ndarray, spec) -> dict:
    """Composition from an analytical prescription.

    Stub for future implementation. The 'spec' argument would define a
    functional form (e.g. constant below a knee, linear ramp above).

    Parameters
    ----------
    height : np.ndarray
        Altitude grid in meters.
    spec : object
        Specification of the functional form (to be defined).

    Returns
    -------
    dict
        Mapping species name -> np.ndarray.
    """
    raise NotImplementedError(
        "analytical_composition is a stub. Implement when a specific "
        "analytical prescription is needed."
    )


def composition_from_file(path: str, height: np.ndarray) -> dict:
    """Read a composition profile from an external file.

    Stub for future implementation. Intended for external model output
    (e.g. Julie Moses photochemical model) interpolated onto the given
    height grid.

    Expected file format: CSV or similar with columns for altitude and
    species mole fractions. To be defined when a real file is in hand.

    Parameters
    ----------
    path : str
        Path to the composition file.
    height : np.ndarray
        Target altitude grid in meters.

    Returns
    -------
    dict
        Mapping species name -> np.ndarray.
    """
    raise NotImplementedError(
        "composition_from_file is a stub. Implement when an external "
        "composition file (e.g. Moses photochemical model) is available."
    )
