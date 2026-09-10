"""Refractivity from number density and composition.

The pure physical relation between number density, composition, and
refractivity. This module produces the fundamental data product:
a refractivity profile from an occultation-derived number density profile
and an assumed atmospheric composition.
"""

from pathlib import Path
from typing import Optional

import numpy as np

from .composition_utils import mean_refractivity_per_particle_profile
from .io_standard import ProfileData, read_profile, write_profile
from .number_density import ensure_number_density


def refractivity_from_number_density(
    number_density: np.ndarray, composition: dict
) -> np.ndarray:
    """Compute refractivity from number density and composition.

    N(h) = number_density(h) * mean_refractivity_per_particle(composition, h)

    Pure physical relation. Agnostic about the source of its inputs.
    Feed Lindal's n(h) and Lindal's composition to recover measured N(h).
    Feed any synthetic n(h) and any composition for that atmosphere's N(h).

    Parameters
    ----------
    number_density : np.ndarray
        Total number density in m^-3 at each level.
    composition : dict
        Mapping species name -> np.ndarray of mole fractions at each level.

    Returns
    -------
    np.ndarray
        Refractivity N = (n_refractive - 1) * 1e6, dimensionless (scaled).
    """
    r_mean = mean_refractivity_per_particle_profile(composition)
    return number_density * r_mean


def number_density_from_refractivity(
    refractivity: np.ndarray, composition: dict
) -> np.ndarray:
    """Compute number density from refractivity and composition.

    n(h) = refractivity(h) / mean_refractivity_per_particle(composition, h)

    The inverse relation. Used to reinterpret a fixed measured refractivity
    under an alternative composition.

    Parameters
    ----------
    refractivity : np.ndarray
        Refractivity N = (n_refractive - 1) * 1e6 at each level.
    composition : dict
        Mapping species name -> np.ndarray of mole fractions at each level.

    Returns
    -------
    np.ndarray
        Total number density in m^-3 at each level.
    """
    r_mean = mean_refractivity_per_particle_profile(composition)
    return refractivity / r_mean


def generate_refractivity_product(
    profile_path: str,
    composition: dict,
    output_path: str,
    plot_path: Optional[str] = None,
) -> None:
    """High-level driver for the fundamental refractivity data product.

    1. Read the profile from the standard netCDF.
    2. Ensure number density is available.
    3. Compute refractivity from number density and composition.
    4. Write an output netCDF with height, number_density, refractivity,
       and the composition mole fractions, with full provenance.
    5. Optionally produce a diagnostic figure.

    Parameters
    ----------
    profile_path : str
        Path to the input standard-format netCDF.
    composition : dict
        Mapping species name -> np.ndarray of mole fractions at each level.
        Must match the number of levels in the input profile.
    output_path : str
        Path for the output netCDF.
    plot_path : str, optional
        If given, save a diagnostic figure to this path.
    """
    profile = read_profile(profile_path)
    n = ensure_number_density(profile)
    refractivity = refractivity_from_number_density(n, composition)

    # Build composition note from species names and representative values
    species_parts = []
    for sp, x_arr in composition.items():
        x_mid = x_arr[len(x_arr) // 2]
        if x_mid > 0.01:
            species_parts.append(f"{sp}: {x_mid * 100:.1f}%")
        else:
            species_parts.append(f"{sp}: {x_mid * 1e6:.1f} ppm")
    comp_note = "Composition: " + ", ".join(species_parts)

    # Build output data
    data = {
        "height": (
            profile.height,
            "m",
            profile.provenance.get("height", "measured"),
            "Geometric altitude relative to 1-bar pressure level",
        ),
        "number_density": (
            n,
            "m-3",
            "derived",
            "Total number density from ideal gas law (p / kB T)",
        ),
        "refractivity": (
            refractivity,
            "dimensionless",
            "derived",
            "Gas refractivity N = (n_refractive - 1) * 1e6",
        ),
    }

    # Include composition mole fractions
    for sp, x_arr in composition.items():
        data[f"x_{sp}"] = (
            x_arr,
            "mol/mol",
            "assumed",
            f"Mole fraction of {sp}",
        )

    # Carry forward source attributes, update with product info
    attrs = dict(profile.attributes)
    attrs["title"] = attrs.get("title", "Refractivity product") + " -- refractivity product"
    attrs["composition_note"] = comp_note
    attrs["history"] = (
        attrs.get("history", "")
        + "; refractivity computed by generate_refractivity_product"
    )

    write_profile(output_path, data, attrs)

    if plot_path is not None:
        _make_plot(profile.height, n, refractivity, comp_note, plot_path)


def _make_plot(
    height: np.ndarray,
    number_density: np.ndarray,
    refractivity: np.ndarray,
    comp_note: str,
    plot_path: str,
) -> None:
    """Produce a two-panel diagnostic figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    height_km = height / 1000.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7), sharey=True)

    ax1.semilogx(refractivity, height_km)
    ax1.set_xlabel("Refractivity N")
    ax1.set_ylabel("Altitude (km rel. 1-bar)")
    ax1.set_title("Refractivity")
    ax1.grid(True, alpha=0.3)

    ax2.semilogx(number_density, height_km)
    ax2.set_xlabel("Number density (m$^{-3}$)")
    ax2.set_title("Number density")
    ax2.grid(True, alpha=0.3)

    fig.suptitle(comp_note, fontsize=10)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
