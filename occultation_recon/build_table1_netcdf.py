"""One-time digitization of Lindal et al. 1985 Table I into standard netCDF.

Reads the CSV transcription of Table I (Voyager 2 ingress), converts units,
and writes a standard-format netCDF file.

Usage:
    python -m occultation_recon.build_table1_netcdf
"""

import csv
import math
from pathlib import Path

import numpy as np

from .io_standard import write_profile


_CSV_PATH = Path(__file__).resolve().parent.parent / "data" / "lindal1985_table1.csv"
_NC_PATH = Path(__file__).resolve().parent.parent / "data" / "lindal1985_voyager2.nc"


def build_netcdf(csv_path: str = None, nc_path: str = None) -> None:
    """Read the Table I CSV and write the standard netCDF.

    Parameters
    ----------
    csv_path : str, optional
        Path to the input CSV. Defaults to data/lindal1985_table1.csv.
    nc_path : str, optional
        Path to the output netCDF. Defaults to data/lindal1985_voyager2.nc.
    """
    if csv_path is None:
        csv_path = str(_CSV_PATH)
    if nc_path is None:
        nc_path = str(_NC_PATH)

    # Read CSV
    pressure_mbar = []
    temperature_K = []
    nh3_ppm = []
    altitude_km = []

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pressure_mbar.append(float(row["pressure_mbar"]))
            temperature_K.append(float(row["temperature_K"]))
            nh3_val = row["nh3_ppm"].strip()
            nh3_ppm.append(float(nh3_val) if nh3_val else float("nan"))
            altitude_km.append(float(row["altitude_km"]))

    # Convert to numpy arrays
    pressure_mbar = np.array(pressure_mbar)
    temperature_K = np.array(temperature_K)
    nh3_ppm = np.array(nh3_ppm)
    altitude_km = np.array(altitude_km)

    # Convert units
    pressure_pa = pressure_mbar * 100.0       # mbar -> Pa
    altitude_m = altitude_km * 1000.0          # km -> m
    nh3_mole_fraction = nh3_ppm * 1.0e-6       # ppm -> mole fraction

    # Order by increasing pressure (decreasing altitude)
    sort_idx = np.argsort(pressure_pa)
    pressure_pa = pressure_pa[sort_idx]
    temperature_K = temperature_K[sort_idx]
    nh3_mole_fraction = nh3_mole_fraction[sort_idx]
    altitude_m = altitude_m[sort_idx]

    # Build data dict for write_profile
    data = {
        "height": (
            altitude_m,
            "m",
            "measured",
            "Geometric altitude relative to 1-bar pressure level",
        ),
        "pressure": (
            pressure_pa,
            "Pa",
            "derived",
            "Atmospheric pressure from Lindal hydrostatic integration",
        ),
        "temperature": (
            temperature_K,
            "K",
            "derived",
            "Gas temperature from number density and pressure via ideal gas law",
        ),
        "nh3_mixing_ratio": (
            nh3_mole_fraction,
            "mol/mol",
            "derived",
            "Ammonia mole fraction from microwave absorption (NaN where not measured)",
        ),
    }

    attributes = {
        "title": "Voyager 2 Radio Occultation Ingress Profile, Saturn",
        "source": "Lindal, Sweetnam, Eshleman 1985, AJ 90, 1136, Table I",
        "height_datum": "1-bar pressure level as defined in source",
        "latitude_deg": (
            "36.5 N planetographic (ingress, evening terminator); "
            "see source Section on Voyager 2"
        ),
        "longitude_deg": "186.8 E system III (start of ingress swath)",
        "history": "Created from digitized Table I by build_table1_netcdf.py",
        "composition_note": (
            "Source p and T derived assuming 94% H2 / 6% He by number, "
            "mean molecular mass 2.135 amu, per Lindal et al. 1985"
        ),
    }

    write_profile(nc_path, data, attributes)
    print(f"Wrote {nc_path} ({len(pressure_pa)} levels)")


if __name__ == "__main__":
    build_netcdf()
