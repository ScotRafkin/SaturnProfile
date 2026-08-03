"""Read and write the standard netCDF profile format.

The netCDF schema is the interchange contract between all modules.
See the build spec Section 2 for the full schema definition.
"""

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import netCDF4 as nc


@dataclass
class ProfileData:
    """Container for a profile dataset read from the standard netCDF format.

    Attributes
    ----------
    height : np.ndarray
        Geometric altitude in meters, monotonic.
    variables : dict
        Mapping variable name -> np.ndarray for all data variables
        (excluding height, which is the coordinate).
    provenance : dict
        Mapping variable name -> provenance string ("measured", "derived",
        or "assumed").
    units : dict
        Mapping variable name -> units string.
    long_names : dict
        Mapping variable name -> long_name string.
    attributes : dict
        Global attributes from the netCDF file.
    height_direction : str
        "increasing" or "decreasing", describing the monotonic direction
        of the height coordinate.
    """

    height: np.ndarray
    variables: Dict[str, np.ndarray] = field(default_factory=dict)
    provenance: Dict[str, str] = field(default_factory=dict)
    units: Dict[str, str] = field(default_factory=dict)
    long_names: Dict[str, str] = field(default_factory=dict)
    attributes: Dict[str, str] = field(default_factory=dict)
    height_direction: str = "increasing"

    def has(self, name: str) -> bool:
        """Check whether a named variable is present in the dataset."""
        return name in self.variables


def write_profile(path: str, data: dict, attributes: dict) -> None:
    """Write a profile to netCDF in the standard format.

    Parameters
    ----------
    path : str
        Output file path.
    data : dict
        Mapping variable name -> (array, units, provenance, long_name).
        Must include 'height'.
    attributes : dict
        Global attributes (title, source, height_datum, etc.).
    """
    if "height" not in data:
        raise ValueError("data must include 'height'")

    height_arr, height_units, height_prov, height_long = data["height"]
    n_levels = len(height_arr)

    Path(path).parent.mkdir(parents=True, exist_ok=True)

    with nc.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("level", n_levels)

        # Write global attributes
        for key, val in attributes.items():
            setattr(ds, key, val)

        # Write height coordinate
        h_var = ds.createVariable("height", "f8", ("level",))
        h_var[:] = height_arr
        h_var.units = height_units
        h_var.provenance = height_prov
        h_var.long_name = height_long

        # Write data variables
        for name, (arr, units, provenance, long_name) in data.items():
            if name == "height":
                continue
            v = ds.createVariable(name, "f8", ("level",))
            v[:] = arr
            v.units = units
            v.provenance = provenance
            v.long_name = long_name


def read_profile(path: str) -> ProfileData:
    """Read a standard-format netCDF profile.

    Parameters
    ----------
    path : str
        Path to the netCDF file.

    Returns
    -------
    ProfileData
        The loaded profile data.

    Raises
    ------
    ValueError
        If the file does not contain a 'height' variable.
    """
    with nc.Dataset(path, "r") as ds:
        if "height" not in ds.variables:
            raise ValueError(
                f"File {path} does not contain a 'height' variable; "
                "this is required by the standard format."
            )

        height = ds.variables["height"][:].data.copy()

        # Check monotonicity
        diffs = np.diff(height)
        if np.all(diffs > 0):
            direction = "increasing"
        elif np.all(diffs < 0):
            direction = "decreasing"
        else:
            warnings.warn(
                f"Height in {path} is not monotonic. This may cause "
                "problems in downstream processing."
            )
            direction = "non-monotonic"

        variables = {}
        provenance = {}
        units = {}
        long_names = {}

        for vname in ds.variables:
            if vname == "height":
                continue
            variables[vname] = ds.variables[vname][:].data.copy()

            if hasattr(ds.variables[vname], "provenance"):
                provenance[vname] = ds.variables[vname].provenance
            else:
                warnings.warn(
                    f"Variable '{vname}' in {path} lacks a 'provenance' "
                    "attribute."
                )

            if hasattr(ds.variables[vname], "units"):
                units[vname] = ds.variables[vname].units
            else:
                warnings.warn(
                    f"Variable '{vname}' in {path} lacks a 'units' attribute."
                )

            if hasattr(ds.variables[vname], "long_name"):
                long_names[vname] = ds.variables[vname].long_name

        # Also capture height metadata
        if hasattr(ds.variables["height"], "provenance"):
            provenance["height"] = ds.variables["height"].provenance
        if hasattr(ds.variables["height"], "units"):
            units["height"] = ds.variables["height"].units
        if hasattr(ds.variables["height"], "long_name"):
            long_names["height"] = ds.variables["height"].long_name

        # Global attributes
        attrs = {}
        for attr in ds.ncattrs():
            attrs[attr] = getattr(ds, attr)

    return ProfileData(
        height=height,
        variables=variables,
        provenance=provenance,
        units=units,
        long_names=long_names,
        attributes=attrs,
        height_direction=direction,
    )
