"""Tests for io_standard.py and the Table I netCDF output."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from occultation_recon.io_standard import write_profile, read_profile


def test_roundtrip():
    """Write a small synthetic profile, read it back, check arrays and attrs."""
    height = np.array([0.0, 1000.0, 2000.0])
    pressure = np.array([1e5, 5e4, 2e4])

    data = {
        "height": (height, "m", "measured", "Geometric altitude"),
        "pressure": (pressure, "Pa", "derived", "Atmospheric pressure"),
    }
    attrs = {"title": "Test profile", "source": "unit test"}

    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
        tmp_path = f.name

    write_profile(tmp_path, data, attrs)
    profile = read_profile(tmp_path)

    assert np.allclose(profile.height, height)
    assert np.allclose(profile.variables["pressure"], pressure)
    assert profile.provenance["pressure"] == "derived"
    assert profile.units["pressure"] == "Pa"
    assert profile.attributes["title"] == "Test profile"
    assert profile.height_direction == "increasing"

    Path(tmp_path).unlink()


def test_missing_height_raises():
    """Reading a file without height should raise a clear error."""
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
        tmp_path = f.name

    import netCDF4 as nc
    with nc.Dataset(tmp_path, "w") as ds:
        ds.createDimension("level", 3)
        v = ds.createVariable("pressure", "f8", ("level",))
        v[:] = [1e5, 5e4, 2e4]

    with pytest.raises(ValueError, match="height"):
        read_profile(tmp_path)

    Path(tmp_path).unlink()


def test_has_method():
    height = np.array([0.0, 1000.0])
    data = {
        "height": (height, "m", "measured", "Altitude"),
        "temperature": (np.array([300.0, 250.0]), "K", "derived", "Temp"),
    }
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
        tmp_path = f.name

    write_profile(tmp_path, data, {"title": "test"})
    profile = read_profile(tmp_path)

    assert profile.has("temperature")
    assert not profile.has("pressure")

    Path(tmp_path).unlink()


def test_table1_netcdf():
    """Verify the Lindal Table I netCDF after it has been built."""
    nc_path = (
        Path(__file__).resolve().parent.parent
        / "data"
        / "lindal1985_voyager2.nc"
    )
    if not nc_path.exists():
        pytest.skip("lindal1985_voyager2.nc not yet built; run build_table1_netcdf first")

    profile = read_profile(str(nc_path))

    # Pressure should be monotonic increasing
    assert np.all(np.diff(profile.variables["pressure"]) > 0), \
        "Pressure not monotonic increasing"

    # Find the 1000 mbar = 100000 Pa level
    p = profile.variables["pressure"]
    idx_1bar = np.argmin(np.abs(p - 100000.0))

    # Altitude at 1000 mbar should be 0.0
    assert abs(profile.height[idx_1bar]) < 1.0, \
        f"Altitude at 1000 mbar should be ~0, got {profile.height[idx_1bar]}"

    # Temperature at 1000 mbar should be 134.8 K
    t = profile.variables["temperature"]
    assert abs(t[idx_1bar] - 134.8) < 0.1, \
        f"T at 1000 mbar should be 134.8 K, got {t[idx_1bar]}"
