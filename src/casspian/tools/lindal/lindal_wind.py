"""The Lindal-only wind: the cloud wind with the shear of the IRIS temperatures. SPEC_11 Step 1.

Inputs: the cloud wind (kind W, as `casspian-wind-from-curve` writes it), the digitized IRIS
temperatures (`occul_data/lindal/iris_temperatures.csv`, SPEC_09), and the run's composition,
gravity and rotation. Output: a kind W dataset with all three parts. `construct` touches no file;
a run's build writes the wind through the `lindal_iris` case of `casspian-wind-shear`, which calls
it (SPEC_11 v0.4 Step 2).

**The reference wind** is the cloud wind's reference wind, assigned to `REFERENCE_PRESSURE_Pa`
(398 mbar, a node of the output grid), where the shear is zero.

**The gradient** at each of the three levels is `iris_temperatures.fit` at its defaults, evaluated
at the planetographic latitudes of the cloud wind's planetocentric grid (the file's own
`latitude_planetographic_deg`, which the wind tool computed by its rule), and converted to per
planetocentric radian with `dphi_g/dphi_c` differenced on that grid.

**The shear from the model's balance.** Along an isobar the transfer integrates
`d ln N / dphi = K = S / g` (Eq. A27 with the composition uniform in latitude, which the run's is),
and at fixed pressure `ln N = const - ln T`, so the relation the transfer uses is

    (d ln T / dphi)_p = -S / g,    S = 2 Omega_abs r (du/dZ)_R     (Eq. A15)

with `(du/dZ)_R = sin(phi) (du/dr)_phi + (cos(phi) / r) (du/dphi)_r` and the derivatives taken to
fixed pressure as `lib.kernel` takes them. Solved for the shear `s = (du/dln p)_phi`:

    s = [ -g (d ln T/dphi)_p / (2 Omega_abs r) - (cos(phi) / r) (du/dphi)_p ] / D
    D = (dln p/dr) [ sin(phi) - (cos(phi) / r) (dr/dphi)_p ],   dln p/dr = -g m_bar / (R T)

`D` vanishes like `sin(phi)` at the equator, the singularity of SPEC_11 section 2 ruling 2. At each
level the geometry is the model's own: `r(phi)` the Eq. B3 surface under that level's wind
(`lib.geoid.wind_geoid`, anchored at `equatorial_radius_m`, by default the gravity file's
normalization radius), `g` the effective radial gravity and `Omega_abs` the absolute rate under that
wind (`lib.gravity`), `T` the IRIS fit's value, and `m_bar` the run's composition at the level.

**In the vertical**, between the three levels linear in `ln p`; above the top level held at its
value or relaxed as `s_top exp(-x / ABOVE_SCALE_HEIGHTS)`, `x = ln(p_top / p)`; below the bottom
level relaxed as `s_bottom exp(-y / BELOW_SCALE_HEIGHTS)`, `y = ln(p / p_bottom)`. The wind is the
reference wind plus the shear integrated in `ln p` from the reference level, in closed form.

**In latitude**, the shear is computed where the fit is defined (between the level's northernmost
and southernmost dots) and poleward of `EQUATORIAL_BAND_DEG` planetocentric. Across the band it is
PCHIP through the computed values, which stays between the two edge values; beyond the data's ends
it tapers linearly to zero at each pole.

**The iteration.** The shear changes the wind, and the wind enters the balance (`Omega_abs`, `g`,
the surface, `(du/dphi)_p`), so the shear, its shaping and the integration are repeated from the
cloud wind until the wind changes nowhere by more than `TOLERANCE_MS`. Not converging in
`MAX_ITERATIONS` is refused, named.

The uncertainty is the cloud wind's, carried unchanged at every pressure. `value_provenance` is the
cloud wind's at the reference level and `parameterized` elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import xarray as xr

from casspian.lib import geoid as gd
from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.lib.control import ControlFileError
from casspian.lib.gravity import g_eff_radial, omega_abs
from casspian.lib.reduction import mean_over_species
from casspian.tools.lindal import iris_temperatures as iris

REFERENCE_PRESSURE_Pa = 39810.7
#: The printed levels of the IRIS file (mbar) and the pressures they are placed at (Pa); the top
#: level's placement is the case's (`top_level_Pa`).
LOWER_LEVELS_Pa = {290: 29000.0, 730: 73000.0}
TOP_LEVEL_MBAR = 110
NODES_PER_DECADE = 20
GRID_DECADES = (0, 6)
ABOVE_CASES = ("held", "relaxed")
ABOVE_SCALE_HEIGHTS = 2.0
BELOW_SCALE_HEIGHTS = 1.0
EQUATORIAL_BAND_DEG = 5.0
TOLERANCE_MS = 0.1
MAX_ITERATIONS = 50
PARAMETERIZED = 2


@dataclass
class Construction:
    """The wind and what went into it, for the acceptance figures."""

    dataset: xr.Dataset
    level_Pa: np.ndarray                    #: the three levels, top first
    shear_ms: np.ndarray                    #: (level, latitude), du/dln p after shaping
    shear_computed_ms: np.ndarray           #: (level, latitude), the balance's, NaN where not computed
    temperature_K: np.ndarray               #: (level, latitude), the fit's value
    gradient_K_per_rad: np.ndarray          #: (level, latitude), dT/dphi_c
    iterations: int
    changes_ms: list = field(default_factory=list)


def pressure_grid() -> np.ndarray:
    """1 Pa to 1 MPa at `NODES_PER_DECADE`, each node rounded to six significant digits as the
    build files write them, so that 39810.7 Pa is a node exactly."""
    low, high = GRID_DECADES
    k = np.arange(low * NODES_PER_DECADE, high * NODES_PER_DECADE + 1)
    return np.array([float(f"{10.0 ** (i / NODES_PER_DECADE):.6g}") for i in k])


def mean_molar_mass(composition, pressure_Pa, latitude_deg) -> np.ndarray:
    """`m_bar` (kg/mol) of the run's composition at each pressure and latitude, log-linear in
    pressure and linear in latitude on the file's own grid. Shape `(pressure, latitude)`."""
    root = composition.to_dataset(inherit=False)
    species = composition["species"].to_dataset(inherit=False)
    names = [str(name) for name in species["species_name"].values]
    x = np.stack([np.asarray(root[f"x_{name}"].values, dtype="float64") for name in names])
    m_bar = mean_over_species(x, np.asarray(species["molar_mass_kg_mol"].values, dtype="float64"))
    levels = np.asarray(root["pressure_Pa"].values, dtype="float64")
    lat = np.asarray(root["latitude_planetocentric_deg"].values, dtype="float64")
    rising, north = np.argsort(np.log(levels)), np.argsort(lat)
    m_bar = m_bar[np.ix_(rising, north)]
    on_levels = np.stack([np.interp(np.log(pressure_Pa), np.log(levels[rising]), m_bar[:, j])
                          for j in range(lat.size)], axis=1)
    return np.stack([np.interp(latitude_deg, lat[north], row) for row in on_levels])


def integrated_shear(shear, levels_Pa, pressure_Pa, above):
    """`int s dln p` from the reference level, closed form, shape `(latitude, pressure)`.

    `shear` is `(3, latitude)` at `levels_Pa` (top first). Linear in `ln p` between the levels,
    held or relaxed above the top, relaxed below the bottom (module docstring)."""
    x1, x2, x3 = np.log(levels_Pa)
    s1, s2, s3 = (shear[i][:, None] for i in range(3))

    def antiderivative(x):
        x = np.asarray(x, dtype="float64")[None, :]
        f2 = 0.5 * (s1 + s2) * (x2 - x1)
        f3 = f2 + 0.5 * (s2 + s3) * (x3 - x2)
        if above == "held":
            top = s1 * (x - x1)
        else:
            top = -s1 * ABOVE_SCALE_HEIGHTS * (1.0 - np.exp(-(x1 - x) / ABOVE_SCALE_HEIGHTS))
        upper = s1 * (x - x1) + (s2 - s1) * (x - x1) ** 2 / (2.0 * (x2 - x1))
        lower = f2 + s2 * (x - x2) + (s3 - s2) * (x - x2) ** 2 / (2.0 * (x3 - x2))
        bottom = f3 + s3 * BELOW_SCALE_HEIGHTS * (1.0 - np.exp(-(x - x3) / BELOW_SCALE_HEIGHTS))
        return np.where(x < x1, top, np.where(x <= x2, upper,
                                              np.where(x <= x3, lower, bottom)))

    return antiderivative(np.log(pressure_Pa)) - antiderivative([np.log(REFERENCE_PRESSURE_Pa)])


def shape_in_latitude(computed, latitude_deg) -> np.ndarray:
    """The shear on the whole grid from its computed values (NaN elsewhere): the equatorial band
    bridged by PCHIP, the polar ends tapered linearly to zero."""
    from scipy.interpolate import PchipInterpolator

    s = np.array(computed, dtype="float64")
    lat = np.asarray(latitude_deg, dtype="float64")
    known = np.isfinite(s)
    band = np.abs(lat) < EQUATORIAL_BAND_DEG
    south, north = np.flatnonzero(known & (lat < 0)), np.flatnonzero(known & (lat > 0))
    if south.size < 2 or north.size < 2:
        raise ValueError("the shear is computed on fewer than two latitudes on one side of the "
                         "equatorial band; nothing to bridge from")
    i0, i1 = south.max(), north.min()
    if not (np.all(~known[i0 + 1:i1]) and band[i0 + 1:i1].all()):
        raise ValueError("the computed shear is not two pieces either side of the equatorial band")
    # SPEC_11 v0.4 section 2 item 4: PCHIP through the computed values, which across the band is
    # monotone between the two edge values and so never leaves them.
    s[i0 + 1:i1] = PchipInterpolator(lat[known], s[known])(lat[i0 + 1:i1])
    first, last = south.min(), north.max()
    s[:first] = s[first] * (lat[:first] + 90.0) / (lat[first] + 90.0)
    s[last + 1:] = s[last] * (90.0 - lat[last + 1:]) / (90.0 - lat[last])
    return s


def construct(cloud, table, composition, gravity, rotation, top_level_Pa, above,
              equatorial_radius_m=None) -> Construction:
    """The wind for one case, from in-memory inputs. No file is read or written.

    `table` is the IRIS CSV as a structured array (`numpy.genfromtxt(..., names=True)`)."""
    if above not in ABOVE_CASES:
        raise ControlFileError(f"above_top = {above!r}; the cases are {list(ABOVE_CASES)}")
    lat_c = np.asarray(cloud["latitude_planetocentric_deg"].values, dtype="float64")
    lat_g = np.asarray(cloud["latitude_planetographic_deg"].values, dtype="float64")
    phi = np.radians(lat_c)
    u_ref = np.asarray(cloud["u_reference_ms"].values, dtype="float64")
    dg_dc = np.gradient(lat_g, lat_c)

    Omega = float(rotation["angular_rate_rad_s"])
    GM = float(gravity["GM_m3s2"])
    J = np.asarray(gravity["J"].values, dtype="float64")
    degrees = np.asarray(gravity["degree"].values)
    R_norm = float(gravity["normalization_radius_m"])
    r_eq = R_norm if equatorial_radius_m is None else float(equatorial_radius_m)

    levels = np.array([float(top_level_Pa), LOWER_LEVELS_Pa[290], LOWER_LEVELS_Pa[730]])
    printed = (TOP_LEVEL_MBAR, 290, 730)
    temperature = np.full((3, lat_c.size), np.nan)
    gradient = np.full((3, lat_c.size), np.nan)
    for k, level in enumerate(printed):
        rows = table["pressure_mbar"] == level
        result = iris.fit(table["latitude_planetographic_deg"][rows], table["temperature_K"][rows],
                          grid=lat_g)
        temperature[k] = result.value_K
        gradient[k] = result.gradient_K_per_deg * dg_dc * (180.0 / np.pi)
    m_bar = mean_molar_mass(composition, levels, lat_c)
    computed_mask = np.isfinite(gradient) & (np.abs(lat_c) >= EQUATORIAL_BAND_DEG)

    pressure = pressure_grid()
    level_columns = np.array([0, 1, 2])
    shear = np.zeros((3, lat_c.size))
    computed = np.full((3, lat_c.size), np.nan)
    u = np.repeat(u_ref[:, None], pressure.size, axis=1)
    u_levels = np.repeat(u_ref[None, :], 3, axis=0)
    changes = []
    for iteration in range(1, MAX_ITERATIONS + 1):
        for k in level_columns:
            u_k = u_levels[k]
            surface = gd.wind_geoid(phi, r_eq, "equatorial_radius",
                                    lambda x, u_k=u_k: np.interp(x, phi, u_k),
                                    Omega, GM, J, degrees, R_norm).radius
            r = np.asarray(surface, dtype="float64")
            g = g_eff_radial(u_k, r, phi, Omega, GM, J, degrees, R_norm)
            rate = omega_abs(u_k, r, phi, Omega)
            d_ln_p_d_r = -g * m_bar[k] / (MOLAR_GAS_CONSTANT * temperature[k])
            d = d_ln_p_d_r * (np.sin(phi) - np.cos(phi) / r * np.gradient(r, phi))
            d_ln_T = gradient[k] / temperature[k]
            with np.errstate(divide="ignore", invalid="ignore"):
                value = ((-g * d_ln_T / (2.0 * rate * r) - np.cos(phi) / r * np.gradient(u_k, phi))
                         / d)
            computed[k] = np.where(computed_mask[k], value, np.nan)
            shear[k] = shape_in_latitude(computed[k], lat_c)
        u_new = u_ref[:, None] + integrated_shear(shear, levels, pressure, above)
        change = float(np.max(np.abs(u_new - u)))
        changes.append(change)
        u = u_new
        u_levels = u_ref[None, :] + integrated_shear(shear, levels, levels, above).T
        if change < TOLERANCE_MS:
            break
    else:
        raise ValueError(f"the wind did not settle in {MAX_ITERATIONS} iterations: the last "
                         f"changes {changes[-3:]} m/s against a tolerance of {TOLERANCE_MS} m/s")
    dataset = _dataset(cloud, u, u_ref, pressure, levels, above)
    return Construction(dataset=dataset, level_Pa=levels, shear_ms=shear, shear_computed_ms=computed,
                        temperature_K=temperature, gradient_K_per_rad=gradient,
                        iterations=iteration, changes_ms=changes)


def _dataset(cloud, u, u_ref, pressure, levels, above) -> xr.Dataset:
    """The output kind W dataset, the cloud wind's attributes and auxiliaries replaced where false."""
    lat_dim, p_dim = "latitude_planetocentric", "pressure"
    reference_column = int(np.flatnonzero(pressure == REFERENCE_PRESSURE_Pa)[0])
    cloud_p = np.asarray(cloud["pressure_Pa"].values, dtype="float64")
    cloud_column = int(np.flatnonzero(cloud_p == float(cloud["reference_level_pressure_Pa"]))[0])
    sigma = np.asarray(cloud["u_total_uncertainty_ms"].values, dtype="float64")[:, cloud_column]
    flags = np.full(u.shape, PARAMETERIZED, dtype="int8")
    flags[:, reference_column] = np.asarray(cloud["value_provenance"].values)[:, cloud_column]
    keep = [n for n in cloud.data_vars if lat_dim not in cloud[n].dims or cloud[n].dims == (lat_dim,)]
    out = cloud[[n for n in keep if n not in ("u_reference_ms", "reference_level_pressure_Pa")]]
    out = out.drop_dims(p_dim, errors="ignore").copy()
    p_attrs = dict(cloud["pressure_Pa"].attrs)
    p_attrs["long_name"] = f"pressure grid, {NODES_PER_DECADE} nodes per decade"
    out = out.assign_coords(pressure_Pa=((p_dim,), pressure, p_attrs))
    two = (lat_dim, p_dim)
    top_mbar = levels[0] / 100.0
    out["u_total_ms"] = (two, u, dict(cloud["u_total_ms"].attrs))
    out["u_total_uncertainty_ms"] = (two, np.repeat(sigma[:, None], pressure.size, axis=1),
                                     dict(cloud["u_total_uncertainty_ms"].attrs))
    out["u_reference_ms"] = ((lat_dim,), u_ref, dict(cloud["u_reference_ms"].attrs))
    out["u_shear_ms"] = (two, u - u_ref[:, None], dict(cloud["u_shear_ms"].attrs))
    out["value_provenance"] = (two, flags, dict(cloud["value_provenance"].attrs))
    out["reference_level_pressure_Pa"] = ((), REFERENCE_PRESSURE_Pa, {
        "units": "Pa", "long_name": "level the cloud top wind is assigned to, where the IRIS shear "
        "is zero", "provenance": "assumed"})
    structure = (f"IRIS thermal wind shear (SPEC_11): the top level at {top_mbar:g} mbar, the shear "
                 f"{'held at its top level value' if above == 'held' else f'relaxed over {ABOVE_SCALE_HEIGHTS:g} scale heights'} "
                 f"above it and relaxed over {BELOW_SCALE_HEIGHTS:g} below 730 mbar")
    out.attrs.update({
        "title": f"{cloud.attrs.get('title', 'cloud top zonal wind')}, with the IRIS thermal wind "
                 "shear (SPEC_11)",
        "vertical_structure": structure,
        "method": f"{cloud.attrs.get('method', 'cloud tracking')}; the shear from Voyager IRIS "
                  "temperature gradients by the model's balance",
        "source": f"{cloud.attrs.get('source', '')}; shear from Conrath, B. J., and Pirraglia, J. A. "
                  "1983, Icarus 53, 286, Fig. 1, digitized (SPEC_09) and fitted (SPEC_10)",
        "observation_level_Pa": REFERENCE_PRESSURE_Pa,
        "observation_level_justification": "the cloud wind is assigned to 398 mbar, where the "
                                           "IRIS-derived shear is taken to be zero (SPEC_11 Step 1)",
        "coverage_pressure_Pa": np.array([pressure.min(), pressure.max()]),
    })
    return out
