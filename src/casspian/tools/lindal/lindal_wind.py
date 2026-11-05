"""The Lindal-only wind: the cloud wind with the shear of the IRIS temperatures. SPEC_11 Step 1.

Inputs: the cloud wind (kind W, as `casspian-wind-from-curve` writes it), the digitized IRIS
temperatures (`occul_data/lindal/iris_temperatures.csv`, SPEC_09), and the run's composition,
gravity and rotation. Output: a kind W dataset with all three parts. `construct` touches no file;
a run's build writes the wind through the `lindal_iris` case of `casspian-wind-shear`, which calls
it (SPEC_11 v0.4 Step 2).

**The reference wind** is the cloud wind's reference wind, assigned to `REFERENCE_PRESSURE_Pa`
(398 mbar, a node of the output grid), where the shear is zero.

**The gradient** at each of the three levels is from `iris_temperatures.fit` at its defaults, evaluated
at the planetographic latitudes of the cloud wind's planetocentric grid (the file's own
`latitude_planetographic_deg`, which the wind tool computed by its rule), and converted to per
planetocentric radian with `dphi_g/dphi_c` differenced on that grid. `gradient` chooses which
(SPEC_13 Step 1 item 1): `"slope"`, the fit's local slope `gradient_K_per_deg` (SPEC_11, the
default); or `"value"`, the derivative of the fit's `value_K` in planetographic latitude, by
central differences on that grid and one-sided at the ends of the fit's span (`np.gradient`), so
that the shear stands for the temperatures the fit gives. `window` is `fit`'s: SPEC_10's adaptive
width (the default) or the smooth one (SPEC_13 v0.5 section 2a item 1).

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

**In the vertical**, between the three levels linear in `ln p`; below the bottom level relaxed as
`s_bottom exp(-y / BELOW_SCALE_HEIGHTS)`, `y = ln(p / p_bottom)`. Above the top level, by case:
held at its value (i); relaxed as `s_top exp(-x / ABOVE_SCALE_HEIGHTS)`, `x = ln(p_top / p)` (ii);
or reversed (iii, SPEC_14): `s_top f(p)` with `f = 1` from the top level up to `hold_Pa`, then
`f = 1 - X / X_c`, `X = ln(hold_Pa / p)`, `X_c = ln(hold_Pa / cross_Pa)`, through zero at
`cross_Pa`, up to `REVERSED_TOP_Pa` (0.2 mbar, the anchor's top), and held at its value there above.
The wind is the reference wind plus the shear integrated in `ln p` from the reference level, in
closed form.

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
ABOVE_CASES = ("held", "relaxed", "reversed")
GRADIENTS = ("slope", "value")
ABOVE_SCALE_HEIGHTS = 2.0
#: Case (iii): the level above which f is held at its value (SPEC_14 section 2 ruling 1).
REVERSED_TOP_Pa = 20.0
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


def reversed_factor(pressure_Pa, hold_Pa, cross_Pa) -> np.ndarray:
    """Case (iii)'s `f(p)`: 1 at and below `hold_Pa`, linear in `ln p` through zero at `cross_Pa`,
    held at its `REVERSED_TOP_Pa` value above that (module docstring)."""
    x = np.log(hold_Pa / np.maximum(np.asarray(pressure_Pa, dtype="float64"), REVERSED_TOP_Pa))
    return 1.0 - np.maximum(x, 0.0) / np.log(hold_Pa / cross_Pa)


def integrated_shear(shear, levels_Pa, pressure_Pa, above, hold_Pa=None, cross_Pa=None):
    """`int s dln p` from the reference level, closed form, shape `(latitude, pressure)`.

    `shear` is `(3, latitude)` at `levels_Pa` (top first). Linear in `ln p` between the levels,
    held, relaxed or reversed above the top, relaxed below the bottom (module docstring)."""
    x1, x2, x3 = np.log(levels_Pa)
    s1, s2, s3 = (shear[i][:, None] for i in range(3))

    def antiderivative(x):
        x = np.asarray(x, dtype="float64")[None, :]
        f2 = 0.5 * (s1 + s2) * (x2 - x1)
        f3 = f2 + 0.5 * (s2 + s3) * (x3 - x2)
        if above == "held":
            top = s1 * (x - x1)
        elif above == "relaxed":
            top = -s1 * ABOVE_SCALE_HEIGHTS * (1.0 - np.exp(-(x1 - x) / ABOVE_SCALE_HEIGHTS))
        else:
            # f = 1 down to x_h = ln hold_Pa, then 1 - (x_h - x) / x_c down to x_t, then constant.
            x_h, x_t = np.log(hold_Pa), np.log(REVERSED_TOP_Pa)
            x_c = np.log(hold_Pa / cross_Pa)
            depth = x_h - np.maximum(x, x_t)                   # X, at most X_t = x_h - x_t
            ramp = s1 * (x_h - x1) - s1 * (depth - depth ** 2 / (2.0 * x_c))
            f_t = 1.0 - (x_h - x_t) / x_c
            beyond = ramp + s1 * f_t * np.minimum(x - x_t, 0.0)
            top = np.where(x >= x_h, s1 * (x - x1), beyond)
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


def derivative_of_value(value_K, latitude_g_deg) -> np.ndarray:
    """`dT/dphi_g` (K per degree) of a fit's values on its grid, NaN where the value is NaN.

    The fit is defined on one contiguous span of the grid (between the data's ends); the
    derivative is central inside it and one-sided at its two ends."""
    value_K = np.asarray(value_K, dtype="float64")
    latitude_g_deg = np.asarray(latitude_g_deg, dtype="float64")
    out = np.full(value_K.shape, np.nan)
    inside = np.flatnonzero(np.isfinite(value_K))
    if inside.size < 2:
        return out
    if inside[-1] - inside[0] + 1 != inside.size:
        raise ValueError("the fit's values are not one contiguous span of the grid")
    out[inside] = np.gradient(value_K[inside], latitude_g_deg[inside])
    return out


def construct(cloud, table, composition, gravity, rotation, top_level_Pa, above,
              equatorial_radius_m=None, gradient="slope", window="adaptive",
              hold_Pa=None, cross_Pa=None) -> Construction:
    """The wind for one case, from in-memory inputs. No file is read or written.

    `table` is the IRIS CSV as a structured array (`numpy.genfromtxt(..., names=True)`).
    `hold_Pa` and `cross_Pa` are case (iii)'s, required by it and refused by the others."""
    if above not in ABOVE_CASES:
        raise ControlFileError(f"above_top = {above!r}; the cases are {list(ABOVE_CASES)}")
    reversal = (hold_Pa, cross_Pa)
    if above == "reversed":
        if None in reversal:
            raise ControlFileError("above_top = 'reversed' requires hold_Pa and cross_Pa")
        hold_Pa, cross_Pa = float(hold_Pa), float(cross_Pa)
        if not float(top_level_Pa) >= hold_Pa > cross_Pa > REVERSED_TOP_Pa:
            raise ControlFileError(
                f"above_top = 'reversed' needs top_level_Pa >= hold_Pa > cross_Pa > "
                f"{REVERSED_TOP_Pa:g} Pa; given {float(top_level_Pa):g}, {hold_Pa:g} and {cross_Pa:g}")
    elif reversal != (None, None):
        raise ControlFileError(f"hold_Pa and cross_Pa are case (iii)'s; above_top = {above!r} "
                               "does not use them")
    if gradient not in GRADIENTS:
        raise ControlFileError(f"gradient = {gradient!r}; the choices are {list(GRADIENTS)}")
    if window not in iris.WINDOWS:
        raise ControlFileError(f"window = {window!r}; the windows are {list(iris.WINDOWS)}")
    gradient_choice = gradient
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
                          grid=lat_g, window=window)
        temperature[k] = result.value_K
        per_deg = (result.gradient_K_per_deg if gradient_choice == "slope"
                   else derivative_of_value(result.value_K, lat_g))
        gradient[k] = per_deg * dg_dc * (180.0 / np.pi)
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
        u_new = u_ref[:, None] + integrated_shear(shear, levels, pressure, above, hold_Pa, cross_Pa)
        change = float(np.max(np.abs(u_new - u)))
        changes.append(change)
        u = u_new
        u_levels = u_ref[None, :] + integrated_shear(shear, levels, levels, above, hold_Pa, cross_Pa).T
        if change < TOLERANCE_MS:
            break
    else:
        raise ValueError(f"the wind did not settle in {MAX_ITERATIONS} iterations: the last "
                         f"changes {changes[-3:]} m/s against a tolerance of {TOLERANCE_MS} m/s")
    dataset = _dataset(cloud, u, u_ref, pressure, levels, above, gradient_choice, window,
                       hold_Pa, cross_Pa)
    return Construction(dataset=dataset, level_Pa=levels, shear_ms=shear, shear_computed_ms=computed,
                        temperature_K=temperature, gradient_K_per_rad=gradient,
                        iterations=iteration, changes_ms=changes)


def balance_weight(wind, composition, gravity, rotation, pressure_Pa, temperature_K,
                   equatorial_radius_m=None) -> np.ndarray:
    """The balance `construct` solves, linearized in the shear with the wind held fixed (SPEC_17
    Step 1 item 3): on each isobar `Delta(d ln T / dphi)_p = w Delta s`, with

        w = -2 Omega_abs r D / g = 2 Omega_abs r m_bar (sin(phi) - (cos(phi) / r)(dr/dphi)_p) / (R T)

    and `D` as in the module docstring. The geometry at each pressure is `construct`'s: the Eq. B3
    surface under that level's wind, `g` and `Omega_abs` under it. Held fixed, the wind's own terms
    (`(du/dphi)_p`, and the wind in `g`, `Omega_abs` and the surface) do not enter `w`.

    `wind` is a kind W dataset, read at `pressure_Pa` linearly in `ln p`; `temperature_K` is
    `(pressure, latitude)` on the wind's planetocentric grid, NaN where unknown. Returns `w` in
    s m-1 rad-1, `(pressure, latitude)`, which changes sign at the equator with `D`."""
    lat_c = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    phi = np.radians(lat_c)
    grid = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    u_grid = np.asarray(wind["u_total_ms"].values, dtype="float64")
    pressure_Pa = np.atleast_1d(np.asarray(pressure_Pa, dtype="float64"))
    temperature_K = np.asarray(temperature_K, dtype="float64").reshape(pressure_Pa.size, lat_c.size)
    rising = np.argsort(np.log(grid))
    Omega = float(rotation["angular_rate_rad_s"])
    GM = float(gravity["GM_m3s2"])
    J = np.asarray(gravity["J"].values, dtype="float64")
    degrees = np.asarray(gravity["degree"].values)
    R_norm = float(gravity["normalization_radius_m"])
    r_eq = R_norm if equatorial_radius_m is None else float(equatorial_radius_m)
    m_bar = mean_molar_mass(composition, pressure_Pa, lat_c)
    w = np.full((pressure_Pa.size, lat_c.size), np.nan)
    for k, level in enumerate(pressure_Pa):
        u_k = np.array([np.interp(np.log(level), np.log(grid[rising]), row[rising]) for row in u_grid])
        r = np.asarray(gd.wind_geoid(phi, r_eq, "equatorial_radius",
                                     lambda x, u_k=u_k: np.interp(x, phi, u_k),
                                     Omega, GM, J, degrees, R_norm).radius, dtype="float64")
        rate = omega_abs(u_k, r, phi, Omega)
        w[k] = (2.0 * rate * r * m_bar[k] * (np.sin(phi) - np.cos(phi) / r * np.gradient(r, phi))
                / (MOLAR_GAS_CONSTANT * temperature_K[k]))
    return w


def _dataset(cloud, u, u_ref, pressure, levels, above, gradient="slope", window="adaptive",
             hold_Pa=None, cross_Pa=None) -> xr.Dataset:
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
    above_text = {
        "held": "held at its top level value",
        "relaxed": f"relaxed over {ABOVE_SCALE_HEIGHTS:g} scale heights",
        "reversed": (f"reversed (case iii, SPEC_14): held to {(hold_Pa or 0) / 100:g} mbar, then "
                     f"linear in ln p through zero at {(cross_Pa or 0) / 100:g} mbar, held above "
                     f"{REVERSED_TOP_Pa / 100:g} mbar; hold_Pa and cross_Pa tuned to the Voyager 2 "
                     "egress, which then checks consistency, not independence"),
    }[above]
    structure = (f"IRIS thermal wind shear (SPEC_11): the top level at {top_mbar:g} mbar, the shear "
                 f"{above_text} above it and relaxed over {BELOW_SCALE_HEIGHTS:g} below 730 mbar")
    out.attrs.update({
        "title": f"{cloud.attrs.get('title', 'cloud top zonal wind')}, with the IRIS thermal wind "
                 "shear (SPEC_11)",
        "vertical_structure": structure,
        "method": f"{cloud.attrs.get('method', 'cloud tracking')}; the shear from Voyager IRIS "
                  "temperature gradients by the model's balance"
                  + ("" if gradient == "slope" else
                     ", the gradient the derivative of the fitted temperatures (SPEC_13)")
                  + ("" if window == "adaptive" else ", the fit's window smooth in latitude (SPEC_13)"),
        "source": f"{cloud.attrs.get('source', '')}; shear from Conrath, B. J., and Pirraglia, J. A. "
                  "1983, Icarus 53, 286, Fig. 1, digitized (SPEC_09) and fitted (SPEC_10)",
        "observation_level_Pa": REFERENCE_PRESSURE_Pa,
        "coverage_pressure_Pa": np.array([pressure.min(), pressure.max()]),
    })
    return out
