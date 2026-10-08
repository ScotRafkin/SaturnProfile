"""Diagnosis for SPEC_11 Step 2: where the round trip's differences come from. Report only.

The round trip compares the delivered temperature change from the anchor to the target with the
IRIS fit's change between the same latitudes. This splits the latitude path into the equatorial
band, where the shear is the PCHIP bridge and not the data's, and the rest, where the tool
inverted the balance from the data. On each piece it integrates two gradients of `ln T` along the
level:

* the IRIS fit's own (`iris_temperatures.fit` at its defaults), the data;
* the gradient the run's wind implies by the same balance the tool inverted, the wind read from the
  run's own file (`inputs/lindal_iris_ii_110_wind.nc`), `u` and `du/dln p` at the level by its
  interpolant, on the tool's geometry (`lindal_wind`'s: the Eq. B3 surface under the level's wind).

It also gives the anchor's own lapse rate against the same dry adiabat as the acceptance, to tell a
property of the transfer from a property of the adiabat's heat capacity.

Writes `reports/step11_2/diagnose_round_trip.txt`. Run from the repository root after the
acceptance.
"""

import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import geoid as gd
from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.lib.gravity import g_eff_radial, omega_abs
from casspian.lib.windfield import WindField
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.lindal import lindal_wind as lw

warnings.simplefilter("ignore")
HERE = Path("reports/step11_2")
RUN = Path("forward/lindal_iris_ii_110")
LEVELS = {110: 11000.0, 290: 29000.0, 730: 73000.0}


def load(path):
    handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


product = load(RUN / "output/lindal_iris_ii_110_profile.nc")
wind = product["inputs/wind"].to_dataset(inherit=False)
gravity = product["inputs/gravity"].to_dataset(inherit=False)
rotation = product["inputs/rotation"].to_dataset(inherit=False)
composition = product["inputs/composition"]
anchor_root = product["anchors/lindal"].to_dataset(inherit=False)
thermo = product["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)
root = product.to_dataset(inherit=False)
table = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)

lat_c = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
lat_g = np.asarray(wind["latitude_planetographic_deg"].values, dtype="float64")
phi = np.radians(lat_c)
dg_dc = np.gradient(lat_g, lat_c)
Omega = float(rotation["angular_rate_rad_s"])
GM, R_norm = float(gravity["GM_m3s2"]), float(gravity["normalization_radius_m"])
J, degrees = np.asarray(gravity["J"].values), np.asarray(gravity["degree"].values)
field = WindField(wind)
anchor_c = float(anchor_root["latitude_planetocentric_deg"])
target_c = float(root["latitude_planetocentric_deg"])
lines = [f"path {anchor_c:.3f} to {target_c:.2f} deg planetocentric; the band |phi_c| < "
         f"{lw.EQUATORIAL_BAND_DEG:g} deg"]

for level, pressure in LEVELS.items():
    rows = table["pressure_mbar"] == level
    fit = iris.fit(table["latitude_planetographic_deg"][rows], table["temperature_K"][rows], grid=lat_g)
    T = fit.value_K
    d_ln_T_iris = fit.gradient_K_per_deg * dg_dc * (180.0 / np.pi) / T
    u = field.wind_at(phi, np.full(phi.shape, pressure))
    _, s = field.wind_derivatives(phi, np.full(phi.shape, pressure))
    r = np.asarray(gd.wind_geoid(phi, R_norm, "equatorial_radius",
                                 lambda x: np.interp(x, phi, u), Omega, GM, J, degrees, R_norm).radius)
    g = g_eff_radial(u, r, phi, Omega, GM, J, degrees, R_norm)
    m_bar = lw.mean_molar_mass(composition, np.array([pressure]), lat_c)[0]
    d = -g * m_bar / (MOLAR_GAS_CONSTANT * T) * (np.sin(phi) - np.cos(phi) / r * np.gradient(r, phi))
    d_u_d_Z = s * d + np.cos(phi) / r * np.gradient(u, phi)
    d_ln_T_wind = -2.0 * omega_abs(u, r, phi, Omega) * r * d_u_d_Z / g

    on_path = (lat_c <= anchor_c) & (lat_c >= target_c)
    band = on_path & (np.abs(lat_c) < lw.EQUATORIAL_BAND_DEG)
    outside = on_path & ~band

    def integral(values, mask):
        x, y = phi[mask], values[mask]
        return float(np.trapezoid(y, x)) if x.size > 1 else 0.0

    T_anchor = float(np.interp(anchor_c, lat_c, T))
    parts = {}
    for name, mask in (("outside the band", outside), ("inside the band", band)):
        north = mask & (lat_c >= 0)
        south = mask & (lat_c <= 0)
        iris_part = -(integral(d_ln_T_iris, north) + integral(d_ln_T_iris, south))
        wind_part = -(integral(d_ln_T_wind, north) + integral(d_ln_T_wind, south))
        parts[name] = (T_anchor * iris_part, T_anchor * wind_part)
    lines.append(
        f"{level} mbar, the change from the anchor to the target in K, the IRIS fit's against the "
        f"run's wind's implied, as T_anchor x the integral of d ln T:")
    for name, (a, b) in parts.items():
        lines.append(f"    {name:17s}: IRIS {a:+.2f}, implied by the wind {b:+.2f}, difference {b - a:+.2f}")
    total_iris = sum(a for a, _ in parts.values())
    total_wind = sum(b for _, b in parts.values())
    lines.append(f"    {'whole path':17s}: IRIS {total_iris:+.2f}, implied by the wind {total_wind:+.2f}, "
                 f"difference {total_wind - total_iris:+.2f}")

# The anchor's own lapse rate against the acceptance's dry adiabat
p_a = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
T_a = np.asarray(thermo["temperature_K"].values, dtype="float64")
z_a = np.asarray(thermo["height_m"].values, dtype="float64")
comp = composition.to_dataset(inherit=False)
column = int(np.argmin(np.abs(np.asarray(comp["latitude_planetocentric_deg"].values) - anchor_c)))
p_comp = np.asarray(comp["pressure_Pa"].values, dtype="float64")
x = {n: np.interp(np.log(p_a), np.log(p_comp), np.asarray(comp[f"x_{n}"].values)[:, column])
     for n in ("H2", "He", "NH3")}
m_bar = np.asarray(anchor_root["mean_molar_mass_kg_mol"].values, dtype="float64")
c_p = MOLAR_GAS_CONSTANT * (3.5 * x["H2"] + 2.5 * x["He"] + 4.0 * x["NH3"]) / m_bar
order = np.argsort(z_a)
lapse = -np.gradient(T_a[order], z_a[order])
p_sorted = p_a[order]
g_a = 10.4
deep = p_sorted > 6.0e4
lines.append(
    f"the anchor's own lapse rate deeper than 600 mbar: {np.round(lapse[deep] * 1e3, 3).tolist()} K/km "
    f"at {np.round(p_sorted[deep] / 100.0, 0).tolist()} mbar; the dry adiabat with frozen rotation "
    f"g / c_p about {g_a / float(np.median(c_p[order][deep])) * 1e3:.3f} K/km (g about {g_a} m/s2)")

text = "\n".join(lines)
print(text)
(HERE / "diagnose_round_trip.txt").write_text(text + "\n", encoding="utf-8")
