"""Acceptance checks for SPEC_03 v0.12 Step 4: `forward.production`, `casspian-forward`, the closure
run, the product, F5 and F6; and the kind C `source_statement` rule (v0.11).

Every check prints its measured value; checks beyond the specification are labeled so. The anchor
and the run's inputs are the clean products of the Step 3 sweep (commit 6ae113e), so nothing is
relaxed. The product and its figures are written under `forward/lindal_closure/output/` (ignored by
git); the negative controls are computed in this script only, with no code path in `forward`.
"""

import math
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr
from matplotlib.image import imread

from casspian.forward import production as fp
from casspian.lib import composition as comp
from casspian.lib import control as ctl
from casspian.lib import geopotential as gp
from casspian.lib import hydrostatic as hs
from casspian.lib import io as cio
from casspian.lib.schema import CasspianSchemaError
from casspian.refrac.anchor import wind_of_latitude
from casspian.tools.plots import render, style

RUN = Path("forward/lindal_closure")
NAMELIST = RUN / "lindal_closure.toml"
PRODUCT = RUN / "output/lindal_closure_profile.nc"
FIGURES = RUN / "output/figures"
D = Path("occul_data/lindal")
HERE = Path("reports/step03_4")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def nodes(tree, base=""):
    out = {}
    root = tree.path.rstrip("/")
    for node in tree.subtree:
        relative = node.path[len(root):] if root else ("" if node.path == "/" else node.path)
        out[base + relative] = node.to_dataset(inherit=False)
    return out


HERE.mkdir(parents=True, exist_ok=True)
if (RUN / "output").exists():
    shutil.rmtree(RUN / "output")

# ---------------------------------------------------------------------------
# 1. casspian-forward writes the product; it reads back with the anchor and inputs embedded
# ---------------------------------------------------------------------------
done = subprocess.run(["casspian-forward", str(NAMELIST)], capture_output=True, text=True)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    tree = cio.read(PRODUCT, "profile")
    product = tree.load()
    tree.close()
hash_warnings = [str(w.message) for w in caught if "input_hashes" in str(w.message)]
embedded = nodes(product["anchors/lindal"])
on_disk_anchor = cio.read(D / "lindal_refractivity.nc", "refractivity")
anchor_disk = nodes(on_disk_anchor.load())
on_disk_anchor.close()
lines, ok = [], done.returncode == 0 and isinstance(product, xr.DataTree)
same = sorted(embedded) == sorted(anchor_disk) and all(embedded[g].identical(anchor_disk[g]) for g in anchor_disk)
ok = ok and same
lines.append(f"casspian-forward exit {done.returncode}: {done.stdout.strip()}{(' ' + done.stderr.strip()) if done.stderr.strip() else ''}")
lines.append(f"reads back as kind profile, a DataTree; anchors/lindal: {len(anchor_disk)} groups identical to "
             f"lindal_refractivity.nc on disk: {same}")
for key, kind in ctl.RUN_INPUT_KINDS.items():
    handle = cio.read(RUN / f"inputs/lindal_closure_{key}.nc", kind)
    loaded = handle.load()
    handle.close()
    disk = nodes(loaded) if isinstance(loaded, xr.DataTree) else {"": loaded}
    inside = nodes(product[f"inputs/{key}"])
    good = sorted(inside) == sorted(disk) and all(inside[g].identical(disk[g]) for g in disk)
    ok = ok and good
    lines.append(f"inputs/{key}: {len(disk)} group(s) identical to lindal_closure_{key}.nc on disk: {good}")
lines.append(f"input_hashes warnings on read: {hash_warnings or 'none'}")
record(1, "casspian-forward writes output/lindal_closure_profile.nc, which reads back as kind profile "
          "with the anchor and the four inputs embedded identical on every group", ok and not hash_warnings,
       "\n".join(lines))

root = product.to_dataset(inherit=False)
anchor_root = anchor_disk[""]
thermo = anchor_disk["/inputs/thermo"]
p_tab = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
T_tab = np.asarray(thermo["temperature_K"].values, dtype="float64")
p = np.asarray(root["pressure_Pa"].values, dtype="float64")
T = np.asarray(root["temperature_K"].values, dtype="float64")
Phi = np.asarray(root["geopotential_m2s2"].values, dtype="float64")
rec = dict(product["production_record"].attrs)

# ---------------------------------------------------------------------------
# 2. R_bar, m_bar and n from the run's kind C equal the anchor's to 1e-14
# ---------------------------------------------------------------------------
lines, ok = [], True
for name in ("mean_refractivity_m3", "mean_molar_mass_kg_mol", "number_density_m3"):
    a = np.asarray(root[name].values, dtype="float64")
    b = np.asarray(anchor_root[name].values, dtype="float64")
    worst = float(np.max(np.abs(a / b - 1.0)))
    ok = ok and worst <= 1e-14
    lines.append(f"{name}: largest relative difference from the anchor's {worst:.2e}")
record(2, "mean_refractivity_m3, mean_molar_mass_kg_mol and number_density_m3, formed from the run's "
          "kind C, equal the anchor's to 1e-14 relative", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 3. p at the top equals p_b exactly, and p_b is the anchor's top tabulated pressure; 4. Phi = 0 at gauge
# ---------------------------------------------------------------------------
p_b = float(root["boundary_pressure_Pa"].values)
record(3, "pressure_Pa at the top level equals p_b exactly, and p_b equals the anchor's top tabulated "
          "pressure", p[0] == p_b and p_b == p_tab[0] and float(root.attrs["boundary_pressure_Pa"]) == p_b,
       f"pressure_Pa[0] = {p[0]!r}; p_b = {p_b!r}; the anchor's top tabulated pressure {p_tab[0]!r}; the "
       f"global boundary_pressure_Pa {root.attrs['boundary_pressure_Pa']!r}; rule {rec['p_b_rule']!r}")
gauge = int(root["gauge_level_index"].values)
record(4, "Phi is zero at the gauge level", Phi[gauge] == 0.0 and p_tab[gauge] == 1.0e4,
       f"gauge level {gauge} at tabulated {p_tab[gauge]!r} Pa; Phi there {Phi[gauge]!r}; Phi from "
       f"{Phi[0]:.6e} to {Phi[-1]:.6e} m2/s2")

# ---------------------------------------------------------------------------
# 5. The residual bounds
# ---------------------------------------------------------------------------
res = p / p_tab - 1.0
bottom = int(np.argmax(p_tab))
not_bottom = np.arange(p.size) != bottom
# SPEC_03 v0.13 (REVIEW_03_step4 finding 1): the bins in tabulated pressure, a level on an edge in
# the deeper bin.
# The 2 mbar edge is the grid level Table I prints as 2.00 mbar (10**2.3 Pa); on the edge, deeper.
P2 = 10.0 ** 2.3
above_2 = p_tab < P2
from_2 = (p_tab >= P2) & not_bottom
from_10 = p_tab >= 1000.0
b_2_10 = (p_tab >= P2) & (p_tab < 1000.0)
b_10_100 = (p_tab >= 1000.0) & (p_tab < 10000.0)
b_100 = (p_tab >= 10000.0) & not_bottom
w = lambda m: float(np.max(np.abs(res[m])))
mean_10 = float(np.mean(res[from_10]))
ok = (w(above_2) <= 1.5e-3 and w(from_2) <= 3e-3 and abs(res[bottom]) <= 4e-3
      and abs(mean_10) <= 1.5e-3)
record(5, "the residual is within 1.5e-3 above 2 mbar, 3e-3 from 2 mbar down excluding the bottom row, "
          "4e-3 at the bottom row, and its mean below 10 mbar within 1.5e-3 (SPEC_03 v0.13)",
       ok,
       f"above 2 mbar ({int(above_2.sum())} levels): max |r| {w(above_2):.3e} (bound 1.5e-3)\n"
       f"from 2 mbar down, bottom row excluded: max |r| {w(from_2):.3e} (bound 3e-3)\n"
       f"bottom row ({p_tab[bottom] / 100:g} mbar): {res[bottom]:+.3e} (bound 4e-3)\n"
       f"mean below 10 mbar: {mean_10:+.3e} (bound 1.5e-3)\n"
       f"against the expected values: above 2 mbar {w(above_2):.2e} (9.8e-4); 2 to 10 mbar {w(b_2_10):.2e} "
       f"(2.2e-3); 10 to 100 mbar {w(b_10_100):.2e} (2.6e-3); below 100 mbar excluding the bottom row "
       f"{w(b_100):.2e} (1.7e-3); bottom row {res[bottom]:+.2e} (-3.3e-3); mean below 10 mbar {mean_10:+.2e} (-4.9e-4)\n"
       f"production_record agrees: above 2 mbar {rec['residual_max_abs_above_2mbar']:.3e}, 2 to 10 mbar "
       f"{rec['residual_max_abs_2_to_10mbar']:.3e}, from 2 mbar down {rec['residual_max_abs_from_2mbar_down_excluding_bottom_row']:.3e}, "
       f"bottom row {rec['residual_bottom_row']:+.3e}, mean below 10 mbar {rec['residual_mean_below_10mbar']:+.3e}")

# ---------------------------------------------------------------------------
# 6. The temperature residual equals the pressure residual
# ---------------------------------------------------------------------------
res_T = T / T_tab - 1.0
record(6, "the temperature residual equals the pressure residual to 1e-12 at every level",
       float(np.max(np.abs(res_T - res))) <= 1e-12,
       f"max |T/T_tab - 1 - (p/p_tab - 1)| = {float(np.max(np.abs(res_T - res))):.2e}")

# ---------------------------------------------------------------------------
# 7. The negative controls, in this script only
# ---------------------------------------------------------------------------
inputs = ctl.load_run_inputs(ctl.read_run_namelist(NAMELIST))
r = np.asarray(anchor_root["radius_m"].values, dtype="float64")
h = np.asarray(anchor_root["height_above_anchor_isobar_m"].values, dtype="float64")
N = np.asarray(anchor_root["refractivity"].values, dtype="float64")
phi_c = math.radians(float(anchor_root["latitude_planetocentric_deg"].values))
constants = (float(inputs.rotation["angular_rate_rad_s"]), float(inputs.gravity["GM_m3s2"]),
             np.asarray(inputs.gravity["J"].values), np.asarray(inputs.gravity["degree"].values),
             float(inputs.gravity["normalization_radius_m"]))
u = float(wind_of_latitude(inputs.wind)(np.array([phi_c]))[0])
_, g_radial, _, _ = gp.effective_gravity_magnitude(u, r, phi_c, *constants)
# SPEC_04 decision O: kind C is a field on (level, latitude) for every use, so `mean_properties`
# is given the latitude whose column to take, by decision L, exactly as `produce` gives it. The
# Lindal composition is uniform in latitude by construction, asserted here so that taking a
# column hides nothing.
_field = inputs.composition.to_dataset(inherit=False)
_varying = [name for name, v in _field.data_vars.items()
            if comp.LATITUDE_DIM in v.dims and not np.all(v.values == v.values[:, :1])]
assert not _varying, f"the composition varies with latitude in {_varying}; this control reads "\
                     "the column `produce` reads and assumes the Lindal field is uniform"
R_bar, m_bar = fp.mean_properties(inputs.composition, math.degrees(phi_c))
rho = hs.density(N, R_bar, m_bar)


def control_residual(increment):
    g_layer = 0.5 * (g_radial[:-1] + g_radial[1:])
    Phi_c = gp.geopotential(g_layer * increment, gauge)
    p_c = hs.pressure_from_top(p_tab[0], hs.layer_mass(rho, Phi_c))
    return p_c / p_tab - 1.0


one = control_residual(np.diff(h))
radial = control_residual(np.diff(r))
m_one, m_radial = float(np.mean(one[from_10])), float(np.mean(radial[from_10]))
record(7, "the negative control (radial component times the altitude increment) and the fully radial "
          "reading both have a mean below 10 mbar below -4e-3, failing the mean bound",
       m_one < -4e-3 and m_radial < -4e-3,
       f"one-factor g_k (h_(k+1) - h_k): mean below 10 mbar {m_one:+.3e} (about -5.1e-3), bottom row "
       f"{one[bottom]:+.3e} (about -7.8e-3)\n"
       f"fully radial g_k (r_(k+1) - r_k): mean below 10 mbar {m_radial:+.3e} (about -9.6e-3), bottom row "
       f"{radial[bottom]:+.3e} (about -1.2e-2)\n"
       f"the production for comparison: mean {mean_10:+.3e}, bottom row {res[bottom]:+.3e}")

# ---------------------------------------------------------------------------
# 8. F5 and F6 by the driver and by casspian-plots by hand, identical apart from the footer time
# ---------------------------------------------------------------------------
by_hand = HERE / "by_hand"
if by_hand.exists():
    shutil.rmtree(by_hand)
hand = subprocess.run(["casspian-plots", str(PRODUCT), "--out", str(by_hand)], capture_output=True, text=True)
names = ["lindal_closure_diag_F5_geopotential.png", "lindal_closure_diag_F6_hydrostatic.png"]
lines, ok = [], hand.returncode == 0
for name in names:
    a, b = imread(FIGURES / name), imread(by_hand / name)
    band = math.ceil(style.FOOTER_TIME_BAND * a.shape[0])
    body = a.shape == b.shape and bool(np.array_equal(a[:-band], b[:-band]))
    ok = ok and body
    lines.append(f"{name}: {a.shape[1]}x{a.shape[0]} px; identical above the {band} px time band: {body}; "
                 f"whole file byte identical: {(FIGURES / name).read_bytes() == (by_hand / name).read_bytes()}")
record(8, "F5 and F6 rendered by the driver and by casspian-plots by hand, identical apart from the footer "
          "time", ok, f"by hand: exit {hand.returncode}: {hand.stdout.strip()}\n" + "\n".join(lines))

# ---------------------------------------------------------------------------
# 9. F6: the residual against the envelope band, levels outside counted
# ---------------------------------------------------------------------------
drawn = render(PRODUCT, HERE / "render_data")
F6 = drawn.data["F6"]
envelope = np.asarray(F6["envelope"])
outside = np.flatnonzero(np.abs(res) > envelope)
record(9, "the residual is drawn in F6 against the envelope band; the levels outside the band are counted "
          "and reported with their values (a description, not a bound; v0.11)",
       len(drawn.written) == 2 and envelope.size == p.size,
       f"{outside.size} of {p.size} levels outside the band:\n"
       + "\n".join(f"  level {k:2d} {p_tab[k] / 100:9.4g} mbar residual {res[k]:+.3e} band {envelope[k]:.3e}"
                   for k in outside)
       + f"\nband from {envelope.min():.3e} to {envelope.max():.3e}")

# ---------------------------------------------------------------------------
# 10. The kind C source_statement rule (v0.11)
# ---------------------------------------------------------------------------
cases = HERE / "source_statement"
if cases.exists():
    shutil.rmtree(cases)
cases.mkdir(parents=True)
reduction_copy, forward_copy = cases / "lindal_composition.nc", cases / "lindal_closure_composition.nc"
shutil.copy2(D / "lindal_composition.nc", reduction_copy)
shutil.copy2(RUN / "inputs/lindal_closure_composition.nc", forward_copy)
for path in (reduction_copy, forward_copy):
    with netCDF4.Dataset(path, "a") as handle:
        handle.delncattr("source_statement")
lines, ok = [], True
try:
    cio.read(reduction_copy, "composition").close()
    ok = False
    lines.append("reduction composition without source_statement: not refused")
except CasspianSchemaError as exc:
    lines.append(f"reduction composition without source_statement: {exc}")
try:
    cio.read(forward_copy, "composition").close()
    lines.append("forward composition without source_statement: accepted")
except CasspianSchemaError as exc:
    ok = False
    lines.append(f"forward composition without source_statement: REFUSED {exc}")
record(10, "the kind C reader enforces source_statement for role reduction and not for forward (SPEC_03 "
           "v0.11 Step 4)", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 11. Beyond the specification: produce is pure and reproduces the product bit for bit
# ---------------------------------------------------------------------------
profile = fp.profile_from_anchor(inputs.anchor)
again = fp.produce(profile, inputs, inputs.gauge_level_index, float(p_tab[0]))
pairs = {"geopotential_m2s2": again.geopotential.geopotential_m2s2, "pressure_Pa": again.pressure_Pa,
         "temperature_K": again.temperature_K, "number_density_m3": again.number_density_m3,
         "mean_refractivity_m3": again.mean_refractivity_m3, "mean_molar_mass_kg_mol": again.mean_molar_mass_kg_mol}
equal = {k: bool(np.array_equal(np.asarray(root[k].values), v)) for k, v in pairs.items()}
record(11, "beyond the specification: produce() on the loaded arrays reproduces every modeled variable "
           "of the product bit for bit", all(equal.values()),
       f"{equal}; u(phi_c) {again.u_ms!r} m/s, record {rec['u_at_phi_c_ms']!r}")

# ---------------------------------------------------------------------------
# 12. Beyond the specification: the product's season, record, paths and companions
# ---------------------------------------------------------------------------
companions = [n for n in root.variables if "_uncertainty" in str(n)]
all_nan = all(np.all(np.isnan(np.asarray(root[n].values))) for n in companions)
paths = str(root.attrs["input_hashes"]).splitlines()
relative = all(not (line.startswith("/") or line[1:3] in (":/", ":\\")) for line in paths)
keys = ["mode", "gauge_level_index", "boundary_pressure_Pa", "p_b_rule", "geopotential_rule", "hydrostatic_rule",
        "closure_comparison_composition", "closure_comparison_wind", "anchor_radius_projection_rule",
        "residual_pressure", "residual_mean_below_10mbar", "casspian_git_commit"]
record(12, "beyond the specification: the product carries its season, the anchor's season, NaN companions, "
           "relative input_hashes and the production record",
       # SPEC_03 v0.13 (REVIEW_03_step4 finding 4): the closure namelist declares the occultation date.
       root.attrs["epoch"] == "1981-08-26"
       and float(root.attrs["solar_longitude_deg"]) == 18.2
       # A one-element array attribute reads back from netCDF as a scalar.
       and np.atleast_1d(root.attrs["anchor_solar_longitudes_deg"]).tolist() == [18.2]
       and all_nan and relative and all(k in rec for k in keys),
       f"epoch {root.attrs['epoch']!r}; solar_longitude_deg {root.attrs['solar_longitude_deg']!r}; "
       f"anchor_solar_longitudes_deg {np.asarray(root.attrs['anchor_solar_longitudes_deg']).tolist()}\n"
       f"{len(companions)} companions, all NaN: {all_nan}; terms unstated "
       f"{root['pressure_uncertainty_Pa'].attrs.get('uncertainty_terms_unstated')!r}\n"
       f"input_hashes:\n  " + "\n  ".join(paths) + "\n"
       f"production_record keys present: {[k for k in keys if k in rec]}\n"
       f"closure comparisons: {[rec[f'closure_comparison_{k}'] for k in ctl.RUN_INPUT_KINDS]}\n"
       f"namelist group resolved form:\n  " + str(product['namelist'].attrs['resolved']).replace('\n', '\n  '))

# ---------------------------------------------------------------------------
# 13. Beyond the specification: hashes
# ---------------------------------------------------------------------------
record(13, "beyond the specification: commit and SHA-256 of the product, provisional until the sweep", True,
       f"{PRODUCT.as_posix()} {cio.sha256(PRODUCT)} {root.attrs['casspian_git_commit']}")

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
